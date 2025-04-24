version 1.0

struct SampleInfo {
    String omics_sample_name
    File R1
    File R2
}

struct RefGenome {
    String name
    File fasta
    File gtf
}

workflow STAR2Pass {
  input {
    Array[SampleInfo] samples
    RefGenome? reference_genome
  }

  if (!defined(reference_genome)) {
    call DownloadReference {}
  }

  RefGenome ref_genome_final = select_first([reference_genome, DownloadReference.genome])

  call CollapseGTF {
    input:
      reference_gtf = ref_genome_final.gtf
  }

  call BuildSTARIndex {
    input:
      reference_fasta = ref_genome_final.fasta,
      reference_gtf = ref_genome_final.gtf
  }

  scatter (sample in samples) {
    call STARalignTwoPass {
      input:
        base_file_name = sample.omics_sample_name,
        star_genome_tar = BuildSTARIndex.star_index_tar,
        r1fastq = sample.R1,
        r2fastq = sample.R2,
        ref_genome_name = ref_genome_final.name
    }

    call RNASeQC {
      input:
        base_file_name = sample.omics_sample_name,
        bam_file = STARalignTwoPass.bam,
        bam_index = STARalignTwoPass.bai,
        ref_gtf = CollapseGTF.collapsed_gtf
    }
  }

  call CombineCountMatricesPython {
    input:
      gene_count_files = STARalignTwoPass.geneCounts,
      sample_names = samples.omics_sample_name
  }

  output {
    Array[File] output_bam = STARalignTwoPass.bam
    Array[File] output_bai = STARalignTwoPass.bai
    Array[File] output_geneCounts = STARalignTwoPass.geneCounts
    Array[File] output_log_final = STARalignTwoPass.log_final
    Array[File] output_log_progress = STARalignTwoPass.log_progress
    Array[File] output_log = STARalignTwoPass.log
    Array[File] output_SJ = STARalignTwoPass.SJout
    Array[File] output_rnaseqc = RNASeQC.rnaseqc_metrics
    File combined_counts_matrix = CombineCountMatricesPython.counts_matrix
    File sample_metadata_template = CombineCountMatricesPython.sample_metadata
  }
}

task DownloadReference {
  input {}

  command <<<
    wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
    wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
    gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
    gunzip GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
  >>>

  output {
    RefGenome genome = object {
      name: "hg38",
      fasta: "GCF_000001405.40_GRCh38.p14_genomic.fna",
      gtf: "GCF_000001405.40_GRCh38.p14_genomic.gtf"
    }
  }

  runtime {
    docker: "getwilds/gtf-smash:latest"
    memory: "4 GB"
    cpu: "1"
  }
}

task BuildSTARIndex {
  input {
    File reference_fasta
    File reference_gtf
    Int sjdbOverhang = 100
    Int genomeSAindexNbases = 14
    Int memory_gb = 64
    Int cpu_cores = 8
  }

  command <<<
    set -eo pipefail
    
    mkdir star_index

    echo "Building STAR index..."
    STAR \
      --runMode genomeGenerate \
      --runThreadN ~{cpu_cores} \
      --genomeDir star_index \
      --genomeFastaFiles ~{reference_fasta} \
      --sjdbGTFfile ~{reference_gtf} \
      --sjdbOverhang ~{sjdbOverhang} \
      --genomeSAindexNbases ~{genomeSAindexNbases}

    tar -czf star_index.tar.gz star_index/
  >>>

  output {
    File star_index_tar = "star_index.tar.gz"
  }

  runtime {
    docker: "getwilds/star:2.7.6a"
    memory: "~{memory_gb} GB"
    cpu: "~{cpu_cores}"
  }
}

task CollapseGTF {
  input {
    File reference_gtf
  }

  command <<<
    set -eo pipefail
    
    echo "Processing GTF file..."
    collapse_annotation.py \
      ~{reference_gtf} \
      collapsed.gtf
  >>>

  output {
    File collapsed_gtf = "collapsed.gtf"
  }

  runtime {
    docker: "getwilds/gtf-smash:latest"
    memory: "4 GB"
    cpu: "1"
  }
}

task STARalignTwoPass {
  input {
    File star_genome_tar
    File r1fastq
    File r2fastq
    String base_file_name
    String ref_genome_name
    Int memory_gb = 62
    Int cpu_cores = 8
    Int star_threads = 6
  }

  command <<<
    set -eo pipefail

    echo "Extracting STAR reference..."
    tar -xvf "~{star_genome_tar}"

    echo "Starting STAR alignment..."
    STAR \
      --genomeDir star_index \
      --readFilesIn "~{r1fastq}" "~{r2fastq}" \
      --runThreadN ~{star_threads} \
      --readFilesCommand zcat \
      --sjdbOverhang 100 \
      --outSAMtype BAM SortedByCoordinate \
      --twopassMode Basic \
      --outTmpDir _STARtmp \
      --outFileNamePrefix "./" \
      --quantMode GeneCounts \
      --quantTranscriptomeBAMcompression 5 

    rm -r star_index _STARtmp

    mv Aligned.sortedByCoord.out.bam "~{base_file_name}.~{ref_genome_name}.Aligned.sortedByCoord.out.bam"
    mv ReadsPerGene.out.tab "~{base_file_name}.~{ref_genome_name}.ReadsPerGene.out.tab"
    mv Log.final.out "~{base_file_name}.~{ref_genome_name}.Log.final.out"
    mv Log.progress.out "~{base_file_name}.~{ref_genome_name}.Log.progress.out"
    mv Log.out "~{base_file_name}.~{ref_genome_name}.Log.out"
    mv SJ.out.tab "~{base_file_name}.~{ref_genome_name}.SJ.out.tab"

    samtools index "~{base_file_name}.~{ref_genome_name}.Aligned.sortedByCoord.out.bam"
  >>>

  output {
    File bam = "~{base_file_name}.~{ref_genome_name}.Aligned.sortedByCoord.out.bam"
    File bai = "~{base_file_name}.~{ref_genome_name}.Aligned.sortedByCoord.out.bam.bai"
    File geneCounts = "~{base_file_name}.~{ref_genome_name}.ReadsPerGene.out.tab"
    File log_final = "~{base_file_name}.~{ref_genome_name}.Log.final.out"
    File log_progress = "~{base_file_name}.~{ref_genome_name}.Log.progress.out"
    File log = "~{base_file_name}.~{ref_genome_name}.Log.out"
    File SJout = "~{base_file_name}.~{ref_genome_name}.SJ.out.tab"
  }

  runtime {
    docker: "getwilds/star:2.7.6a"
    memory: "~{memory_gb} GB"
    cpu: "~{cpu_cores}"
  }
}

task RNASeQC {
  input {
    File bam_file
    File bam_index
    File ref_gtf
    String base_file_name
    Int memory_gb = 4
    Int cpu_cores = 2
  }

  command <<<
    echo "Running RNA-SeQC..."
    rnaseqc "~{ref_gtf}" "~{bam_file}" OUTPUT \
      --sample="~{base_file_name}" \
      --coverage 

    tar -cvzf "~{base_file_name}.QC.tar.gz" OUTPUT/*
  >>>

  output {
    File rnaseqc_metrics = "~{base_file_name}.QC.tar.gz"
  }

  runtime {
    docker: "getwilds/rnaseqc:2.4.2"
    memory: "~{memory_gb} GB"
    cpu: "~{cpu_cores}"
  }
}

task CombineCountMatricesPython {
  input {
    Array[File] gene_count_files
    Array[String] sample_names
    Int memory_gb = 4
    Int cpu_cores = 1
    # Column to extract from ReadsPerGene.out.tab files:
    # 2 = unstranded counts
    # 3 = stranded counts, first read forward
    # 4 = stranded counts, first read reverse
    Int count_column = 2
  }

  command <<<
    #!/usr/bin/env python3
    
    import pandas as pd
    import os
    import sys
    
    # Get the list of gene count files and sample names
    count_files = ["~{sep='","' gene_count_files}"]
    sample_names = ["~{sep='","' sample_names}"]
    count_column = ~{count_column}
    
    print(f"Processing {len(count_files)} count files...")
    
    # Function to read STAR gene count file
    def read_star_counts(file_path, sample_name, count_col):
        # Skip the first 4 lines (summary statistics)
        df = pd.read_csv(file_path, sep='\t', skiprows=4, header=None)
        
        # Select only gene ID column and the requested count column
        df = df.iloc[:, [0, count_col-1]]
        
        # Name the columns
        df.columns = ['gene_id', sample_name]
        
        return df
    
    # Read the first file to get the gene list
    print(f"Reading first file: {os.path.basename(count_files[0])}")
    combined = read_star_counts(count_files[0], sample_names[0], count_column)
    
    # Add the rest of the samples
    if len(count_files) > 1:
        for i in range(1, len(count_files)):
            print(f"Reading file {i+1}/{len(count_files)}: {os.path.basename(count_files[i])}")
            sample_counts = read_star_counts(count_files[i], sample_names[i], count_column)
            combined = pd.merge(combined, sample_counts, on='gene_id')
    
    # Write out the combined matrix
    print(f"Writing combined matrix to {output_name}...")
    combined.to_csv("combined_counts_matrix.txt", sep='\t', index=False)
    
    # Create a sample metadata template for DESeq2
    metadata = pd.DataFrame({
        'sample_id': sample_names,
        'condition': ['condition'] * len(sample_names)
    })
    metadata.to_csv("sample_metadata_template.txt", sep='\t', index=False)
    
    # Print summary
    print(f"Combined {len(sample_names)} samples into a single count matrix")
    print(f"Total genes: {len(combined)}")
    print("Output files:")
    print(f"  - combined_counts_matrix.txt (main counts matrix)")
    print("  - sample_metadata_template.txt (template for DESeq2 sample metadata)")
  >>>

  output {
    File counts_matrix = "~{output_name}"
    File sample_metadata = "sample_metadata_template.txt"
  }

  runtime {
    docker: "python:3.9-slim"
    memory: "~{memory_gb} GB"
    cpu: "~{cpu_cores}"
  }
}

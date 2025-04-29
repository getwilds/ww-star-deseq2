version 1.0

struct SampleInfo {
    String omics_sample_name
    File R1
    File R2
    String condition
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
    String reference_level = ""
    String contrast = ""
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

  call CombineCountMatrices {
    input:
      gene_count_files = STARalignTwoPass.geneCounts,
      samples = samples
  }

  call RunDESeq2 {
    input:
      counts_matrix = CombineCountMatrices.counts_matrix,
      sample_metadata = CombineCountMatrices.sample_metadata,
      reference_level = reference_level,
      contrast = contrast
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
    File combined_counts_matrix = CombineCountMatrices.counts_matrix
    File sample_metadata_template = CombineCountMatrices.sample_metadata
    File deseq2_all_results = RunDESeq2.deseq2_results
    File deseq2_significant_results = RunDESeq2.deseq2_significant
    File deseq2_normalized_counts = RunDESeq2.deseq2_normalized_counts
    File deseq2_pca_plot = RunDESeq2.deseq2_pca_plot
    File deseq2_volcano_plot = RunDESeq2.deseq2_volcano_plot
    File deseq2_heatmap = RunDESeq2.deseq2_heatmap
  }
}

task DownloadReference {
  input {}

  command <<<
    set -eo pipefail
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
    set -eo pipefail

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

task CombineCountMatrices {
  input {
    Array[File] gene_count_files
    Array[SampleInfo] samples
    Int memory_gb = 4
    Int cpu_cores = 1
    # Column to extract from ReadsPerGene.out.tab files:
    # 2 = unstranded counts
    # 3 = stranded counts, first read forward
    # 4 = stranded counts, first read reverse
    Int count_column = 2
  }

  command <<<
    set -eo pipefail

    combine_star_counts.py \
      --input ~{sep=' ' gene_count_files} \
      --output combined_counts_matrix.txt \
      --metadata sample_metadata.txt \
      --count_column ~{count_column} \
      --samples "~{write_json(samples)}"
  >>>

  output {
    File counts_matrix = "combined_counts_matrix.txt"
    File sample_metadata = "sample_metadata.txt"
  }

  runtime {
    docker: "ghcr.io/tefirman/combine-counts:latest"
    memory: "~{memory_gb} GB"
    cpu: "~{cpu_cores}"
  }
}

task RunDESeq2 {
  input {
    File counts_matrix
    File sample_metadata
    String condition_column = "condition"
    String reference_level = ""
    String contrast = ""
    Int memory_gb = 8
    Int cpu_cores = 2
  }

  command <<<
    set -eo pipefail
    
    Rscript /deseq2_analysis.R \
      --counts_file="~{counts_matrix}" \
      --metadata_file="~{sample_metadata}" \
      --condition_column="~{condition_column}" \
      --reference_level="~{reference_level}" \
      --contrast="~{contrast}" \
      --output_prefix="deseq2_results"
  >>>

  output {
    File deseq2_results = "deseq2_results_all_genes.csv"
    File deseq2_significant = "deseq2_results_significant.csv"
    File deseq2_normalized_counts = "deseq2_normalized_counts.csv"
    File deseq2_pca_plot = "deseq2_results_pca.pdf"
    File deseq2_volcano_plot = "deseq2_results_volcano.pdf"
    File deseq2_heatmap = "deseq2_results_heatmap.pdf"
  }

  runtime {
    docker: "ghcr.io/tefirman/deseq2:latest"
    memory: "~{memory_gb} GB"
    cpu: "~{cpu_cores}"
  }
}

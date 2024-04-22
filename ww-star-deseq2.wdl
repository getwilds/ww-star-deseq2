version 1.0
# Perform alignment via STAR two-pass methodology and analyze via DESeq2

#### WORKFLOW DEFINITION

workflow STAR2Pass {
  input {
    File batch_file
    File star_genome_tar 
    String ref_genome 
    File rnaseqc_genes_gtf
  }

  Array[Object] batch_info = read_objects(batch_file)

  scatter (job in batch_info){
    String sample_name = job.omics_sample_name
    String molecular_id = job.molecular_id
    String fastq_r1_locs = job.R1
    String fastq_r2_locs = job.R2
    String base_file_name = sample_name + "_" + molecular_id

    call FindFastqs {
      input:
        fastq_r1_str=fastq_r1_locs, 
        fastq_r2_str=fastq_r2_locs
    }

    call ConcatenateFastQs {
      input: 
        fastq_r1_array=FindFastqs.r1_locs, 
        fastq_r2_array=FindFastqs.r2_locs, 
        base_file_name=base_file_name
    }

    call STARalignTwoPass {
      input:
        base_file_name=base_file_name,
        star_genome_refs_zipped=star_genome_tar,
        r1fastq=ConcatenateFastQs.r1fastq,
        r2fastq=ConcatenateFastQs.r2fastq,
        ref_genome=ref_genome
    }

    call RNASeQC {
      input:
        base_file_name=base_file_name,
        bam_file=STARalignTwoPass.bam,
        bam_index=STARalignTwoPass.bai,
        ref_gtf=rnaseqc_genes_gtf
    }

  } # End scatter 

  # Outputs that will be retained when execution is complete
  output {
    Array[File] output_bam = STARalignTwoPass.bam
    Array[File] output_bai = STARalignTwoPass.bai
    Array[File] output_geneCounts = STARalignTwoPass.geneCounts
    Array[File] output_log_final = STARalignTwoPass.log_final
    Array[File] output_log_progress = STARalignTwoPass.log_progress
    Array[File] output_log = STARalignTwoPass.log
    Array[File] output_SJ = STARalignTwoPass.SJout
    Array[File] output_rnaseqc = RNASeQC.rnaseqc_metrics
  }

  parameter_meta {
    batch_file: "input tsv describing the samples to be analyzed and the locations of their fastq files"
    star_genome_tar: "reference genome files necessary for STAR analysis in tar file format"
    ref_genome: "name of the reference genome"
    rnaseqc_genes_gtf: "gene-based gtf annotation file providing the genomic location of each gene"

    output_bam: "array of aligned bam files for each sample"
    output_bai: "array of corresponding index files for each aligned bam file"
    output_gene_counts: "array of text files containing the number of reads in each gene for each sample"
    output_log_final: "array of text files containing an overarching summary of the analysis performed for each sample"
    output_log_progress: "array of text files containing a detailed progress report for each sample"
    output_log: "array of text files containing STAR's raw command line output for each sample"
    output_SJ: "array of text files containing splice junction details for each sample being analyzed"
    output_rnaseqc: "array of tar files containing RNA QC data for each sample being analyzed"
  }
} # End Workflow

#### TASK DEFINITIONS

# Locates the specified fastq files
task FindFastqs {
  input {
    String fastq_r1_str
    String fastq_r2_str
  }

  command <<<
    IFS="," read -ra ARR1 <<< "~{fastq_r1_str}"
    for item in "${ARR1[@]}"; do  
      echo "$item" >> R1out
    done
    IFS="," read -ra ARR2 <<< "~{fastq_r2_str}"
    for item in "${ARR2[@]}"; do  
      echo "$item" >> R2out
    done
  >>>

  output {
    Array[File] r1_locs = read_lines("R1out")
    Array[File] r2_locs = read_lines("R2out")
  }

  runtime {
    docker: "ubuntu:latest"
    memory: "2 GB"
    cpu: "2"
  }

  parameter_meta {
    fastq_r1_str: "comma-separated string of R1 fastq locations for the sample in question"
    fastq_r2_str: "comma-separated string of R2 fastq locations for the sample in question"

    r1_locs: "array of file objects corresponding to the provided R1 fastq's"
    r2_locs: "array of file objects corresponding to the provided R2 fastq's"
  }
}

# Concatenates the provided fastq files into a single fastq
task ConcatenateFastQs {
  input {
    String base_file_name
    Array[File] fastq_r1_array
    Array[File] fastq_r2_array
  }

  command <<<
    cat ~{sep=' ' fastq_r1_array} > "~{base_file_name}.R1.fastq.gz"
    cat ~{sep=' ' fastq_r2_array} > "~{base_file_name}.R2.fastq.gz"
  >>>

  output {
    File r1fastq = "~{base_file_name}.R1.fastq.gz"
    File r2fastq = "~{base_file_name}.R2.fastq.gz"
  }

  runtime {
    docker: "ubuntu:latest"
    memory: "2GB"
    cpu: "2"
  }

  parameter_meta {
    base_file_name: "base file name to use when saving results"
    fastq_r1_array: "array of R1 fastq files to concatenate"
    fastq_r2_array: "array of R2 fastq files to concatenate"

    r1fastq: "final concatenated R1 fastq file"
    r2fastq: "final concatenated R2 fastq file"
  }
}

# Aligns reads using STAR's two-pass methodology
task STARalignTwoPass {
  input {
    File star_genome_refs_zipped
    File r1fastq
    File r2fastq
    String base_file_name
    String ref_genome
  }

  String star_db_dir = basename(star_genome_refs_zipped, ".tar.gz")

  command <<<
    set -eo pipefail
    tar -xzf "~{star_genome_refs_zipped}"
    STAR \
      --genomeDir "~{star_db_dir}" \
      --readFilesIn "~{r1fastq}" "~{r2fastq}" \
      --runThreadN 6 \
      --readFilesCommand zcat \
      --sjdbOverhang 100 \
      --outSAMtype BAM SortedByCoordinate \
      --twopassMode Basic \
      --quantMode GeneCounts \
      --quantTranscriptomeBAMcompression 5 
    mv Aligned.sortedByCoord.out.bam "~{base_file_name}.~{ref_genome}.Aligned.sortedByCoord.out.bam"
    mv ReadsPerGene.out.tab "~{base_file_name}.~{ref_genome}.ReadsPerGene.out.tab"
    mv Log.final.out "~{base_file_name}.~{ref_genome}.Log.final.out"
    samtools index "~{base_file_name}.~{ref_genome}.Aligned.sortedByCoord.out.bam"
  >>>

  output {
    File bam = "~{base_file_name}.~{ref_genome}.Aligned.sortedByCoord.out.bam"
    File bai = "~{base_file_name}.~{ref_genome}.Aligned.sortedByCoord.out.bam.bai"
    File geneCounts = "~{base_file_name}.~{ref_genome}.ReadsPerGene.out.tab"
    File log_final = "~{base_file_name}.~{ref_genome}.Log.final.out"
    File log_progress = "Log.progress.out"
    File log = "Log.out"
    File SJout = "SJ.out.tab"
  }

  runtime {
    docker: "ghcr.io/getwilds/star:2.7.6a"
    memory: "62 GB"
    cpu: "8"
  }

  parameter_meta {
    star_genome_refs_zipped: "compressed zip file containing the reference files necessary for STAR"
    r1fastq: "R1 fastq containing raw reads to align via STAR"
    r2fastq: "R2 fastq containing raw reads to align via STAR"
    base_file_name: "base file name to use when saving results"
    ref_genome: "name of the reference genome being used"

    bam: "aligned bam files for the sample in question produced by STAR"
    bai: "corresponding index file for the bam file"
    geneCounts: "text file containing the number of reads in each gene"
    log_final: "text file containing an overarching summary of the analysis performed"
    log_progress: "text file containing a detailed progress report"
    log: "text file containing STAR's raw command line output"
    SJout: "text file containing splice junction details"
  }
}

# Calculates relevant RNA QC statistics using RNASeQC
# Make sure to collapse your exon-based gtf into a gene-based gtf via the script below:
# https://github.com/broadinstitute/gtex-pipeline/blob/master/gene_model/collapse_annotation.py
task RNASeQC {
  input {
    File bam_file
    File bam_index
    File ref_gtf
    String base_file_name
  }

  command <<<
    rnaseqc "~{ref_gtf}" "~{bam_file}" OUTPUT \
      --sample="~{base_file_name}" --coverage 
    tar -cvzf "~{base_file_name}.QC.tar.gz" OUTPUT/*
  >>>

  output {
    File rnaseqc_metrics = "~{base_file_name}.QC.tar.gz"
  }

  runtime {
    docker: "ghcr.io/getwilds/rnaseqc:2.4.2"
    memory: "4 GB"
    cpu: "2"
  }

  parameter_meta {
    bam_file: "aligned bam file to be QC-analyzed"
    bam_index: "corresponding index file for the bam file"
    ref_gtf: "gene-based gtf annotation file providing the genomic location of each gene"
    base_file_name: "base file name to use when saving results"

    rnaseqc_metrics: "tar file containing RNA QC data for the sample in question"
  }
}

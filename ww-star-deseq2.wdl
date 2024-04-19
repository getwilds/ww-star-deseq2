version 1.0
# Perform alignment via STAR two-pass methodology and analyze via DESeq2

#### WORKFLOW DEFINITION

workflow STAR2Pass {
  input {
    File batchFile
    File STARgenomeTAR 
    String referenceGenome 
    File RNASeQC_genesGtf
  }

  Array[Object] batchInfo = read_objects(batchFile)

  scatter (job in batchInfo){
    String sampleName = job.omics_sample_name
    String molecular_id = job.molecular_id
    String fastqR1Locations = job.R1
    String fastqR2Locations = job.R2
    String base_file_name = sampleName + "_" + molecular_id

    call FindFastqs {
      input:
        fastqR1String=fastqR1Locations, 
        fastqR2String=fastqR2Locations
    }

    call ConcatenateFastQs {
      input: 
        fastqR1Array=FindFastqs.R1Locations, 
        fastqR2Array=FindFastqs.R2Locations, 
        base_file_name=base_file_name
    }

    call STARalignTwoPass {
      input:
        base_file_name=base_file_name,
        star_genome_refs_zipped=STARgenomeTAR,
        r1fastq=ConcatenateFastQs.R1fastq,
        r2fastq=ConcatenateFastQs.R2fastq,
        referenceGenome=referenceGenome
    }

    call RNASeQC {
      input:
        base_file_name=base_file_name,
        bam_file=STARalignTwoPass.bam,
        bam_index=STARalignTwoPass.bai,
        refGtf=RNASeQC_genesGtf
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
    batchFile: "input tsv describing the samples to be analyzed and the locations of their fastq files"
    STARgenomeTAR: "reference genome files necessary for STAR analysis in tar file format"
    referenceGenome: "name of the reference genome"
    RNASeQC_genesGtf: "gene-based gtf annotation file providing the genomic location of each gene"

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
    String fastqR1String
    String fastqR2String
  }

  command <<<
    IFS="," read -ra ARR1 <<< "~{fastqR1String}"
    for item in "${ARR1[@]}"; do  
      echo "$item" >> R1out
    done
    IFS="," read -ra ARR2 <<< "~{fastqR2String}"
    for item in "${ARR2[@]}"; do  
      echo "$item" >> R2out
    done
  >>>

  output {
    Array[File] R1Locations = read_lines("R1out")
    Array[File] R2Locations = read_lines("R2out")
  }

  runtime {
    docker: "ubuntu:latest"
    memory: "2 GB"
    cpu: "2"
  }

  parameter_meta {
    fastqR1String: "comma-separated string of R1 fastq locations for the sample in question"
    fastqR2String: "comma-separated string of R2 fastq locations for the sample in question"

    R1Locations: "array of file objects corresponding to the provided R1 fastq's"
    R2Locations: "array of file objects corresponding to the provided R2 fastq's"
  }
}

# Concatenates the provided fastq files into a single fastq
task ConcatenateFastQs {
  input {
    String base_file_name
    Array[File] fastqR1Array
    Array[File] fastqR2Array
  }

  command <<<
    cat ~{sep=' ' fastqR1Array} > "~{base_file_name}.R1.fastq.gz"
    cat ~{sep=' ' fastqR2Array} > "~{base_file_name}.R2.fastq.gz"
  >>>

  output {
    File R1fastq = "~{base_file_name}.R1.fastq.gz"
    File R2fastq = "~{base_file_name}.R2.fastq.gz"
  }

  runtime {
    docker: "ubuntu:latest"
    memory: "2GB"
    cpu: "2"
  }

  parameter_meta {
    base_file_name: "base file name to use when saving results"
    fastqR1Array: "array of R1 fastq files to concatenate"
    fastqR2Array: "array of R2 fastq files to concatenate"

    R1fastq: "final concatenated R1 fastq file"
    R2fastq: "final concatenated R2 fastq file"
  }
}

# Aligns reads using STAR's two-pass methodology
task STARalignTwoPass {
  input {
    File star_genome_refs_zipped
    File r1fastq
    File r2fastq
    String base_file_name
    String referenceGenome
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
    mv Aligned.sortedByCoord.out.bam "~{base_file_name}.~{referenceGenome}.Aligned.sortedByCoord.out.bam"
    mv ReadsPerGene.out.tab "~{base_file_name}.~{referenceGenome}.ReadsPerGene.out.tab"
    mv Log.final.out "~{base_file_name}.~{referenceGenome}.Log.final.out"
    samtools index "~{base_file_name}.~{referenceGenome}.Aligned.sortedByCoord.out.bam"
  >>>

  output {
    File bam = "~{base_file_name}.~{referenceGenome}.Aligned.sortedByCoord.out.bam"
    File bai = "~{base_file_name}.~{referenceGenome}.Aligned.sortedByCoord.out.bam.bai"
    File geneCounts = "~{base_file_name}.~{referenceGenome}.ReadsPerGene.out.tab"
    File log_final = "~{base_file_name}.~{referenceGenome}.Log.final.out"
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
    referenceGenome: "name of the reference genome being used"

    bam: ""
    bai: ""
    geneCounts: ""
    log_final: ""
    log_progress: ""
    log: ""
    SJout: ""
  }
}

# Calculates relevant RNA QC statistics using RNASeQC
# Make sure to collapse your exon-based gtf into a gene-based gtf via the script below:
# https://github.com/broadinstitute/gtex-pipeline/blob/master/gene_model/collapse_annotation.py
task RNASeQC {
  input {
    File bam_file
    File bam_index
    File refGtf
    String base_file_name
  }

  command <<<
    rnaseqc "~{refGtf}" "~{bam_file}" OUTPUT \
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
    refGtf: "gene-based gtf annotation file providing the genomic location of each gene"
    base_file_name: "base file name to use when saving results"

    rnaseqc_metrics: "tar file containing RNA QC data for the sample in question"
  }
}

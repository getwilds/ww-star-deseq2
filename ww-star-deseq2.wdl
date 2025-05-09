version 1.0

struct SampleInfo {
    String name
    File r1
    File r2
    String condition
}

struct RefGenome {
    String name
    File fasta
    File gtf
}

workflow star_deseq2 {
  meta {
    author: "Taylor Firman"
    email: "tfirman@fredhutch.org"
    description: "WDL workflow for RNA-seq alignment via STAR and DESeq2 differential expression analysis"
    url: "https://github.com/getwilds/ww-star-deseq2"
    outputs: {
        star_bam: "",
        star_bai: "",
        star_gene_counts: "",
        star_log_final: "",
        star_log_progress: "",
        star_log: "",
        star_sj: "",
        rnaseqc_metrics: "",
        combined_counts_matrix: "",
        sample_metadata: "",
        deseq2_all_results: "",
        deseq2_significant_results: "",
        deseq2_normalized_counts: "",
        deseq2_pca_plot: "",
        deseq2_volcano_plot: "",
        deseq2_heatmap: ""
    }
  }

  parameter_meta {
    samples: "list of samples to be analyzed by STAR and DESeq2"
    reference_genome: "reference genome to use during STAR alignment"
    reference_level: ""
    contrast: ""
  }

  input {
    Array[SampleInfo] samples
    RefGenome? reference_genome
    String reference_level = ""
    String contrast = ""
  }

  if (!defined(reference_genome)) {
    call download_reference {}
  }

  RefGenome ref_genome_final = select_first([reference_genome, download_reference.genome])

  call collapse_gtf { input:
      reference_gtf = ref_genome_final.gtf
  }

  call build_star_index { input:
      reference_fasta = ref_genome_final.fasta,
      reference_gtf = ref_genome_final.gtf
  }

  scatter (sample in samples) {
    call star_align_two_pass { input:
        sample_data = sample,
        star_genome_tar = build_star_index.star_index_tar,
        ref_genome_name = ref_genome_final.name
    }

    call rnaseqc_cov { input:
        base_file_name = sample.name,
        bam_file = star_align_two_pass.bam,
        bam_index = star_align_two_pass.bai,
        ref_gtf = collapse_gtf.collapsed_gtf
    }
  }

  call combine_count_matrices { input:
      gene_count_files = star_align_two_pass.gene_counts,
      sample_names = star_align_two_pass.name,
      sample_conditions = star_align_two_pass.condition
  }

  call run_deseq2 { input:
      counts_matrix = combine_count_matrices.counts_matrix,
      sample_metadata = combine_count_matrices.sample_metadata,
      reference_level = reference_level,
      contrast = contrast
  }

  output {
    Array[File] star_bam = star_align_two_pass.bam
    Array[File] star_bai = star_align_two_pass.bai
    Array[File] star_gene_counts = star_align_two_pass.gene_counts
    Array[File] star_log_final = star_align_two_pass.log_final
    Array[File] star_log_progress = star_align_two_pass.log_progress
    Array[File] star_log = star_align_two_pass.log
    Array[File] star_sj = star_align_two_pass.sj_out
    Array[File] rnaseqc_metrics = rnaseqc_cov.rnaseqc_metrics
    File combined_counts_matrix = combine_count_matrices.counts_matrix
    File sample_metadata = combine_count_matrices.sample_metadata
    File deseq2_all_results = run_deseq2.deseq2_results
    File deseq2_significant_results = run_deseq2.deseq2_significant
    File deseq2_normalized_counts = run_deseq2.deseq2_normalized_counts
    File deseq2_pca_plot = run_deseq2.deseq2_pca_plot
    File deseq2_volcano_plot = run_deseq2.deseq2_volcano_plot
    File deseq2_heatmap = run_deseq2.deseq2_heatmap
  }
}

task download_reference {
  meta {
    description: "Task for pulling down the GRCh38 reference fasta and gtf."
    outputs: {
        genome: "reference genome details to be used throughout the workflow"
    }
  }

  parameter_meta {}

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
    docker: "getwilds/gtf-smash:v8"
    memory: "4 GB"
    cpu: 1
  }
}

task build_star_index {
  meta {
    description: "Task for building the STAR index files from fasta/gtf."
    outputs: {
        star_index_tar: ""
    }
  }

  parameter_meta {
    reference_fasta: ""
    reference_gtf: ""
    sjdb_overhang: ""
    genome_sa_index_nbases: ""
    memory_gb: ""
    cpu_cores: ""
  }

  input {
    File reference_fasta
    File reference_gtf
    Int sjdb_overhang = 100
    Int genome_sa_index_nbases = 14
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
      --genomeFastaFiles "~{reference_fasta}" \
      --sjdbGTFfile "~{reference_gtf}" \
      --sjdbOverhang ~{sjdb_overhang} \
      --genomeSAindexNbases ~{genome_sa_index_nbases}

    tar -czf star_index.tar.gz star_index/
  >>>

  output {
    File star_index_tar = "star_index.tar.gz"
  }

  runtime {
    docker: "getwilds/star:2.7.6a"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task collapse_gtf {
  meta {
    description: "Task for collapsing the provided gtf into one transcript per gene."
    outputs: {
        collapsed_gtf: ""
    }
  }

  parameter_meta {
    reference_gtf: ""
  }

  input {
    File reference_gtf
  }

  command <<<
    set -eo pipefail
    
    echo "Processing GTF file..."
    collapse_annotation.py \
      "~{reference_gtf}" \
      collapsed.gtf
  >>>

  output {
    File collapsed_gtf = "collapsed.gtf"
  }

  runtime {
    docker: "getwilds/gtf-smash:v8"
    memory: "4 GB"
    cpu: 1
  }
}

task star_align_two_pass {
  meta {
    description: "Task for aligning RNA-seq reads using STAR's two-pass technique."
    outputs: {
        name: "",
        condition: "",
        bam: "",
        bai: "",
        gene_counts: "",
        log_final: "",
        log_progress: "",
        log: "",
        sj_out: ""
    }
  }

  parameter_meta {
    sample_data: ""
    star_genome_tar: ""
    ref_genome_name: ""
    sjdb_overhang: ""
    memory_gb: ""
    cpu_cores: ""
    star_threads: ""
  }

  input {
    SampleInfo sample_data
    File star_genome_tar
    String ref_genome_name
    Int sjdb_overhang = 100
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
      --readFilesIn "~{sample_data.r1}" "~{sample_data.r2}" \
      --runThreadN ~{star_threads} \
      --readFilesCommand zcat \
      --sjdbOverhang ~{sjdb_overhang} \
      --outSAMtype BAM SortedByCoordinate \
      --twopassMode Basic \
      --outTmpDir _STARtmp \
      --outFileNamePrefix "./" \
      --quantMode GeneCounts \
      --quantTranscriptomeBAMcompression 5 

    rm -r star_index _STARtmp

    mv Aligned.sortedByCoord.out.bam "~{sample_data.name}.~{ref_genome_name}.Aligned.sortedByCoord.out.bam"
    mv ReadsPerGene.out.tab "~{sample_data.name}.~{ref_genome_name}.ReadsPerGene.out.tab"
    mv Log.final.out "~{sample_data.name}.~{ref_genome_name}.Log.final.out"
    mv Log.progress.out "~{sample_data.name}.~{ref_genome_name}.Log.progress.out"
    mv Log.out "~{sample_data.name}.~{ref_genome_name}.Log.out"
    mv SJ.out.tab "~{sample_data.name}.~{ref_genome_name}.SJ.out.tab"

    samtools index "~{sample_data.name}.~{ref_genome_name}.Aligned.sortedByCoord.out.bam"
  >>>

  output {
    String name = sample_data.name
    String condition = sample_data.condition
    File bam = "~{sample_data.name}.~{ref_genome_name}.Aligned.sortedByCoord.out.bam"
    File bai = "~{sample_data.name}.~{ref_genome_name}.Aligned.sortedByCoord.out.bam.bai"
    File gene_counts = "~{sample_data.name}.~{ref_genome_name}.ReadsPerGene.out.tab"
    File log_final = "~{sample_data.name}.~{ref_genome_name}.Log.final.out"
    File log_progress = "~{sample_data.name}.~{ref_genome_name}.Log.progress.out"
    File log = "~{sample_data.name}.~{ref_genome_name}.Log.out"
    File sj_out = "~{sample_data.name}.~{ref_genome_name}.SJ.out.tab"
  }

  runtime {
    docker: "getwilds/star:2.7.6a"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task rnaseqc_cov {
  meta {
    description: "Task for aligning RNA-seq reads using STAR's two-pass technique."
    outputs: {
        rnaseqc_metrics: ""
    }
  }

  parameter_meta {
    bam_file: ""
    bam_index: ""
    ref_gtf: ""
    base_file_name: ""
    memory_gb: ""
    cpu_cores: ""
  }

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
    cpu: cpu_cores
  }
}

task combine_count_matrices {
  input {
    Array[File] gene_count_files
    Array[String] sample_names
    Array[String] sample_conditions
    Int memory_gb = 4
    Int cpu_cores = 1
    Int count_column = 2
        # Column to extract from ReadsPerGene.out.tab files:
        # 2 = unstranded counts
        # 3 = stranded counts, first read forward
        # 4 = stranded counts, first read reverse
  }

  command <<<
    set -eo pipefail

    combine_star_counts.py \
      --input "~{sep=' ' gene_count_files}" \
      --names "~{sep=' ' sample_names}" \
      --conditions "~{sep=' ' sample_conditions}" \
      --output combined_counts_matrix.txt \
      --count_column ~{count_column}
  >>>

  output {
    File counts_matrix = "combined_counts_matrix.txt"
    File sample_metadata = "sample_metadata.txt"
  }

  runtime {
    docker: "getwilds/combine-counts:0.1.0"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

task run_deseq2 {
  meta {
    description: "Task for analyzing differential expression via DESeq2."
    outputs: {
        deseq2_results: "",
        deseq2_significant: "",
        deseq2_normalized_counts: "",
        deseq2_pca_plot: "",
        deseq2_volcano_plot: "",
        deseq2_heatmap: ""
    }
  }

  parameter_meta {
    counts_matrix: ""
    sample_metadata: ""
    condition_column: ""
    reference_level: ""
    contrast: ""
    memory_gb: ""
    cpu_cores: ""
  }

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
    File deseq2_normalized_counts = "deseq2_results_normalized_counts.csv"
    File deseq2_pca_plot = "deseq2_results_pca.pdf"
    File deseq2_volcano_plot = "deseq2_results_volcano.pdf"
    File deseq2_heatmap = "deseq2_results_heatmap.pdf"
  }

  runtime {
    docker: "getwilds/deseq2:1.40.2"
    memory: "~{memory_gb} GB"
    cpu: cpu_cores
  }
}

workflow STAR2Pass {
  # Note read length for STAR reference was 100bp so currently this assumes your sequencing is that long too. 
  File batchFile
  Array[Object] batchInfo = read_objects(batchFile)
  File STARgenomeTAR 
  String referenceGenome 
  File RNASeQC_genesGtf

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

    # call RNASeQC {
    #   input:
    #     base_file_name=base_file_name,
    #     bam_file=STARalignTwoPass.bam,
    #     bam_index=STARalignTwoPass.bai,
    #     refGtf=RNASeQC_genesGtf
    # }

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
    ##Array[File] output_rnaseqc = RNASeQC.rnaseqc_metrics
  }
} # End Workflow

## TASK DEFINITIONS
task FindFastqs {
  String fastqR1String
  String fastqR2String
  String tempR1array = "{ARR1[@]}"
  String tempR2array = "{ARR2[@]}"
  command {
    IFS="," read -ra ARR1 <<< "${fastqR1String}"
    for item in "$${tempR1array}"; do  
      echo "$item" >> R1out
    done

    IFS="," read -ra ARR2 <<< "${fastqR2String}"
    for item in "$${tempR2array}"; do  
      echo "$item" >> R2out
    done
  }
  runtime {
    docker: "ubuntu:latest"
    memory: "2 GB"
    cpu: "2"
  }
  output {
    Array[File] R1Locations = read_lines("R1out")
    Array[File] R2Locations = read_lines("R2out")
  }
}

task ConcatenateFastQs {
  String base_file_name
  Array[File] fastqR1Array
  Array[File] fastqR2Array
  command {
    cat ${sep=' ' fastqR1Array} > ${base_file_name}.R1.fastq.gz
    cat ${sep=' ' fastqR2Array} > ${base_file_name}.R2.fastq.gz
   }
  runtime {
    docker: "ubuntu:latest"
    memory: "2GB"
    cpu: "2"
  }
  output {
    File R1fastq = "${base_file_name}.R1.fastq.gz"
    File R2fastq = "${base_file_name}.R2.fastq.gz"
  }
}

task STARalignTwoPass {
	File star_genome_refs_zipped
	File r1fastq
	File r2fastq
	String base_file_name
  String referenceGenome
	command {
	  set -e

    mkdir genomeDir
	  tar -xzf ${star_genome_refs_zipped} -C genomeDir

    STAR \
      --genomeDir genomeDir \
      --readFilesIn ${r1fastq} ${r2fastq} \
      --runThreadN 6 \
      --readFilesCommand zcat \
      --sjdbOverhang 100 \
      --outSAMtype BAM SortedByCoordinate \
      --twopassMode Basic \
      --quantMode GeneCounts \
      --quantTranscriptomeBAMcompression 5 

    mv Aligned.sortedByCoord.out.bam ${base_file_name}.${referenceGenome}.Aligned.sortedByCoord.out.bam
    mv ReadsPerGene.out.tab ${base_file_name}.${referenceGenome}.ReadsPerGene.out.tab
    mv Log.final.out ${base_file_name}.${referenceGenome}.Log.final.out
    samtools index ${base_file_name}.${referenceGenome}.Aligned.sortedByCoord.out.bam
	}
	output {
		File bam = "${base_file_name}.${referenceGenome}.Aligned.sortedByCoord.out.bam"
    File bai = "${base_file_name}.${referenceGenome}.Aligned.sortedByCoord.out.bam.bai"
    File geneCounts = "${base_file_name}.${referenceGenome}.ReadsPerGene.out.tab"
		File log_final = "${base_file_name}.${referenceGenome}.Log.final.out"
    File log_progress = "Log.progress.out"
		File log = "Log.out"
    File SJout = "SJ.out.tab"
	}
	runtime {
		docker: "fredhutch/star:2.7.1a"
		memory: "62 GB"
		cpu: "8"
	}
}
### r5.2xlarge = 8 cpu, 64GiB memory

task RNASeQC {
  File bam_file
  File bam_index
  File refGtf
  String base_file_name
  command {
    java -Xmx4g -jar /path/to/RNA-SeQC.jar -n 1000 \
      ${refGtf} ${bam_file} OUTPUT \
      --sample=${base_file_name} \
      --coverage 

      tar -cvzf ${base_file_name}.QC.tar.gz OUTPUT/*
  }
  output {
    File rnaseqc_metrics = "${base_file_name}.QC.tar.gz"
  }
  runtime {
    docker: "broadinstitute/gtex_rnaseq"
    memory: "4 GB"
    cpu: "2"
  }
}

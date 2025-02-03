/*
   * Save the seedfile into the straindb database
   cp $seed ${seed.baseName}.copy.tsv
   */
  process save_seedfile {

      publishDir "s3://genomics-workflow-core/aws-miti-straindb-us-west-2/aws_glue/assembly_seedfiles/"
      container params.container

      input:
      path seed
      output:
      path "${seed}"

     script:
     """
      head $seed
     """
  }

  /*
   * Parse software version numbers
   */
  process get_software_versions {

      errorStrategy 'ignore'
      container params.container
      cpus 2
      memory 8.GB

      publishDir "${params.output_path}/${sample}/${prefix}/software_info", mode: 'copy', pattern: "*.{txt}"

      input:
      tuple  val(sample), val(prefix), val(reads1), val(reads2), val(long_reads)
      output:
      file "*.txt"

      script:
      // TODO nf-core: Get all tools to print their version number here
      // echo "Git info: $workflow.repository - $workflow.revision [$workflow.commitId]" > v_repo.txt
      """
      echo $workflow.manifest.version > v_pipeline.txt
      echo $workflow.nextflow.version > v_nextflow.txt
      """
  }


  /*
   * Run Hybridassembly Pipeline -long reads only
   */
  process hybridassembly_long {
      tag "$sample"
      //container "xianmeng/nf-hybridassembly:latest"
      container params.container
      cpus 32
      memory 64.GB

      //publishDir "${output_path}/${sample}/${prefix}/UNICYCLER_LONG", mode:'copy', pattern: "*.{log,gfa,fasta}"

      input:
      tuple  val(sample), val(prefix), val(reads1), val(reads2), val(long_reads)


      output:
      //path "tmp_*/Sync/UNICYCLER/*"
      //file "tmp_*/Sync/UNICYCLER/assembly.fasta"
      tuple val(sample), val(prefix), val(reads1), val(reads2), val(long_reads), path("tmp_*/Sync/UNICYCLER/assembly.fasta")

      script:
      """
      export sampleRate="${params.sampleRate}"
      export longreads="${long_reads}"
      export S3OUTPUTPATH="${params.output_path}/${sample}/${prefix}/UNICYCLER_LONG"
      run_unicycler_long_only2.sh
      """
  }



  //seedfile_ch.view()
  /*
   * Run Hybridassembly Pipeline
   */
  process hybridassembly {

      tag "$sample"
      //container "xianmeng/nf-hybridassembly:latest"
      container params.container
      cpus 32
      memory 128.GB

      input:
    	//tuple val(sample), val(reads1), val(reads2), val(long_reads) from seedfile_ch2
      tuple val(sample), val(prefix), val(reads1), val(reads2), val(long_reads), path(exist_long)


      output:
      tuple val(sample), val(prefix), path("tmp_*/raw_fastq/read1.fastq.gz"), path("tmp_*/raw_fastq/read2.fastq.gz"), path("tmp_*/Sync/UNICYCLER/assembly.fasta")

      script:
      """
      export SAMPLE_NAME="${sample}"
      export coverage="${params.coverage}"
      export fastq1="${reads1}"
      export fastq2="${reads2}"
      export longreads="${long_reads}"
      export assembly="${exist_long}"
      export S3OUTPUTPATH="${params.output_path}/${sample}/${prefix}"
      run_unicycler_existing_long2.sh

      """
  }

  process hybridassembly_bwa {
      tag "$sample"
      //container "xianmeng/nf-hybridassembly:latest"
      container params.container_bwa
      cpus 16
      memory 32.GB

      //publishDir "${output_path}/${sample}/${prefix}/UNICYCLER_POLISHED", mode:'copy', pattern: "*.{fasta}"

      input:
      tuple val(sample), val(prefix), path(reads1), path(reads2), path(draft_assembly)

      output:
      tuple val(sample), val(prefix), path("alignments_1.sam"),  path("alignments_2.sam"), path(draft_assembly)

      script:
      """
      bwa index ${draft_assembly}
      bwa mem -t 16 -a ${draft_assembly} ${reads1} > alignments_1.sam
      bwa mem -t 16 -a ${draft_assembly} ${reads2} > alignments_2.sam
      rm *.amb *.ann *.bwt *.pac *.sa
      """
  }


  process hybridassembly_polish {
      tag "$sample"
      //container "xianmeng/nf-hybridassembly:latest"
      container params.container_polish
      cpus 16
      memory 32.GB

      publishDir "${params.output_path}/${sample}/${prefix}/UNICYCLER_POLISHED", mode:'copy', pattern: "*.{fasta,log}"

      input:
      tuple val(sample), val(prefix),path(alignments_1), path(alignments_2), path(draft_assembly)

      output:
      file "${sample}.fasta"
      file "Polypolish.log"

      script:
      """
      polypolish filter --in1 ${alignments_1} --in2 ${alignments_2} --out1 filtered_1.sam --out2 filtered_2.sam
      polypolish polish ${draft_assembly} filtered_1.sam filtered_2.sam > ${sample}.fasta 2>Polypolish.log
      """
  }

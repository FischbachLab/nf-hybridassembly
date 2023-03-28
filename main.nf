#!/usr/bin/env nextflow
nextflow.enable.dsl=1
// If the user uses the --help flag, print the help text below
params.help = false

// Function which prints help message text
def helpMessage() {
    log.info"""
    Run the Hybrid assembly pipeline for a given short and long read dataset

    Required Arguments:
      --seedfile      file        a file contains reads1, reads2 and long reads for assembly
      --output_path   path        Output s3 path

    Options:
      --coverage      num         Reads coverage depth (e.g. 200)
      -profile        docker      run locally


    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

if (params.output_path == "null") {
	exit 1, "Missing the output path"
}

Channel
  .fromPath(params.seedfile)
  .ifEmpty { exit 1, "Cannot find the input seedfile" }

/*
 * Defines the pipeline inputs parameters (giving a default value for each for them)
 * Each of the following parameters can be specified as command line options
 */

def output_path = "${params.output_path}"
//def output_path=s3://genomics-workflow-core/Pipeline_Results/${params.output_prefix}"

//println output_path
//	.map{ row -> tuple(row.sample, row.reads1, row.reads2, row.long_reads, row.exist_long)}
Channel
	.fromPath(params.seedfile)
	.ifEmpty { exit 1, "Cannot find any seed file matching: ${params.seedfile}." }
  .splitCsv(header: ['sample', 'prefix', 'reads1', 'reads2', 'long_reads'], sep: '\t')
	.map{ row -> tuple(row.sample, row.prefix, row.reads1, row.reads2, row.long_reads)}
	.set { seedfile_ch }

  seedfile_ch.into {seedfile_ch1; seedfile_ch2}

  /*
   * Save the seedfile into the straindb database
   cp $seed ${seed.baseName}.copy.tsv
   */
  process save_seedfile {

      publishDir "s3://genomics-workflow-core/aws-miti-straindb-us-west-2/aws_glue/assembly_seedfiles/"
      container params.container
      
      input:
      file seed from  Channel.fromPath(params.seedfile)
      output:
      file "${seed}"

     script:
     """
      ls $seed
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

      publishDir "${output_path}/${sample}/${prefix}/software_info", mode: 'copy'

      input:
      tuple  val(sample), val(prefix), val(reads1), val(reads2), val(long_reads) from seedfile_ch1
      output:
      file "*"

      script:
      // TODO nf-core: Get all tools to print their version number here
      """
      echo $workflow.manifest.version > v_pipeline.txt
      echo $workflow.nextflow.version > v_nextflow.txt
      echo "Git info: $workflow.repository - $workflow.revision [$workflow.commitId]" > v_repo.txt
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
      memory 128.GB

      publishDir "${output_path}/${sample}/${prefix}/UNICYCLER_LONG", mode:'copy', pattern: "*.{log,gfa,fasta}"

      input:
    	tuple  val(sample), val(prefix), val(reads1), val(reads2), val(long_reads) from seedfile_ch2


      output:
      path "tmp_*/Sync/UNICYCLER/*"
      //file "tmp_*/Sync/UNICYCLER/assembly.fasta"
      tuple val(sample), val(prefix), val(reads1), val(reads2), val(long_reads), path("tmp_*/Sync/UNICYCLER/assembly.fasta") into long_ch

      script:
      """
      export sampleRate="${params.sampleRate}"
      export longreads="${long_reads}"
      export S3OUTPUTPATH="${output_path}/${sample}/${prefix}"
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
      tuple val(sample), val(prefix), val(reads1), val(reads2), val(long_reads), path(exist_long) from long_ch


      output:
      //path "*"

      script:
      """
      export SAMPLE_NAME="${sample}"
      export coverage="${params.coverage}"
      export fastq1="${reads1}"
      export fastq2="${reads2}"
      export longreads="${long_reads}"
      export assembly="${exist_long}"
      export S3OUTPUTPATH="${output_path}/${sample}/${prefix}"
      run_unicycler_existing_long2.sh

      """
  }

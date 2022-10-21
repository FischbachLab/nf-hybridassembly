#!/usr/bin/env nextflow

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

Channel
	.fromPath(params.seedfile)
	.ifEmpty { exit 1, "Cannot find any seed file matching: ${params.seedfile}." }
  .splitCsv(header: ['sample', 'long_reads'], sep: '\t')
	.map{ row -> tuple(row.sample, row.long_reads)}
	.set { seedfile_ch }

  //seedfile_ch.view()
  /*
   * Run Hybridassembly Pipeline
   */
  process hybridassembly {

      //container "xianmeng/nf-hybridassembly:latest"
      container params.container
      cpus 32
      memory 128.GB

      publishDir "${output_path}", mode:'copy'

      input:
    	tuple val(sample), val(long_reads) from seedfile_ch


      output:
      //path "*"

      script:
      """
      export sampleRate="${params.sampleRate}"
      export longreads="${long_reads}"
      export S3OUTPUTPATH="${output_path}/${sample}"
      run_unicycler_long_only.sh
      """
  }

/*
 * Run Hybridassembly Pipeline

process hybridassembly_2 {

    //container "xianmeng/nf-hybridassembly:latest"
    container "fischbachlab/nf-hybridassembly:latest"
    cpus 32
    memory 128.GB

    publishDir "${output_path}", mode:'copy'

    input:
    //file read1 from read1_ch
    //file read2 from read2_ch

    output:
    //path "*"

    script:
    """
    export coverage="${params.coverage}"
    export fastq1="${params.reads1}"
    export fastq2="${params.reads2}"
    export longreads="${params.long_reads}"
    export S3OUTPUTPATH="${output_path}"
    run_unicycler_hybrid_nanopore.sh
    """
}

 */

#!/usr/bin/env nextflow

// If the user uses the --help flag, print the help text below
params.help = false

// Function which prints help message text
def helpMessage() {
    log.info"""
    Run the Hybrid assembly pipeline for a given short and long read dataset

    Required Arguments:
      --reads1        R1          Forward reads file path ( paired-end library )
      --reads2        R2          Reverse reads file path ( paired-end library )
      --long_reads    long reads  Long reads file path ( Nanopore or PacBio )
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
  .fromPath(params.reads1)
  .ifEmpty { exit 1, "Cannot find fastq R1 file" }

Channel
  .fromPath(params.reads2)
  .ifEmpty { exit 1, "Cannot find fastq R2 file" }

Channel
    .fromPath(params.long_reads)
    .ifEmpty { exit 1, "Cannot find fastq long reads file" }


/*
 * Defines the pipeline inputs parameters (giving a default value for each for them)
 * Each of the following parameters can be specified as command line options
 */

def output_path = "${params.output_path}"
//def output_path=s3://genomics-workflow-core/Pipeline_Results//${params.output_prefix}"

//println output_path
/*
Channel
    .fromPath(params.reads1)
    .set { read1_ch }

Channel
    .fromPath(params.reads2)
    .set { read2_ch }
*/


/*
 * Run Hybridassembly Pipeline
 */
process hybridassembly {

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
    export coverge="${params.coverage}"
    export fastq1="${params.reads1}"
    export fastq2="${params.reads2}"
    export longreads="${params.long_reads}"
    export S3OUTPUTPATH="${output_path}"
    run_unicycler_hybrid_nanopore.sh
    """
}

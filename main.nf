#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include {save_seedfile; get_software_versions; hybridassembly_long; hybridassembly; hybridassembly_bwa; hybridassembly_polish} from './modules/hybrid_assembly'


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
  

  Channel
    .fromPath(params.seedfile)
    .ifEmpty { exit 1, "Cannot find any seed file matching: ${params.seedfile}." }
    .splitCsv(header: ['sample', 'prefix', 'reads1', 'reads2', 'long_reads'], sep: '\t')
    .map{ row -> tuple(row.sample, row.prefix, row.reads1, row.reads2, row.long_reads)}
    .set { seedfile_ch }

  workflow {

    seedfile_ch | get_software_versions
    Channel.fromPath(params.seedfile) | save_seedfile
    seedfile_ch | hybridassembly_long

    hybridassembly_long.out | hybridassembly
    hybridassembly.out | hybridassembly_bwa
    hybridassembly_bwa.out | hybridassembly_polish

  }

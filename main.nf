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

Channel
	.fromPath(params.seedfile)
	.ifEmpty { exit 1, "Cannot find any seed file matching: ${params.seedfile}." }
  .splitCsv(header: ['sample', 'reads1', 'reads2', 'long_reads'], sep: '\t')
	.map{ row -> tuple(row.sample, row.reads1, row.reads2, row.long_reads)}
	.set { seedfile_ch }

  //seedfile_ch.view()
  /*
   * Run Hybridassembly Pipeline
   */
  process hybridassembly_1 {

      tag "$sample"
      //container "xianmeng/nf-hybridassembly:latest"
      container params.container
      cpus   32
      memory 128.GB

      //publishDir "${output_path}", mode:'copy'

      input:
    	tuple val(sample), val(reads1), val(reads2), val(long_reads) from seedfile_ch

      output:
      tuple val(sample),  path("tmp_*/trimmed_fastq/read1_sampled.fastq.gz"), path("tmp_*/trimmed_fastq/read2_sampled.fastq.gz"), path("tmp_*/trimmed_fastq/long_trimmed.fastq.gz") into reads_ch

      script:
      """
      export coverage="${params.coverage}"
      export fastq1="${reads1}"
      export fastq2="${reads2}"
      export longreads="${long_reads}"
      export S3OUTPUTPATH="${output_path}/${sample}"
      run_unicycler_hybrid_nanopore_1.sh
      """
  }

/*
 * Run Hybridassembly Pipeline -step 2
 * Assembly
 file tmp_/trimmed_fastq/read1_sampled.fastq.gz into read1_ch2
 *file tmp_/trimmed_fastq/read1_sampled.fastq.gz into read2_ch2
 file tmp_/trimmed_fastq/long_trimmed.fastq.gz  into longread_ch2
 file tmp_/Sync/UNICYCLER/assembly.fasta into assembly_ch2
*/

process hybridassembly_2 {

    tag "${sample}"
    container "quay.io/biocontainers/unicycler:0.5.0--py39h2add14b_2"
    cpus 16
    memory 64.GB

    publishDir "${output_path}/${sample}/UNICYCLER", mode:'copy', pattern: "*.{log,gfa,fasta}"
    publishDir "${output_path}/${sample}/UNICYCLER", mode:'copy', pattern: "assembly.*"
    //publishDir "${output_path}/${sample}/UNICYCLER", mode:'copy', pattern: "*.{log,gfa}"
    //publishDir "${output_path}/${sample}/Logs", mode:'copy', pattern: '*.txt'

    input:
  	tuple val(sample), path(read1), path(read2), path(long_read) from reads_ch

    output:
    file "tmp_*/Sync/UNICYCLER/*"
    file "tmp_*/Sync/UNICYCLER/assembly.gfa"
    //file "tmp_*/Sync/Logs/unicycler_assembly.log.txt"
    tuple val(sample),  path("tmp_*/trimmed_fastq/read1_sampled.fastq.gz"), path("tmp_*/trimmed_fastq/read2_sampled.fastq.gz"), path("tmp_*/trimmed_fastq/long_trimmed.fastq.gz"), path("tmp_*/Sync/UNICYCLER/assembly.fasta"), path("tmp_*/Sync/UNICYCLER/unicycler.log") into reads_ch2

    script:
    """
    export fastq1="${read1}"
    export fastq2="${read2}"
    export longreads="${long_read}"
    export S3OUTPUTPATH="${output_path}/${sample}"
    run_unicycler_hybrid_nanopore_2.sh
    """
}

/*
 * Run Hybridassembly Pipeline Step 3
 * Post-processing
*/

process hybridassembly_3 {
    tag "$sample"
    container "fischbachlab/nf-hybridassembly:latest"
    cpus 16
    memory 64.GB

    //publishDir "${output_path}", mode:'copy'

    input:
    tuple val(sample), path(read_1), path(read_2), path(long_read), path(assembled_reads), path(assembly_log) from reads_ch2


    output:
    //path "*"

    script:
    """
    export fastq1="${read_1}"
    export fastq2="${read_2}"
    export longreads="${long_read}"
    export assembly="${assembled_reads}"
    export log="${assembly_log}"
    export S3OUTPUTPATH="${output_path}/${sample}"
    run_unicycler_hybrid_nanopore_3.sh
    """
}

#!/bin/bash -x

set -e
set -u
set -o pipefail

START_TIME=$SECONDS
export PATH="/opt/conda/bin:${PATH}"

LOCAL=$(pwd)
coreNum=${coreNum:-16}



# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
QC_FASTQ="${OUTPUTDIR}/trimmed_fastq"


LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
ASSEMBLY_OUTPUT="${LOCAL_OUTPUT}/UNICYCLER"
QUAST_OUTPUT="${LOCAL_OUTPUT}/quast"
FASTQC_OUTPUT="${LOCAL_OUTPUT}/fastqc_illumina"
FASTQC_OUTPUT2="${LOCAL_OUTPUT}/fastqc_nanopore"
RAW_NANOPORE="${LOCAL_OUTPUT}/nanopore"
#FASTQ_NAME=${fastq1%/*}
SAMPLE_NAME=$(basename ${S3OUTPUTPATH})

LOCAL_DB_PATH=${LOCAL}/databases
LOCAL_DB_NAME="contigs"
OUTPUT_PREFIX="trimmed_reads_vs_${LOCAL_DB_NAME}"
BAMQC_OUTPUT="${LOCAL_OUTPUT}/bamqc_reads_vs_contigs"
BAMQC_OUTPUT2="${LOCAL_OUTPUT}/bamqc_longreads_vs_contigs"
BWT_OUTPUT="${LOCAL_OUTPUT}/bowtie2"
POST_OUTPUT="${LOCAL_OUTPUT}/post-processing"
ANNO_OUTPUT="${LOCAL_OUTPUT}/annotation"

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}" "${QC_FASTQ}" "${POST_OUTPUT}"
mkdir -p "${ASSEMBLY_OUTPUT}" "${QUAST_OUTPUT}" "${FASTQC_OUTPUT}" "${FASTQC_OUTPUT2}" "${ANNO_OUTPUT}"
mkdir -p "${LOCAL_DB_PATH}" "${BWT_OUTPUT}"  "${BAMQC_OUTPUT}" "${BAMQC_OUTPUT2}" "${RAW_NANOPORE}"
trap '{ rm -rf ${OUTPUTDIR} ; exit 255; }' 1

hash_kmer=${hash_kmer:-51}


mv $fastq1 "${QC_FASTQ}/read1_sampled.fastq.gz"
mv $fastq2 "${QC_FASTQ}/read2_sampled.fastq.gz"
mv $longreads "${QC_FASTQ}/long_trimmed.fastq.gz"

## Run assembly
# --contamination lambda nanopore only
unicycler \
-t ${coreNum} \
--contamination lambda \
--min_fasta_length 1000 \
--mode bold \
-1 "${QC_FASTQ}/read1_sampled.fastq.gz" \
-2 "${QC_FASTQ}/read2_sampled.fastq.gz" \
-l "${QC_FASTQ}/long_trimmed.fastq.gz" \
-o ${ASSEMBLY_OUTPUT} \
&>> ${LOG_DIR}/unicycler_assembly.log.txt


######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'Unicycler took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
#aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
# rm -rf "${OUTPUTDIR}"

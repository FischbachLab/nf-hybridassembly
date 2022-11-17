#!/bin/bash -x

set -e
set -u
set -o pipefail

START_TIME=$SECONDS
export PATH="/opt/conda/bin:${PATH}"

LOCAL=$(pwd)
coreNum=${coreNum:-16};
sampleRate=${sampleRate:-100}
LongTargetBase=${LongTargetBase:-1000000000}
# s3 inputs from env variables
#longreads="${1}" input long reads in fastq
#S3OUTPUTPATH = "${2}"

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
QC_FASTQ="${OUTPUTDIR}/trimmed_fastq"

LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
ASSEMBLY_OUTPUT="${LOCAL_OUTPUT}/UNICYCLER"
QUAST_OUTPUT="${LOCAL_OUTPUT}/quast"
FASTQC_OUTPUT="${LOCAL_OUTPUT}/fastqc"
FASTQC_OUTPUT2="${LOCAL_OUTPUT}/fastqc_pacbio"


mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}" "${QC_FASTQ}"
mkdir -p "${ASSEMBLY_OUTPUT}" "${QUAST_OUTPUT}" "${FASTQC_OUTPUT}" "${FASTQC_OUTPUT2}"
trap '{ rm -rf ${OUTPUTDIR} ; exit 255; }' 1

# Copy fastq.gz files from S3, only 2 files per sample
aws s3 cp --quiet ${longreads} "${RAW_FASTQ}/long.fastq.gz"

# convert bam to fastq for long reads
#samtools bam2fq  "${RAW_FASTQ}/long.bam" >  "${RAW_FASTQ}/long.fastq"

# filter long reads without reference --target_bases 500000000 \  --keep_percent $sampleRate \
filtlong \
--min_length 1000 \
--length_weight 10 \
--target_bases ${LongTargetBase} \
"${RAW_FASTQ}/long.fastq.gz" | gzip > "${QC_FASTQ}/long_trimmed.fastq.gz"


#Run fastqc for long reads
#fastqc \
#-t ${coreNum} \
#-o ${FASTQC_OUTPUT2} \
#"${QC_FASTQ}/long_trimmed.fastq.gz"


## Run assembly
timem unicycler \
-t ${coreNum} \
--mode bold \
-l "${QC_FASTQ}/long_trimmed.fastq.gz" \
-o ${ASSEMBLY_OUTPUT} |\
tee -a ${LOG_DIR}/unicycler_assembly.log.txt

# Run Quast without reference
quast.py \
--glimmer \
-t ${coreNum} \
${ASSEMBLY_OUTPUT}/assembly.fasta \
-o ${QUAST_OUTPUT} | tee -a ${LOG_DIR}/quast.log.txt

#--pacbio "${QC_FASTQ}/long_trimmed.fastq.gz" \
#--glimmer \
#--rna-finding \
######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
# rm -rf "${OUTPUTDIR}"

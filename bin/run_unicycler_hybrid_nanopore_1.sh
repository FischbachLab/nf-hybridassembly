#!/bin/bash -x

set -e
set -u
set -o pipefail

START_TIME=$SECONDS
export PATH="/opt/conda/bin:${PATH}"

LOCAL=$(pwd)
coreNum=${coreNum:-32}
ShortTargetBase=${ShortTargetBase:-500000000}
LongTargetBase=${LongTargetBase:-500000000}
genomeSize=${genomeSize:-5000000}
coverage=${coverage:-400}

# if the genomeSize is give, then overwrite the TargetBase parameters
ShortTargetBase=$((genomeSize*coverage*1))
LongTargetBase=$((genomeSize*coverage*1))

# Make sure neither the short read set nor long read set exceeds the certain limit
if [ $ShortTargetBase -gt 3990000000 ]
then
   ShortTargetBase=3990000000
fi

if [ $LongTargetBase -gt 3990000000 ]
then
   LongTargetBase=3990000000
fi
# s3 inputs from env variables
#fastq1="${1}"
#fastq2="${2}"
#longreads="${3}" input long reads in fastq format
#S3OUTPUTPATH="${4}"
#genomeSize="${5}" genome size in bp


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

# Copy fastq.gz files from S3, only 2 files per sample
aws s3 cp --quiet ${fastq1} "${RAW_FASTQ}/read1.fastq.gz"
aws s3 cp --quiet ${fastq2} "${RAW_FASTQ}/read2.fastq.gz"
aws s3 cp --quiet ${longreads} "${RAW_FASTQ}/long.fastq.gz"

# count the number of PE reads
#zcat "${RAW_FASTQ}/read1.fastq.gz" | echo $((`wc -l`/2)) > ${LOG_DIR}/illumina_count.txt
# count the long reads
#zcat "${RAW_FASTQ}/long.fastq" | echo $((`wc -l`/4)) > ${LOG_DIR}/nanopore_count.txt


## Collate called bases into a single fastq file
#aws s3 sync --quiet ${longreads} ${RAW_NANOPORE}
#find "${RAW_NANOPORE}/*.fastq" | xargs cat > "${RAW_FASTQ}/long.fastq"

#cp  ${fastq1} "${RAW_FASTQ}/read1.fastq.gz"
#cp  ${fastq2} "${RAW_FASTQ}/read2.fastq.gz"


###############################################################
# Pre-processing before assembly
echo "Reads Pre-processing before assembly"
###############################################################

# discard reads that have mismatching lengths of bases and qualities
echo "**************************" >> ${LOG_DIR}/bbtools.log.txt
echo "Reads Reformatting" >> ${LOG_DIR}/bbtools.log.txt
echo "**************************" >> ${LOG_DIR}/bbtools.log.txt
reformat.sh -Xmx16g -eoom \
  in="${RAW_FASTQ}/read1.fastq.gz" \
  in2="${RAW_FASTQ}/read2.fastq.gz" \
  out="${QC_FASTQ}/repaired-interleaved.fastq.gz" \
  tossbrokenreads=t &> ${LOG_DIR}/bbtools.log.txt

#repair.sh \
# in="${RAW_FASTQ}/read1.fastq.gz" \
# in2="${RAW_FASTQ}/read2.fastq.gz" \
# out="${QC_FASTQ}/repaired-interleaved.fastq.gz" |\
# tee -a ${LOG_DIR}/bbtools.log.txt

# Use dedupe to remove the exact PCR duplicates
#in="${RAW_FASTQ}/read1.fastq.gz" \
#in2="${RAW_FASTQ}/read2.fastq.gz" \
#dedupe.sh -Xmx120g \
#    in="${QC_FASTQ}/repaired-interleaved.fastq.gz" \
#    out="${QC_FASTQ}/deduped-interleaved.fastq.gz" \
#    tossbrokenreads=t &> ${LOG_DIR}/bbtools.log.txt


#reformat.sh \
#    in="${QC_FASTQ}/deduped-interleaved.fastq.gz" \
#    out1="${QC_FASTQ}/read1_deduped.fastq.gz" \
#    out2="${QC_FASTQ}/read2_deduped.fastq.gz"
# repair the PE reads

#repair.sh \
# in="${QC_FASTQ}/deduped-interleaved.fastq.gz" \
# out="${QC_FASTQ}/repaired-interleaved.fastq.gz" |\
# tee -a ${LOG_DIR}/bbtools.log.txt

# Constant definitions for bbduk, increased the quailty
adapterFile="adapters,phix"
trimQuality=${trimQuality:-25}
minLength=${minLength:-50}
kmer_value=${kmer_value:-23}
min_kmer_value=${min_kmer_value:-11}

# Use bbduk to trim short reads, -eoom exits when out of memory
#in1="${QC_FASTQ}/read1_deduped.fastq.gz" \
#in2="${QC_FASTQ}/read2_deduped.fastq.gz" \
echo "**************************" >> ${LOG_DIR}/bbtools.log.txt
echo "Reads trimming/filtering" >> ${LOG_DIR}/bbtools.log.txt
echo "**************************" >> ${LOG_DIR}/bbtools.log.txt
bbduk.sh -Xmx16g tbo -eoom hdist=1 qtrim=rl ktrim=r \
    entropy=0.5 entropywindow=50 entropyk=5 entropytrim=rl \
    in="${QC_FASTQ}/repaired-interleaved.fastq.gz" \
    out1="${QC_FASTQ}/read1_trimmed.fastq.gz" \
    out2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
    ref=${adapterFile} \
    k="${kmer_value}" \
    mink="${min_kmer_value}" \
    trimq="${trimQuality}" \
    minlen="${minLength}" \
    tossbrokenreads=t \
    refstats="${LOCAL_OUTPUT}/BBDuk/adapter_trimming_stats_per_ref.txt" \
    >> ${LOG_DIR}/bbtools.log.txt 2>&1


reformat.sh -Xmx16g -eoom \
      in="${QC_FASTQ}/read1_trimmed.fastq.gz" \
      in2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
      out="${QC_FASTQ}/trimmed-interleaved.fastq.gz" \
      >> ${LOG_DIR}/bbtools.log.txt 2>&1

echo "**************************" >> ${LOG_DIR}/bbtools.log.txt
echo "Dereplication" >> ${LOG_DIR}/bbtools.log.txt
echo "**************************" >> ${LOG_DIR}/bbtools.log.txt
timem dedupe.sh -Xmx100g \
        in="${QC_FASTQ}/trimmed-interleaved.fastq.gz" \
        out="${QC_FASTQ}/deduped-interleaved.fastq.gz" \
        >> ${LOG_DIR}/bbtools.log.txt 2>&1


# Normalize uneven coverage, depth of under 3x will be presumed to be errors and discarded.
# target is kmer detpth, not read depth ( read depth = kmer depth * (read length/(read length -kmer size +1 )) )
#in="${QC_FASTQ}/read1_trimmed.fastq.gz" \
#in2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
#bbnorm.sh in=reads.fq out=corrected.fq ecc=t keepall passes=1 bits=16 prefilter
echo "**************************" >> ${LOG_DIR}/bbtools.log.txt
echo "Coverage Normalization" >> ${LOG_DIR}/bbtools.log.txt
echo "**************************" >> ${LOG_DIR}/bbtools.log.txt
timem bbnorm.sh -Xmx100g \
target=${coverage} min=3 \
in="${QC_FASTQ}/deduped-interleaved.fastq.gz" \
out="${QC_FASTQ}/normalized-interleaved.fastq.gz" \
>> ${LOG_DIR}/bbtools.log.txt 2>&1

#Uses kmer counts to error-correct reads.
echo "**************************" >> ${LOG_DIR}/bbtools.log.txt
echo "Error correction" >> ${LOG_DIR}/bbtools.log.txt
echo "**************************" >> ${LOG_DIR}/bbtools.log.txt
timem tadpole.sh -Xmx30g \
mode=correct \
in="${QC_FASTQ}/normalized-interleaved.fastq.gz" \
out="${QC_FASTQ}/normalized-correct-interleaved.fastq.gz" \
>> ${LOG_DIR}/bbtools.log.txt 2>&1

echo "**************************" >> ${LOG_DIR}/bbtools.log.txt
echo "Reads Sampling" >> ${LOG_DIR}/bbtools.log.txt
echo "**************************" >> ${LOG_DIR}/bbtools.log.txt
# Downsample the short reads if Uncycler hang ${sampleRate}
# Downsample if the normalized reads are still over the upper limit
#in="${QC_FASTQ}/read1_trimmed.fastq.gz" \
#in2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
reformat.sh \
samplebasestarget=${ShortTargetBase} \
in="${QC_FASTQ}/normalized-correct-interleaved.fastq.gz" \
out="${QC_FASTQ}/read1_sampled.fastq.gz" \
out2="${QC_FASTQ}/read2_sampled.fastq.gz" \
>> ${LOG_DIR}/bbtools.log.txt 2>&1

#Run fastqc for short reads
fastqc \
-t ${coreNum} \
-o ${FASTQC_OUTPUT} \
"${QC_FASTQ}/read1_sampled.fastq.gz" \
"${QC_FASTQ}/read2_sampled.fastq.gz"


echo "**************************"
echo "Pre-Processing Long reads"
echo "**************************"


# convert bam to fastq for long reads
#samtools bam2fq  "${RAW_FASTQ}/long.bam" >  "${RAW_FASTQ}/long.fastq"
# filter long reads with reference(short reads) --trim \
filtlong \
-1 "${QC_FASTQ}/read1_trimmed.fastq.gz" \
-2 "${QC_FASTQ}/read2_trimmed.fastq.gz" \
--min_length 1000 \
--keep_percent 90 \
--length_weight 10 \
--target_bases ${LongTargetBase} \
"${RAW_FASTQ}/long.fastq.gz" | gzip > "${QC_FASTQ}/long_trimmed.fastq.gz"
# comment out --split option to keep the reads as long as possible in case some regions are not covered by short reads
#--split 500 \

#Run fastqc for long reads
fastqc \
-t ${coreNum} \
-o ${FASTQC_OUTPUT2} \
"${QC_FASTQ}/long_trimmed.fastq.gz"



######################### HOUSEKEEPING #############################
#DURATION=$((SECONDS - START_TIME))
#hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
#printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
#echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
# rm -rf "${OUTPUTDIR}"

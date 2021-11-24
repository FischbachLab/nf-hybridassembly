#!/bin/bash -x

set -e
set -u
set -o pipefail

START_TIME=$SECONDS
export PATH="/opt/conda/bin:${PATH}"

LOCAL=$(pwd)
coreNum=${coreNum:-16};
genomeSize=${genomeSize:-5000000}
coverage=${coverage:-300}

# if a genomeSize is give, then overwrite the TargetBase parameters
ShortTargetBase=$((genomeSize*coverage))
LongTargetBase=$((genomeSize*coverage))

# s3 inputs from env variables
#fastq1="${1}"
#fastq2="${2}"
#longreads="${3}" input Nanopore long reads in bam format
#assembly={4} input existing_long_read_assembly reads in fasta
#S3OUTPUTPATH = "${5}"

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
FASTQ_NAME=${fastq1%/*}
SAMPLE_NAME=$(basename ${FASTQ_NAME})

LOCAL_DB_PATH=${LOCAL}/databases
LOCAL_DB_NAME="contigs"
OUTPUT_PREFIX="trimmed_reads_vs_${LOCAL_DB_NAME}"
BAMQC_OUTPUT="${LOCAL_OUTPUT}/bamqc_reads_vs_contigs"
BWT_OUTPUT="${LOCAL_OUTPUT}/bowtie2"

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}" "${QC_FASTQ}"
mkdir -p "${ASSEMBLY_OUTPUT}" "${QUAST_OUTPUT}" "${FASTQC_OUTPUT}" "${FASTQC_OUTPUT2}"
mkdir -p "${LOCAL_DB_PATH}" "${BWT_OUTPUT}"  "${BAMQC_OUTPUT}" "${BWT_OUTPUT}"

trap '{ rm -rf ${OUTPUTDIR} ; exit 255; }' 1

hash_kmer=${hash_kmer:-51}

# Copy fastq.gz files from S3, only 2 files per sample
aws s3 cp --quiet ${fastq1} "${RAW_FASTQ}/read1.fastq.gz"
aws s3 cp --quiet ${fastq2} "${RAW_FASTQ}/read2.fastq.gz"
aws s3 cp --quiet ${longreads} "${RAW_FASTQ}/long.fastq.gz"
aws s3 cp --quiet ${assembly} "${RAW_FASTQ}/existing_assembly.fa"

#cp  ${fastq1} "${RAW_FASTQ}/read1.fastq.gz"
#cp  ${fastq2} "${RAW_FASTQ}/read2.fastq.gz"

###############################################################
# Pre-processing before assembly
echo "Reads Pre-processing before assembly"
###############################################################

# discard reads that have mismatching lengths of bases and qualities
reformat.sh -Xmx16g -eoom \
  in="${RAW_FASTQ}/read1.fastq.gz" \
  in2="${RAW_FASTQ}/read2.fastq.gz" \
  out="${QC_FASTQ}/repaired-interleaved.fastq.gz" \
  tossbrokenreads=t &> ${LOG_DIR}/bbtools.log.txt

echo "**************************" >> ${LOG_DIR}/bbtools.log.txt
echo "Reads trimming/Filtering" >> ${LOG_DIR}/bbtools.log.txt
echo "**************************" >> ${LOG_DIR}/bbtools.log.txt

# Constant definitions for bbduk
adapterFile="adapters,phix"
trimQuality=${trimQuality:-25}
minLength=${minLength:-50}
kmer_value=${kmer_value:-23}
min_kmer_value=${min_kmer_value:-11}

# Use bbduk to trim short reads, -eoom exits when out of memory
bbduk.sh -Xmx16g tbo -eoom hdist=1 qtrim=rl ktrim=r \
    entropy=0.5 entropywindow=50 entropyk=5 \
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
timem dedupe.sh -Xmx30g \
        in="${QC_FASTQ}/trimmed-interleaved.fastq.gz" \
        out="${QC_FASTQ}/deduped-interleaved.fastq.gz" \
        >> ${LOG_DIR}/bbtools.log.txt 2>&1

echo "**************************" >> ${LOG_DIR}/bbtools.log.txt
echo "Coverage Normalization" >> ${LOG_DIR}/bbtools.log.txt
echo "**************************" >> ${LOG_DIR}/bbtools.log.txt
timem bbnorm.sh -Xmx30g \
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

# Downsample the short reads if Uncycler hang ${sampleRate}
# Downsample if the normalized reads are still over the upper limit
#in="${QC_FASTQ}/read1_trimmed.fastq.gz" \
#in2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
reformat.sh \
samplebasestarget=${ShortTargetBase} \
in="${QC_FASTQ}/normalized-correct-interleaved.fastq.gz" \
out="${QC_FASTQ}/read1_sampled.fastq.gz" \
out2="${QC_FASTQ}/read2_sampled.fastq.gz" \
&>> ${LOG_DIR}/bbtools.log.txt

#Run fastqc for short reads
fastqc \
-t ${coreNum} \
-o ${FASTQC_OUTPUT} \
"${QC_FASTQ}/read1_sampled.fastq.gz" \
"${QC_FASTQ}/read2_sampled.fastq.gz"


# convert bam to fastq for long reads
#samtools bam2fq  "${RAW_FASTQ}/long.bam" >  "${RAW_FASTQ}/long.fastq"

# filter long reads with reference --target_bases 500000000 \ #--trim \
#--split 500 \
filtlong \
-1 "${QC_FASTQ}/read1_trimmed.fastq.gz" \
-2 "${QC_FASTQ}/read2_trimmed.fastq.gz" \
--min_length 1000 \
--keep_percent 90 \
--target_bases ${LongTargetBase} \
"${RAW_FASTQ}/long.fastq.gz" | gzip > "${QC_FASTQ}/long_trimmed.fastq.gz"


# Downsample the short reads if Uncycler hang ${sampleRate}
#reformat.sh \
#samplebasestarget=${ShortTargetBase} \
#in="${QC_FASTQ}/read1_trimmed.fastq.gz" \
#in2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
#out="${QC_FASTQ}/read1_sampled.fastq.gz" \
#out2="${QC_FASTQ}/read2_sampled.fastq.gz" \
#&>> ${LOG_DIR}/bbtools.log.txt

#Run fastqc for long reads
#fastqc \
#-t ${coreNum} \
#-o ${FASTQC_OUTPUT2} \
#"${QC_FASTQ}/long_trimmed.fastq.gz"


## Run assembly --existing_long_read_assembly must be used with hybrid mode only\
timem unicycler \
-t ${coreNum} \
-1 "${QC_FASTQ}/read1_sampled.fastq.gz" \
-2 "${QC_FASTQ}/read2_sampled.fastq.gz" \
-l "${QC_FASTQ}/long_trimmed.fastq.gz" \
--mode bold \
--existing_long_read_assembly "${RAW_FASTQ}/existing_assembly.fa" \
-o ${ASSEMBLY_OUTPUT} \
&>> ${LOG_DIR}/unicycler_assembly.log.txt


## Build the database
# --threads ${coreNum} \
bowtie2-build --seed 42 \
    ${ASSEMBLY_OUTPUT}/assembly.fasta \
    ${LOCAL_DB_PATH}/${LOCAL_DB_NAME} \
    &>> ${LOG_DIR}/bowtie2_build_db_index.log

bowtie2 --sensitive-local -p ${coreNum} --seed 42 -x ${LOCAL_DB_PATH}/${LOCAL_DB_NAME} \
-1 "${QC_FASTQ}/read1_sampled.fastq.gz" -2 "${QC_FASTQ}/read2_sampled.fastq.gz" -S "${BWT_OUTPUT}/contigs_aligned.sam"
# index contigs
samtools faidx "${ASSEMBLY_OUTPUT}/assembly.fasta"
# create bam file
samtools import "${ASSEMBLY_OUTPUT}/assembly.fasta.fai" "${BWT_OUTPUT}/contigs_aligned.sam" "${BWT_OUTPUT}/contigs_aligned.bam"
# sort bam file
samtools sort -@ ${coreNum} "${BWT_OUTPUT}/contigs_aligned.bam"  -o "${BWT_OUTPUT}/contigs_aligned.sorted.bam"
# index bam
samtools index  "${BWT_OUTPUT}/contigs_aligned.sorted.bam"
# generate contigs stats
samtools idxstats "${BWT_OUTPUT}/contigs_aligned.sorted.bam" > "${LOG_DIR}/contigs_aligned.idxstats.txt"
samtools flagstat "${BWT_OUTPUT}/contigs_aligned.sorted.bam" > "${LOG_DIR}/contigs_aligned.flagstat.txt"

#bamqc
qualimap bamqc --java-mem-size=16G -nt ${coreNum} -outdir "${BAMQC_OUTPUT}" -bam "${BWT_OUTPUT}/contigs_aligned.sorted.bam" -c
# clean data
rm "${BWT_OUTPUT}/contigs_aligned.sam" "${BWT_OUTPUT}/contigs_aligned.bam" "${BWT_OUTPUT}/contigs_aligned.sorted.bam"


# Run Quast without reference
quast.py \
-t ${coreNum} \
--glimmer \
--rna-finding \
--contig-thresholds 5000,10000,25000,50000,100000,250000,500000,1000000 \
${ASSEMBLY_OUTPUT}/assembly.fasta \
-o ${QUAST_OUTPUT} &>> ${LOG_DIR}/quast.log.txt

#########################################################
echo "Reads stats after assembly"
#########################################################

# Count input reads, PE reads count once
totalShortReads=$(( $( zcat ${RAW_FASTQ}/read1.fastq.gz | wc -l ) / 4 ))
shortAfterTrim=$(( $( zcat  ${QC_FASTQ}/read1_sampled.fastq.gz | wc -l ) / 4 ))
totalContigs=$(grep -c "^>" ${ASSEMBLY_OUTPUT}/assembly.fasta )
shortUsedRate=`echo "scale=2; $shortAfterTrim*100/$totalShortReads" | bc -l`

if grep -q "  incomplete" ${ASSEMBLY_OUTPUT}/unicycler.log; then
    status="Incomplete"
else
    status="Complete"
fi

short_mapping=`grep "number of mapped reads" "${BAMQC_OUTPUT}/genome_results.txt" | awk -F '[()]' '{print $2}'`


echo 'Sample_Name,Total_ShortReads,Total_contigs,shortAfterFiltering,short_mapping,Assembly_status' > ${LOG_DIR}/read_stats.csv
echo ${SAMPLE_NAME}','${totalShortReads}','${totalContigs}','${shortAfterTrim}'('$shortUsedRate'%),'${short_mapping}','${status} >> ${LOG_DIR}/read_stats.csv

# Grep Assembly summary
tail -n 1 ${LOG_DIR}/unicycler_assembly.log.txt > ${LOG_DIR}/unicycler_summary.txt
grep -A 15 "^Bridged assembly graph" ${ASSEMBLY_OUTPUT}/unicycler.log >> ${LOG_DIR}/unicycler_summary.txt




######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
# rm -rf "${OUTPUTDIR}"

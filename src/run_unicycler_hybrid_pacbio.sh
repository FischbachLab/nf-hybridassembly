#!/bin/bash -x

set -e
set -u
set -o pipefail

START_TIME=$SECONDS
export PATH="/opt/conda/bin:${PATH}"

LOCAL=$(pwd)
coreNum=${coreNum:-16}
ShortTargetBase=${ShortTargetBase:-500000000}
LongTargetBase=${LongTargetBase:-500000000}
genomeSize=${genomeSize:-5000000}
coverage=${coverage:-400}

# if the genomeSize is give, then overwrite the TargetBase parameters
ShortTargetBase=$((genomeSize*coverage*1))
LongTargetBase=$((genomeSize*coverage*1))

# Make sure neither the short read set nor long read set exceeds 2000 Mbp
if [ $ShortTargetBase -gt 2000000000 ]
then
   ShortTargetBase=2000000000
fi

# Keep more long reads
if [ $LongTargetBase -gt 2000000000 ]
then
   LongTargetBase=2000000000
fi
# s3 inputs from env variables
#fastq1="${1}"
#fastq2="${2}"
#longreads="${3}" input long reads in fastq format
#S3OUTPUTPATH="${4}"
#genomeSize="${5}" genome size in bp

# for inputing existing long reads
#assembly={6} input existing_long_read_assembly reads in fasta

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
QC_FASTQ="${OUTPUTDIR}/trimmed_fastq"


LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
ASSEMBLY_OUTPUT="${LOCAL_OUTPUT}/UNICYCLER"
QUAST_OUTPUT="${LOCAL_OUTPUT}/quast"
FASTQC_OUTPUT="${LOCAL_OUTPUT}/fastqc_illumina"
FASTQC_OUTPUT2="${LOCAL_OUTPUT}/fastqc_pacbio"
RAW_PACBIO="${LOCAL_OUTPUT}/pacbio"
#FASTQ_NAME=${fastq1%/*}
SAMPLE_NAME=$(basename ${S3OUTPUTPATH})

LOCAL_DB_PATH=${LOCAL}/databases
LOCAL_DB_NAME="contigs"
OUTPUT_PREFIX="trimmed_reads_vs_${LOCAL_DB_NAME}"
BAMQC_OUTPUT="${LOCAL_OUTPUT}/bamqc_reads_vs_contigs"
BAMQC_OUTPUT2="${LOCAL_OUTPUT}/bamqc_longreads_vs_contigs"
BWT_OUTPUT="${LOCAL_OUTPUT}/bowtie2"
POST_OUTPUT="${LOCAL_OUTPUT}/post-processing"

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}" "${QC_FASTQ}" "${POST_OUTPUT}"
mkdir -p "${ASSEMBLY_OUTPUT}" "${QUAST_OUTPUT}" "${FASTQC_OUTPUT}" "${FASTQC_OUTPUT2}"
mkdir -p "${LOCAL_DB_PATH}" "${BWT_OUTPUT}"  "${BAMQC_OUTPUT}" "${BAMQC_OUTPUT2}" "${RAW_PACBIO}"
trap '{ rm -rf ${OUTPUTDIR} ; exit 255; }' 1

hash_kmer=${hash_kmer:-51}

# Copy fastq.gz files from S3, only 2 files per sample
aws s3 cp --quiet ${fastq1} "${RAW_FASTQ}/read1.fastq.gz"
aws s3 cp --quiet ${fastq2} "${RAW_FASTQ}/read2.fastq.gz"
#aws s3 cp --quiet ${longreads} "${RAW_FASTQ}/long.fastq.gz"

# If long reads in bam format
aws s3 cp --quiet ${longreads} "${RAW_FASTQ}/long.bam"
samtools bam2fq  "${RAW_FASTQ}/long.bam" | gzip > "${RAW_FASTQ}/long.fastq.gz"
# copy existing long reads
#aws s3 cp --quiet ${assembly} "${RAW_FASTQ}/existing_assembly.fa"

## Collate called bases into a single fastq file
#aws s3 sync --quiet ${longreads} ${RAW_PACBIO}
#find "${RAW_PACBIO}/*.fastq" | xargs cat > "${RAW_FASTQ}/long.fastq"

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
trimQuality=${trimQuality:-25} #old 25
minLength=${minLength:-50}  #old 50
kmer_value=${kmer_value:-23}
min_kmer_value=${min_kmer_value:-11}

# Use bbduk to trim short reads, -eoom exits when out of memory
#in1="${QC_FASTQ}/read1_deduped.fastq.gz" \
#in2="${QC_FASTQ}/read2_deduped.fastq.gz" \
echo "**************************" >> ${LOG_DIR}/bbtools.log.txt
echo "Reads trimming/Filtering" >> ${LOG_DIR}/bbtools.log.txt
echo "**************************" >> ${LOG_DIR}/bbtools.log.txt
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

# convert bam to fastq for long reads
#samtools bam2fq  "${RAW_FASTQ}/long.bam" >  "${RAW_FASTQ}/long.fastq"
#Uses kmer counts to error-correct reads.
echo "**************************"
echo "Long Reads Filtering"
echo "**************************"

# get long read length Histogram
readlength.sh in="${RAW_FASTQ}/long.fastq.gz" bin=1000 max=10000 ignorebadquality > ${LOG_DIR}/longreads.LengthHistogram.txt

# filter long reads with reference(short reads) --trim \ --target_bases ${LongTargetBase} \ --length_weight 10 \
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

# get filtered long read length Histogram
readlength.sh in="${QC_FASTQ}/long_trimmed.fastq.gz" bin=1000 max=20000 ignorebadquality > ${LOG_DIR}/Filtered_longreads.LengthHistogram.txt

## Run assembly --mode bold \
# --contamination lambda nanopore only --existing_long_read_assembly "${RAW_FASTQ}/existing_assembly.fa" \
timem unicycler \
-t ${coreNum} \
--min_fasta_length 1000 \
-1 "${QC_FASTQ}/read1_sampled.fastq.gz" \
-2 "${QC_FASTQ}/read2_sampled.fastq.gz" \
-l "${QC_FASTQ}/long_trimmed.fastq.gz" \
-o ${ASSEMBLY_OUTPUT} |\
tee -a ${LOG_DIR}/unicycler_assembly.log.txt


## Build the database
# --threads ${coreNum} \
bowtie2-build --seed 42 \
    ${ASSEMBLY_OUTPUT}/assembly.fasta \
    ${LOCAL_DB_PATH}/${LOCAL_DB_NAME} |\
    tee -a ${LOG_DIR}/bowtie2_build_db_index.log

bowtie2 --sensitive-local -p ${coreNum} --seed 42 -x ${LOCAL_DB_PATH}/${LOCAL_DB_NAME} \
-1 "${QC_FASTQ}/read1_sampled.fastq.gz" -2 "${QC_FASTQ}/read2_sampled.fastq.gz" |\
samtools view -@ ${coreNum} -bh -o "${BWT_OUTPUT}/contigs_aligned.bam" - | \
tee -a ${LOG_DIR}/bowtie2_mapping.log.txt

# -S "${BWT_OUTPUT}/contigs_aligned.sam"
# index contigs
#samtools faidx "${ASSEMBLY_OUTPUT}/assembly.fasta"
# create bam file
#samtools import "${ASSEMBLY_OUTPUT}/assembly.fasta.fai" "${BWT_OUTPUT}/contigs_aligned.sam" "${BWT_OUTPUT}/contigs_aligned.bam"
# sort bam file
samtools sort -@ ${coreNum} "${BWT_OUTPUT}/contigs_aligned.bam"  -o "${BWT_OUTPUT}/contigs_aligned.sorted.bam"
# index bam
samtools index  "${BWT_OUTPUT}/contigs_aligned.sorted.bam"
# generate contigs stats
samtools idxstats "${BWT_OUTPUT}/contigs_aligned.sorted.bam" > "${LOG_DIR}/contigs_aligned.idxstats.txt"
samtools flagstat "${BWT_OUTPUT}/contigs_aligned.sorted.bam" > "${LOG_DIR}/contigs_aligned.flagstat.txt"

#bamqc
qualimap bamqc --java-mem-size=16G -nt ${coreNum} -outdir "${BAMQC_OUTPUT}" -bam "${BWT_OUTPUT}/contigs_aligned.sorted.bam" -c

# align long reads to contigs
minimap2 -t ${coreNum} -a "${ASSEMBLY_OUTPUT}/assembly.fasta"  "${QC_FASTQ}/long_trimmed.fastq.gz" | samtools sort -@${coreNum} -o "${BAMQC_OUTPUT2}/longreads_vs_contigs.sorted.bam" -
samtools index "${BAMQC_OUTPUT2}/longreads_vs_contigs.sorted.bam"
# bam qc
qualimap bamqc --java-mem-size=16G -nt ${coreNum} -outdir "${BAMQC_OUTPUT2}" -bam "${BAMQC_OUTPUT2}/longreads_vs_contigs.sorted.bam" -c


# Run Quast without reference
quast.py \
-t ${coreNum} \
--glimmer \
--rna-finding \
--contig-thresholds 5000,10000,25000,50000,100000,250000,500000,1000000 \
${ASSEMBLY_OUTPUT}/assembly.fasta \
-o ${QUAST_OUTPUT} | tee -a ${LOG_DIR}/quast.log.txt

#########################################################
echo "Reads stats after assembly"
#########################################################

# Count input reads, PE reads count once
totalShortReads=$(( $( zcat ${RAW_FASTQ}/read1.fastq.gz | wc -l ) / 4 ))
totalTrimmedReads=$(( $( zcat ${QC_FASTQ}/read1_trimmed.fastq.gz | wc -l ) / 4 ))
totalLongReads=$(( $( zcat ${RAW_FASTQ}/long.fastq.gz | wc -l ) / 4 ))
sampledReads=$(( $( zcat  ${QC_FASTQ}/read1_sampled.fastq.gz | wc -l ) / 4 ))
longAfterTrim=$(( $( zcat ${QC_FASTQ}/long_trimmed.fastq.gz | wc -l ) / 4 ))
totalContigs=$(grep -c "^>" ${ASSEMBLY_OUTPUT}/assembly.fasta )

shortUsedRate=`echo "scale=2; $sampledReads*100/$totalTrimmedReads" | bc -l`
longUsedRate=`echo "scale=2; $longAfterTrim*100/$totalLongReads" | bc -l`

AvgLongReadLength=`cat ${LOG_DIR}/longreads.LengthHistogram.txt | grep -P "^#Avg:" | cut -f 2-`

if grep -q "  incomplete" ${ASSEMBLY_OUTPUT}/unicycler.log; then
    status="Incomplete"
else
    status="Complete"
fi

short_mapping=`grep "number of mapped reads" "${BAMQC_OUTPUT}/genome_results.txt" | awk -F '[()]' '{print $2}'`
long_mapping=`grep "number of mapped reads" "${BAMQC_OUTPUT2}/genome_results.txt" | awk -F '[()]' '{print $2}'`

echo "*********Get long reads coverage stats by short reads *****************"
minimap2 -t ${coreNum} -ax sr "${QC_FASTQ}/long_trimmed.fastq.gz" "${QC_FASTQ}/read1_sampled.fastq.gz" "${QC_FASTQ}/read2_sampled.fastq.gz" > "${BWT_OUTPUT}/short_vs_long.sam"
pileup.sh in="${BWT_OUTPUT}/short_vs_long.sam" out=${LOG_DIR}/long_vs_short-coverage.txt 2> ${LOG_DIR}/long_vs_short-stats.txt

LongByShortCoverage=`cat ${LOG_DIR}/long_vs_short-stats.txt | grep -P "^Percent of reference bases covered:" | cut -f 2-`

# clean data
rm -r "${BWT_OUTPUT}"

echo 'Sample_Name,Total_LongReads,Total_ShortReads,Total_contigs,longAfterSampling,shortAfterSampling,long_mapping,short_mapping,AvgLongReadLength,LongByShortCoverage,Assembly_status' > ${LOG_DIR}/read_stats.csv
echo ${SAMPLE_NAME}','${totalLongReads}','${totalShortReads}','${totalContigs}','${longAfterTrim}'('$longUsedRate'%),'${sampledReads}'('$shortUsedRate'%),'${long_mapping}','${short_mapping}','${AvgLongReadLength}','${LongByShortCoverage}'%,'${status} >> ${LOG_DIR}/read_stats.csv

# Grep Assembly summary
tail -n 1 ${LOG_DIR}/unicycler_assembly.log.txt > ${LOG_DIR}/unicycler_summary.txt
grep -A 15 "^Bridged assembly graph" ${ASSEMBLY_OUTPUT}/unicycler.log >> ${LOG_DIR}/unicycler_summary.txt



###############################################################
# Post-processing if the assembly is incomplete
###############################################################


if grep -q "  incomplete" ${ASSEMBLY_OUTPUT}/unicycler.log; then

    echo "**************************" >> ${LOG_DIR}/post-processing.log.txt
    echo "Assembly Post-Processing" >> ${LOG_DIR}/post-processing.log.txt
    echo "**************************" >> ${LOG_DIR}/post-processing.log.txt

      Assembly=${ASSEMBLY_OUTPUT}/assembly.fasta
      pacbio_reads="${RAW_FASTQ}/long.fastq.gz"

      genome=${SAMPLE_NAME} #"Blautia-sp-KLE"

      # Align Long reads to drraft assembly

      minimap2 -t 8 $Assembly $pacbio_reads > ${RAW_FASTQ}/aln.mm

      # Form scaffolds
      java -jar /home/gitdir/lrscaf/target/LRScaf-1.1.9.jar -c $Assembly -a ${RAW_FASTQ}/aln.mm -t mm -o ${POST_OUTPUT}/lrscaf_output

      #Convert fastq to fasta
      reformat.sh qin=33 overwrite=t ignorebadquality=t in=$pacbio_reads out=${RAW_FASTQ}/${SAMPLE_NAME}.fasta



      if [ -f "${POST_OUTPUT}/lrscaf_output/scaffolds_summary.info" ]; then
          echo "scaffolds_summary.info" >>  ${LOG_DIR}/post-processing.log.txt
          cat ${POST_OUTPUT}/lrscaf_output/scaffolds_summary.info >> ${LOG_DIR}/post-processing.log.txt
          echo "---------------------------------" >> ${LOG_DIR}/post-processing.log.txt

          mkdir -p ${POST_OUTPUT}/tgsgapcloser_output
          cd ${POST_OUTPUT}/tgsgapcloser_output
          #
          tgsgapcloser --tgstype pb --scaff ${POST_OUTPUT}/lrscaf_output/scaffolds.fasta --reads ${RAW_FASTQ}/${SAMPLE_NAME}.fasta --output gaps --ne >gaps.log 2>gaps.err

          echo "gaps.gap_fill_detail" >> ${LOG_DIR}/post-processing.log.txt
          cat gaps.gap_fill_detail >> ${LOG_DIR}/post-processing.log.txt
          echo " " >> ${LOG_DIR}/post-processing.log.txt

          cd ${POST_OUTPUT}

          cp  ${POST_OUTPUT}/tgsgapcloser_output/gaps.scaff_seqs ${POST_OUTPUT}/gapfilled_scaffolds.fasta
          mkdir -p ${POST_OUTPUT}/circular_output
          cd ${POST_OUTPUT}/circular_output

          # Completeness chcek  Output file: permute_scaffold.fasta
          circ-permute.pl -in ${POST_OUTPUT}/tgsgapcloser_output/gaps.scaff_seqs -split 600

          tgsgapcloser --tgstype pb --scaff ${POST_OUTPUT}/circular_output/permute_scaffold.fasta --reads ${RAW_FASTQ}/${SAMPLE_NAME}.fasta --output circular --ne >circular.log 2>circular.err
          echo "circular.gap_fill_detail" >> ${LOG_DIR}/post-processing.log.txt
          cat circular.gap_fill_detail >> ${LOG_DIR}/post-processing.log.txt

          cd ${POST_OUTPUT}
          cp ${POST_OUTPUT}/circular_output/circular.scaff_seqs ${POST_OUTPUT}/final_circular_scaffolds.fasta

          # Run Quast for comparison between assembled contigs and gap-closed contigs
          quast.py -t 8 --glimmer --rna-finding --contig-thresholds 2500,5000,10000,50000,100000,500000,1000000 \
          -l $genome,${genome}.GapsClosed  \
          -o post_quast $Assembly ${POST_OUTPUT}/tgsgapcloser_output/gaps.scaff_seqs

      else
          echo "NO scaffolds are in the assembly" >> ${LOG_DIR}/post-processing.log.txt
          echo "Circular check ONLY" >> ${LOG_DIR}/post-processing.log.txt

          mkdir -p ${POST_OUTPUT}/circular_output
          cd ${POST_OUTPUT}/circular_output

          # Completeness chcek  Output file: permute_scaffold.fasta
          circ-permute.pl -in $Assembly -split 600

          tgsgapcloser --tgstype pb --scaff ${POST_OUTPUT}/circular_output/permute_scaffold.fasta --reads ${RAW_FASTQ}/${SAMPLE_NAME}.fasta --output circular --ne >circular.log 2>circular.err
          echo "circular.gap_fill_detail" >> ${LOG_DIR}/post-processing.log.txt
          cat circular.gap_fill_detail >> ${LOG_DIR}/post-processing.log.txt
          echo "---------------------------------" >> ${LOG_DIR}/post-processing.log.txt

          cd ${POST_OUTPUT}
          cp ${POST_OUTPUT}/circular_output/circular.scaff_seqs ${POST_OUTPUT}/final_circular_scaffolds.fasta

          quast.py -t 8 --glimmer --rna-finding --contig-thresholds 2500,5000,10000,50000,100000,500000,1000000 \
          -l $genome,${genome}.GapsClosed \
          -o post_quast $Assembly ${POST_OUTPUT}/circular_output/circular.scaff_seqs

      fi

fi


# Generate a yaml file for Annotation
# https://github.com/ncbi/pgap/wiki/Input-Files

cd $LOCAL
if head -n 1 ${ASSEMBLY_OUTPUT}/assembly.fasta | grep -q "circular" ; then
    topology_name="circular"
else
    topology_name="linear"
fi

export ${topology_name}
strain_name=$SAMPLE_NAME
export ${strain_name}

( echo "cat <<EOF >submol.yaml";
  cat template.yaml;
  echo "EOF";
) >temp.yml

. temp.yml
mv submol.yaml ${ANNO_OUTPUT}

rm -f temp.yml


######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
# rm -rf "${OUTPUTDIR}"
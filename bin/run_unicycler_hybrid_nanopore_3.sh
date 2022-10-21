#!/bin/bash -x

set -e
set -u
set -o pipefail

START_TIME=$SECONDS
export PATH="/opt/conda/bin:${PATH}"

LOCAL=$(pwd)
coreNum=${coreNum:-8}
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


mv $fastq1 "${QC_FASTQ}/read1_sampled.fastq.gz"
mv $fastq2 "${QC_FASTQ}/read2_sampled.fastq.gz"
mv $longreads "${QC_FASTQ}/long_trimmed.fastq.gz"
mv $assembly "${ASSEMBLY_OUTPUT}/assembly.fasta"

## Build the database
# --threads ${coreNum} \
bowtie2-build --seed 42 \
    ${ASSEMBLY_OUTPUT}/assembly.fasta \
    ${LOCAL_DB_PATH}/${LOCAL_DB_NAME} |\
    tee -a ${LOG_DIR}/bowtie2_build_db_index.log

bowtie2 --sensitive-local -p ${coreNum} --seed 42 -x ${LOCAL_DB_PATH}/${LOCAL_DB_NAME} \
-1 "${QC_FASTQ}/read1_sampled.fastq.gz" -2 "${QC_FASTQ}/read2_sampled.fastq.gz" | \
samtools view -@ ${coreNum} -bh -o "${BWT_OUTPUT}/contigs_aligned.bam" - | \
tee -a ${LOG_DIR}/bowtie2_mapping.log.txt

#-S "${BWT_OUTPUT}/contigs_aligned.sam"
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
#get coverage cal_depth
RunDate=`cat ${ASSEMBLY_OUTPUT}/unicycler.log | grep -P "^Starting Unicycler" | cut -d "(" -f2 | cut -d " " -f1`
shortDepth=`cat ${BAMQC_OUTPUT}/genome_results.txt | grep -P "mean coverageData =" | cut -d= -f2 | cut -f1 -dX`
longDepth=`cat ${BAMQC_OUTPUT2}/genome_results.txt | grep -P "mean coverageData =" | cut -d= -f2 | cut -f1 -dX`

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


echo 'SampleName,Date,LongDepth,ShortDepth,TotalLongReads,TotalShortReads,TotalContigs,longAfterSamplingPct,shortAfterSamplingPct,longMappingPct,shortMappingPct,AvgLongReadLength,LongByShortCoverage,AssemblyStatus' > ${LOG_DIR}/read_stats.csv
echo ${SAMPLE_NAME}','${RunDate}','${longDepth}','${shortDepth}','${totalLongReads}','${totalShortReads}','${totalContigs}','$longUsedRate','$shortUsedRate','${long_mapping}','${short_mapping}','${AvgLongReadLength}','${LongByShortCoverage}','${status} >> ${LOG_DIR}/read_stats.csv

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

      Assembly=${ASSEMBLY_OUTPUT}/assembly.fasta  #"/data/BrianYu/nanopore_hybird/CloseGaps/Blautia-sp-KLE/Blautia-sp-KLE.assembly.fasta"
      Nanopore_reads="${RAW_FASTQ}/long.fastq.gz" #"/data/BrianYu/nanopore_hybird/CloseGaps/Blautia-sp-KLE/Blautia-sp-KLE-1732-HM-1032__barcode15.fastq.gz"

      genome=${SAMPLE_NAME} #"Blautia-sp-KLE"

      # Align Long reads to drraft assembly

      minimap2 -t 8 $Assembly $Nanopore_reads > ${RAW_FASTQ}/aln.mm

      # Form scaffolds
      java -jar /home/gitdir/lrscaf/target/LRScaf-1.1.9.jar -c $Assembly -a ${RAW_FASTQ}/aln.mm -t mm -o ${POST_OUTPUT}/lrscaf_output

      #Convert fastq to fasta
      reformat.sh qin=33 overwrite=t ignorebadquality=t in=$Nanopore_reads out=${RAW_FASTQ}/${SAMPLE_NAME}.fasta



      if [ -f "${POST_OUTPUT}/lrscaf_output/scaffolds_summary.info" ]; then
          echo "scaffolds_summary.info" >>  ${LOG_DIR}/post-processing.log.txt
          cat ${POST_OUTPUT}/lrscaf_output/scaffolds_summary.info >> ${LOG_DIR}/post-processing.log.txt
          echo "---------------------------------" >> ${LOG_DIR}/post-processing.log.txt

          mkdir -p ${POST_OUTPUT}/tgsgapcloser_output
          cd ${POST_OUTPUT}/tgsgapcloser_output
          #
          tgsgapcloser --scaff ${POST_OUTPUT}/lrscaf_output/scaffolds.fasta --reads ${RAW_FASTQ}/${SAMPLE_NAME}.fasta --output gaps --ne >gaps.log 2>gaps.err

          echo "gaps.gap_fill_detail" >> ${LOG_DIR}/post-processing.log.txt
          cat gaps.gap_fill_detail >> ${LOG_DIR}/post-processing.log.txt
          echo " " >> ${LOG_DIR}/post-processing.log.txt

          cd ${POST_OUTPUT}

          cp  ${POST_OUTPUT}/tgsgapcloser_output/gaps.scaff_seqs ${POST_OUTPUT}/gapfilled_scaffolds.fasta
          mkdir -p ${POST_OUTPUT}/circular_output
          cd ${POST_OUTPUT}/circular_output

          # Completeness chcek  Output file: permute_scaffold.fasta
          circ-permute.pl -in ${POST_OUTPUT}/tgsgapcloser_output/gaps.scaff_seqs -split 600

          tgsgapcloser --scaff ${POST_OUTPUT}/circular_output/permute_scaffold.fasta --reads ${RAW_FASTQ}/${SAMPLE_NAME}.fasta --output circular --ne >circular.log 2>circular.err
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

          tgsgapcloser --scaff ${POST_OUTPUT}/circular_output/permute_scaffold.fasta --reads ${RAW_FASTQ}/${SAMPLE_NAME}.fasta --output circular --ne >circular.log 2>circular.err
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
strain_name=`echo "${SAMPLE_NAME}" | tr "-" "_"`
export ${strain_name}

( echo "cat <<EOF >submol.yaml";
  cat /usr/local/bin/template.yaml;
  echo "EOF";
) >temp.yml

. temp.yml
mv submol.yaml "${ANNO_OUTPUT}/submol.yaml"

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

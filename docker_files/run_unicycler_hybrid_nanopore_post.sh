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
coverage=${coverage:-80}

# if a genomeSize is give, then overwrite the TargetBase parameters
ShortTargetBase=$((genomeSize*coverage*1))
LongTargetBase=$((genomeSize*coverage))

# Make sure neither the short read set nor long read set exceeds 500 Mbp
if [ $ShortTargetBase -gt 500000000 ]
then
   ShortTargetBase=500000000
fi

if [ $LongTargetBase -gt 500000000 ]
then
   LongTargetBase=500000000
fi
# s3 inputs from env variables
#fastq1="${1}"
#fastq2="${2}"
#longreads="${3}" input long reads in fastq format
#S3OUTPUTPATH="${4}"
#genomeSize="${5}" genome size in bp
#UAssembly=${6}

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

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}" "${QC_FASTQ}" "${POST_OUTPUT}"
mkdir -p "${ASSEMBLY_OUTPUT}" "${QUAST_OUTPUT}" "${FASTQC_OUTPUT}" "${FASTQC_OUTPUT2}"
mkdir -p "${LOCAL_DB_PATH}" "${BWT_OUTPUT}"  "${BAMQC_OUTPUT}" "${BAMQC_OUTPUT2}" "${RAW_NANOPORE}"
trap '{ rm -rf ${OUTPUTDIR} ; exit 255; }' 1

hash_kmer=${hash_kmer:-51}

# Copy fastq.gz files from S3, only 2 files per sample
aws s3 cp --quiet ${longreads} "${RAW_FASTQ}/long.fastq.gz"
aws s3 cp --quiet ${UAssembly} "${RAW_FASTQ}/uassembly.fasta"

      Assembly="${RAW_FASTQ}/uassembly.fasta"  #"/data/BrianYu/nanopore_hybird/CloseGaps/Blautia-sp-KLE/Blautia-sp-KLE.assembly.fasta"
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
          cat lrscaf_output/scaffolds_summary.info >> ${LOG_DIR}/post-processing.log.txt
          echo "---------------------------------" >> ${LOG_DIR}/post-processing.log.txt

          mkdir -p ${POST_OUTPUT}/tgsgapcloser_output
          cd ${POST_OUTPUT}/tgsgapcloser_output
          #
          tgsgapcloser --scaff ${POST_OUTPUT}/lrscaf_output/scaffolds.fasta --reads ${RAW_FASTQ}/${SAMPLE_NAME}.fasta --output gaps --ne >gaps.log 2>gaps.err

          echo "gaps.gap_fill_detail" >> ${LOG_DIR}/post-processing.log.txt
          cat gaps.gap_fill_detail >> ${LOG_DIR}/post-processing.log.txt
          echo " " >> ${LOG_DIR}/post-processing.log.txt

          cd ${POST_OUTPUT}

          #ln -s -f ./tgsgapcloser_output/gaps.scaff_seqs gapfilled_scaffolds.fasta

          mkdir -p ${POST_OUTPUT}/circular_output
          cd ${POST_OUTPUT}/circular_output

          # Completeness chcek  Output file: permute_scaffold.fasta
          perl /mnt/circ-permute.pl -in ${POST_OUTPUT}/tgsgapcloser_output/gaps.scaff_seqs

          tgsgapcloser --scaff ${POST_OUTPUT}/circular_output/permute_scaffold.fasta --reads ${RAW_FASTQ}/${SAMPLE_NAME}.fasta --output circular --ne >circular.log 2>circular.err
          echo "circular.gap_fill_detail" >> ${LOG_DIR}/post-processing.log.txt
          cat circular.gap_fill_detail >> ${LOG_DIR}/post-processing.log.txt

          cd ${POST_OUTPUT}
          #ln -s ./circular_output/circular.scaff_seqs circular_scaffolds.fasta

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
          perl /mnt/circ-permute.pl -in $Assembly

          tgsgapcloser --scaff ${POST_OUTPUT}/circular_output/permute_scaffold.fasta --reads ${RAW_FASTQ}/${SAMPLE_NAME}.fasta --output circular --ne >circular.log 2>circular.err
          echo "circular.gap_fill_detail" >> ${LOG_DIR}/post-processing.log.txt
          cat circular.gap_fill_detail >> ${LOG_DIR}/post-processing.log.txt
          echo "---------------------------------" >> ${LOG_DIR}/post-processing.log.txt

          cd ${POST_OUTPUT}
          #ln -s -f ./circular_output/circular.scaff_seqs circular_scaffolds.fasta

          quast.py -t 8 --glimmer --rna-finding --contig-thresholds 2500,5000,10000,50000,100000,500000,1000000 \
          -l $genome,${genome}.GapsClosed \
          -o post_quast $Assembly ${POST_OUTPUT}/circular_output/circular.scaff_seqs

      fi






######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
# rm -rf "${OUTPUTDIR}"

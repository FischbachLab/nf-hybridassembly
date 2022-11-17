Sample scripts
====================

# A simple script showing the hybrid assembly example for the Nextflow framework.


```{bash}
aws batch submit-job \
  --job-name nf-hybrid-assembly \
  --job-queue priority-maf-pipelines \
  --job-definition nextflow-production \
  --container-overrides command="s3://nextflow-pipelines/nf-hybrid-assembly, \
"--reads1","s3://czb-seqbot/fastqs/211110_A00111_0817_AH2CNWDSX3/CZBMI-NICHE_PriorityEffect/LibPlateF4_CJ481_Prevotella-buccae_R1.fastq.gz", \
"--reads2","s3://czb-seqbot/fastqs/211110_A00111_0817_AH2CNWDSX3/CZBMI-NICHE_PriorityEffect/LibPlateF4_CJ481_Prevotella-buccae_R2.fastq.gz", \
"--long_reads","s3://czb-seqbot/nanopore/211105_MN19452_0078_FAQ78612/L-la_C-bi_B-ps_B-fr_A-e_A-s_E-l_St_F-m_C-a_P-b_L-l/20211106_0028_MN19452_FAQ78612_e188e18f/01_BASECALLED/Prevotella-buccae_CJ481__barcode11.fastq.gz", \
"--coverage", "100", \
"--output_path", "s3://genomics-workflow-core/Results/HybridAssembly/biohub" "
```

## The latest version with a seedfile as the input
```{bash}
aws batch submit-job \
  --job-name nf-hybrid-assembly \
  --job-queue priority-maf-pipelines \
  --job-definition nextflow-production \
  --container-overrides command="s3://nextflow-pipelines/nf-hybrid-assembly, \
"--seedfile", "s3://genomics-workflow-core/Results/HybridAssembly/BrianYu_20220311/HotspringMatIsolate_Chloroflexus-Stock92-20220203.tsv", \
"--coverage", "100", \
"--output_path", "s3://genomics-workflow-core/Results/HybridAssembly/BrianYu_20220311" "
```

### long reads only assembly *** NOT WORKING using long.nf rather than default name main.nf ***
#### sampleRate must be greater than 0 and less than 100
```{bash}
aws batch submit-job \
  --job-name nf-hybrid-assembly \
  --job-queue priority-maf-pipelines \
  --job-definition nextflow-production \
  --container-overrides command="s3://nextflow-pipelines/nf-hybrid-assembly, \
"--sampleRate", "50", \
"--seedfile", "s3://genomics-workflow-core/Results/HybridAssembly/Nathan_20220301/1.seedfile.tsv", \
"--output_path", "s3://genomics-workflow-core/Results/HybridAssembly/Nathan_20220301/tmp" "
```

# Seedfile example
## Note that the seedfile is a tab-separated values file without header
## The format is sample_name, prefix, short_R1, short_R2 and long_reads

```{bash}
SH0001651-00109 20221021 s3://maf-sequencing/Illumina/MITI-MCB/G04_MITI001_SH0001651_R1.fastq.gz s3://maf-sequencing/Illumina/MITI-MCB/G04_MITI001_SH0001651_R2.fastq.gz s3://maf-sequencing/nanopore/MITI/Combined/SH0001651-00109.combined.fastq.gz
```
# MITI MCB Hybrid assembly example
### Final output path: s3://genomics-workflow-core/Results/HybridAssembly/MITI-MCB/sample_name/prefix/
```{bash}
aws batch submit-job \
  --job-name nf-hybrid-assembly \
  --job-queue priority-maf-pipelines \
  --job-definition nextflow-production \
  --container-overrides command="FischbachLab/nf-hybridassembly, \
"--seedfile", "s3://genomics-workflow-core/Results/HybridAssembly/20221018.11.seedfile.tsv", \
"--output_path", "s3://genomics-workflow-core/Results/HybridAssembly/MITI-MCB" "
```

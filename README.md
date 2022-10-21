Hello world script
====================

A simple script showing the hybrid assembly example for the Nextflow framework.


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

# updated version with a seedfile as the input
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

# long reads only assembly *** NOT WORKING using long.nf rather than default name main.nf ***
# sampleRate must be greater than 0 and less than 100
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

# MITI MCB Hybrid assembly example
## Note that the seedfile has no header and the format is sample_name, short_R1, short_R2 and long_reads
```{bash}
aws batch submit-job \
  --job-name nf-hybrid-assembly \
  --job-queue priority-maf-pipelines \
  --job-definition nextflow-production \
  --container-overrides command="s3://nextflow-pipelines/nf-hybrid-assembly, \
"--seedfile", "s3://genomics-workflow-core/Results/HybridAssembly/20221018.11.seedfile.tsv", \
"--output_path", "s3://genomics-workflow-core/Results/HybridAssembly/MITI-MCB" "
```

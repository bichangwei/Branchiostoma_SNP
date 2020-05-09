# Branchiostoma_SNP
This repository contains some scripts for variant calling and inferring the demographic history of Branchiostoma.

1. reads_mapping_and_markdup.sh
The purpose of this script is to generate input bam files for variant calling and PSMC inference. Before mapping reads to the reference genome, you should define the software and file paths according to your own directory. Additionally, you should add some indexing for your reference genome if not present.

2. variant_calling.sh
This script is applied to call variants from the generated bam files, which is generated on the basis of GATK Best Practice Workflow. Since there was no known dbSNP dataset of Branchiostoma for BQSR, we generated a SNP dataset as the known dbSNP dataset using this script. If your studied species is not the model organism and there is no known dbSNP dataset, you can use our approach to create your own dbSNP dataset. The parameters used in "VariantFiltration" are recommended by GATK forum.

3. variant_filter_and_annotate.sh
This script is applied to filter low quality variants and annotate the high quality variants for further analysis. Apart from the recommended parameters, we set the lowDP and highDP as 0.1- and 2-fold of the average depth of the sequencing to filter overcalling or miscalling variants. The variant annotation software SNPeff v.4.3 and ANNOVAR v.2016-02-01 were used in this study to annotate the detected variants. Most importantly, you should build the local annotation database before running the annotation software.

4. run.PSMC.sh
This script is applied to infer the historical population size dynamics of Branchiostoma population. First, you should call the consequence for each sample using samtools mpileup and bcftools. Then, you should convert the fastq format files to the input format for PSMC using the "fq2psmcfa" program provided with the PSMC software. Next, appling the main "psmc" program to infer the population size history. The "psmc_plot.pl" script was finally used to visualize the PSMC results. 
The following figure shows the population size dynamics of Branchiostoma belcheri during the Pleistocene.

![image](https://github.com/bichangwei/Branchiostoma_SNP/blob/master/images/B_belcheri_PSMC.jpg)


One may want to perform 100 bootstrap replicates to assess the uncertainty in effective population size estimates. Then, you should randomly sample with replacement 5-Mb sequence segments using the '-b' option in "psmc" program. 
The following figure shows the population size dynamics of Branchiostoma floridae after 100 bootstrap replicates.

![image](https://github.com/bichangwei/Branchiostoma_SNP/blob/master/images/Bootstrap_B_floridae.jpg)


5. Estimate the best sequencing depth for PSMC.
In this workflow, we also provide a method to estimate the best sequencing depth for PSMC if you have high-depth sequencing data.
First, you should use "samtools view" or "sambamba view" to downsample your bam format files into 5x, 10x, 15x, 20x, 25x, 30x, 35x,...
Then, repeat the step provided in the "run.PSMC.sh" script.
The following figure shows the PSMC results of different sequencing depths.

![image](https://github.com/bichangwei/Branchiostoma_SNP/blob/master/images/Diff_depth_PSMC.jpg)


Reference

https://github.com/lh3/psmc/blob/master/README

https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-

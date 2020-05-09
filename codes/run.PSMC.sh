#! /bin/bash
GENOME=/home/bicw/Desktop/work/Branchiostoma/Reference/Branchiostoma_18h27.fa
SAMTOOLS=/home/bicw/bin/samtools_v1.3.1/samtools
BCFTOOLS=/home/bicw/bin/bcftools_v1.3.1/bcftools
SAMBAMBA=/home/bicw/bin/sambamba_v0.6.6/sambamba
BAMDir=/home/bicw/Desktop/work/Branchiostoma/BQSR_recalBAM
OutputDir=/home/bicw/Desktop/work/Branchiostoma/PSMC
#It is recommended to set minDP to a third of the average depth and maxDP to twice.
minDP=$1
maxDP=$2
for sample in BB_ZJ{1..5} BB_XM{1..5}
do
	#Calling consensus sequence for each sample
	${SAMTOOLS} mpileup -C 50 -q 25 -uf ${GENOME} ${BAMDir}/${sample}.BQSR.bam | ${BCFTOOLS} call -c - | /usr/share/samtools/vcfutils.pl vcf2fq -d ${minDP} -D ${maxDP} > ${OutputDir}/${sample}.BQSR.psmc.fq
	#-C: adjust mapping quality; 
	#-q: determine the cutoffs for mapQ
	#-u: generate uncompressed VCF output
	#-f: the reference fasta used (needs to be indexed)

	#Convert fastq files to the input format for PSMC
	/home/bicw/bin/psmc/utils/fq2psmcfa -q20 ${OutputDir}/${sample}.BQSR.psmc.fq > ${OutputDir}/${sample}.BQSR.psmcfa
	/home/bicw/bin/psmc/psmc -N25 -t13 -r5 -p "2+58*1+2*2" -o ${OutputDir}/${sample}.BQSR.psmc ${OutputDir}/${sample}.BQSR.psmcfa
	#-N: maximum number of iterations
	#-t: maximum 2N0 coalescent time
	#-r: initial theta/rho ratio
	#-p: pattern of parameters
	#the `-p' and `-t' options are manually chosen such that after 20 rounds of iterations
done

#visualize the population size history
/home/bicw/bin/psmc/utils/psmc_plot.pl -g 1.5 -x 1000 -X 5000000 -Y 180 -u 3e-09 -M BB_ZJ1,BB_ZJ2,BB_ZJ3,BB_ZJ4,BB_ZJ5,BB_XM1,BB_XM2,BB_XM3,BB_XM4,BB_XM5 BB_Total_3e-09 ${OutputDir}/BB_ZJ1.BQSR.psmc ${OutputDir}/BB_ZJ2.BQSR.psmc ${OutputDir}/BB_ZJ3.BQSR.psmc ${OutputDir}/BB_ZJ4.BQSR.psmc ${OutputDir}/BB_ZJ5.BQSR.psmc ${OutputDir}/BB_XM1.BQSR.psmc ${OutputDir}/BB_XM2.BQSR.psmc ${OutputDir}/BB_XM3.BQSR.psmc ${OutputDir}/BB_XM4.BQSR.psmc ${OutputDir}/BB_XM5.BQSR.psmc

#-g: bunber of years per generation
#-x: minimum generations
#-X: maximum generations
#-Y: maximum popsize
#-u: per-generation mutationb rate
#-M: multiline mode

#100 bootstrap replicates
for sample in BB_ZJ{1..5} BB_XM{1..5}
do
	mkdir ${OutputDir}/${sample}
	#randomly sampling with replacement 5-Mb sequence segments
	/home/bicw/bin/psmc/utils/splitfa ${OutputDir}/${sample}.BQSR.psmcfa > ${OutputDir}/${sample}/${sample}.split.psmc.fa
	
	#Convert fastq files to the input format for 100 bootstrap replicates of PSMC
	seq 100 | xargs -i echo /home/bicw/bin/psmc/psmc -N25 -t13 -r5 -b -p "2+58*1+2*2" -o ${OutputDir}/${sample}/${sample}.round-{}.psmc ${OutputDir}/${sample}/${sample}.split.psmc.fa | sh
	
	#merge all PSMC results of 100 bootstrao replicates for each sample
	cat ${OutputDir}/${sample}.BQSR.psmc ${OutputDir}/${sample}/${sample}.round-*.psmc > ${OutputDir}/${sample}.100bootstrap.psmc
done

#visualize PSMC results of 100 bootstrap replicates for each sample
/home/bicw/bin/psmc/utils/psmc_plot.pl -g 1.5 -x 1000 -X 5000000 -Y 180 -u 3e-09 ${sample} ${OutputDir}/${sample}.100bootstrap.psmc
#merge all PSMC results of 100 bootstrap replicates and draw better pictures in Adobe Illustrator 2020

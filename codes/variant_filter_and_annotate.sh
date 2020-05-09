#! /bin/bash
GATK4=/home/bicw/bin/GATK/gatk-4.0.2.1/gatk
GENOME=/home/bicw/Desktop/work/Branchiostoma/Reference/Branchiostoma_18h27.fa

for sample in BB_ZJ1
do
########1. Variant filteration
	########1.1. select SNPs from the VCF files
	${GATK4} --java-options -Xmx20g SelectVariants -select-type SNP -V ${sample}/${sample}.BQSR.vcf -O ${sample}/${sample}.BQSR.SNP.vcf

	########1.2. filter bad SNPs with strict parameters; the lowDP and highDP were set as 0.1- and 2-fold of the average depth of the sequencing
	${GATK4} --java-options -Xmx20g VariantFiltration -V ${sample}/${sample}.BQSR.SNP.vcf \
	--filter-name "lowMQ" --filter-expression "MQ < 40.0" \
	--filter-name "lowMQRankSum" --filter-expression "MQRankSum < -12.5" \
	--filter-name "highFS" --filter-expression "FS > 60.0" \
	--filter-name "lowReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0" \
	--filter-name "lowQD" --filter-expression "QD < 2.0"  \
	--filter-name "highSOR" --filter-expression "SOR > 3.0" \
	--filter-name "lowDP" --filter-expression "DP < 13.0" \
	--filter-name "highDP" --filter-expression "DP > 245.0" \
	-O ${sample}/${sample}.BQSR.SNP.filter.vcf
	
	########1.3. select SNPs from the VCF files
        ${GATK4} --java-options -Xmx20g SelectVariants -select-type INDEL -V ${sample}/${sample}.BQSR.vcf -O ${sample}/${sample}.BQSR.Indel.vcf

	########1.4. filter bad Indels with strict parameters; the lowDP and highDP were set as 0.1- and 2-fold of the average depth of the sequencing
	${GATK4} --java-options -Xmx20g VariantFiltration -V ${sample}/${sample}.BQSR.Indel.vcf \
	--filter-name "highFS" --filter-expression "FS > 200.0" \
	--filter-name "lowReadPosRankSum" --filter-expression "ReadPosRankSum < -20.0" \
	--filter-name "lowQD" --filter-expression "QD < 2.0" \
	--filter-name "highSOR" --filter-expression "SOR > 10.0" \
	--filter-name "lowDP" --filter-expression "DP < 13.0" \
	--filter-name "highDP" --filter-expression "DP > 245.0" \
	-O ${sample}/${sample}.BQSR.Indel.filter.vcf

	########1.5. exclude the bad variants from the VCF files and merge the SNPs and Indels into a single file
	${GATK4} --java-options -Xmx20g SelectVariants --exclude-filtered  -V ${sample}/${sample}.BQSR.SNP.filter.vcf -O ${sample}/${sample}.BQSR.SNP.filtered.vcf

	${GATK4} --java-options -Xmx20g SelectVariants --exclude-filtered  -V ${sample}/${sample}.BQSR.Indel.filter.vcf -O ${sample}/${sample}.BQSR.Indel.filtered.vcf

	${GATK4} --java-options -Xmx20g MergeVcfs -I ${sample}/${sample}.BQSR.SNP.filtered.vcf -I ${sample}/${sample}.BQSR.Indel.filtered.vcf -O ${sample}/${sample}.merge.BQSR.filtered.vcf

########2. Variant annotation with SNPeff v.4.3
	########2.1. build the annotation database for Branchiostoma_18h27 reference genome
	java -jar /home/bicw/bin/snpEff/snpEff.jar build -gff3 -v Branchiostoma_18h27
	
	########2.2. annotate the VCF files
	java -Xmx30g -jar /home/bicw/bin/snpEff/snpEff.jar Branchiostoma_18h27 ${sample}/${sample}.merge.BQSR.filtered.vcf > ${sample}/${sample}.merge.BQSR.filtered.eff

########3. variant annotation with ANNOVAR v.2016-02-01
	########3.1. build the annotation database for Branchiostoma_18h27 reference genome
	/usr/bin/gffread -T Branchiostoma_18h27.gff -o Branchiostoma_18h27.gtf
	/usr/bin/gtfToGenePred -genePredExt Branchiostoma_18h27.gtf Branchiostoma_18h27_refGene.txt
	perl /home/bicw/bin/Annovar/retrieve_seq_from_fasta.pl --format refGene --seqfile Branchiostoma_18h27.fa Branchiostoma_18h27_refGene.txt -outfile Branchiostoma_18h27_refGeneMrna.fa
	
	########3.2. annotate the VCF files with Annovar
	perl /home/bicw/bin/Annovar/convert2annovar.pl -format vcf4old ${sample}/${sample}.merge.BQSR.filtered.vcf > ${sample}/${sample}.merge.BQSR.filtered.annovar.txt
	java -jar /home/bicw/bin/Annovar/annotate_variation.pl -buildver ${sample}/${sample}.merge.BQSR.filtered.annovar.txt /home/bicw/bin/Annovar/Branchiostoma_18h27/ --geneanno --outfil ${sample}/${sample}.merge.BQSR.filtered.annovar	
done

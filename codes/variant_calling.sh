#! /bin/bash
GATK3=/home/bicw/bin/gatk/GenomeAnalysisTK.jar
GATK4=/home/bicw/bin/GATK/gatk-4.0.2.1/gatk
GENOME=/home/bicw/Desktop/work/Branchiostoma/Reference/Branchiostoma_18h27.fa
SAMTOOLS=/home/bicw/bin/samtools_v1.3.1/samtools
BCFTOOLS=/home/bicw/bin/bcftools_v1.3.1/bcftools

######## Creating Dictionary for reference genome before variant calling
$GATK4 CreateSequenceDictionary -R $GENOME -O /home/bicw/Desktop/work/Branchiostoma/Reference/Branchiostoma_18h27.dict

for sample in BB_ZJ{1..5} BB_XM{1..5}
do
#########1. Variant Calling Round1:
	########1.1. Identify what regions need to be realigned
	time java -Xmx20g -jar $GATK3 -R $GENOME -T RealignerTargetCreator -I $sample/$sample.mem.M.filtered.markdup.bam -o $sample/$sample.realn.intervals  && echo "===== GATK3 RealignerTargetCreator done ====="

	########1.2. Perform the actual realignment
	time java -Xmx20g -jar $GATK3 -R $GENOME -T IndelRealigner -targetIntervals $sample/$sample.realn.intervals -I $sample/$sample.mem.M.filtered.markdup.bam -o $sample/$sample.realn.bam && echo "===== GATK3 IndelRealigner done ====="

	########1.3. Variant call using samtools and bcftools
	time samtools  mpileup -ugf $GENOME  $sample/$sample.realn.bam | bcftools call -vmO v -o $sample/$sample.samtools.round1.vcf && echo "===== samtools call SNP done ===="

	########1.4. Variant call using GATK4 HaplotypeCaller
	time $GATK4 --java-options -Xmx20g HaplotypeCaller -R $GENOME -I $sample/$sample.mem.M.filtered.markdup.bam -O $sample/$sample.gatk.hc.round1.vcf && echo "===== Round1 GATK call SNP done ===="

        ########1.5. Index the generated vcf files
        $GATK4 --java-options -Xmx20g IndexFeatureFile -F $sample/$sample.samtools.round1.vcf
	$GATK4 --java-options -Xmx20g IndexFeatureFile -F $sample/$sample.gatk.hc.round1.vcf

#########2. Recalibrate Base Quality Scores(BQSR)
	########2.1. select the same variation shared by GATK and SAMtools
	$GATK4 --java-options -Xmx20g SelectVariants --variant $sample/$sample.gatk.hc.round1.vcf --concordance $sample/$sample.samtools.round1.vcf -O $sample/$sample.samtools.hc.round1.vcf
	
	########2.2. filter bad variants with strict parameters and select the filtered variants as known dbSNP for BQSR
	$GATK4 --java-options -Xmx20g VariantFiltration \
        -V $sample/$sample.samtools.hc.round1.vcf \
        --filter-name "lowQD" --filter-expression "QD < 10.0" \
        --filter-name "lowMQ" --filter-expression "MQ < 50.0" \
        --filter-name "highFS" --filter-expression "FS > 10.0" \
        --filter-name "lowMQRankSum" --filter-expression "MQRankSum < -5.0" \
        --filter-name "lowReadPosRankSum" --filter-expression "ReadPosRankSum < -8.0" \
        -O $sample/$sample.samtools.hc.round1.filter.vcf

	$GATK4 --java-options -Xmx20g SelectVariants --exclude-filtered  -V $sample/$sample.samtools.hc.round1.filter.vcf -O $sample/$sample.samtools.hc.round1.filtered.vcf
	
	########2.3. recalibrate base quality scores based on the above known dbSNP
	time $GATK4 --java-options -Xmx20g BaseRecalibrator -R $GENOME -I $sample/$sample.mem.M.filtered.markdup.bam --known-sites $sample/$sample.samtools.hc.round1.filtered.vcf -O $sample/$sample.racal.table && echo "===== BaseRecalibrator done ===="
	
	time $GATK4 --java-options -Xmx20g ApplyBQSR --bqsr-recal-file $sample/$sample.racal.table -R $GENOME -I $sample/$sample.mem.M.filtered.markdup.bam -O $sample/$sample.BQSR.bam && echo "===== ApplyBQSR done ===="

########3. Variant Calling Round2:
	time $GATK4 --java-options -Xmx20g HaplotypeCaller -R $GENOME -I $sample/$sample.BQSR.bam -O $sample/$sample.BQSR.vcf && echo "===== Round2 GATK call SNP done ===="
done

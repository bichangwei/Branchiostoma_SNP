#! /bin/bash
GATK4=/home/bicw/bin/GATK/gatk-4.0.2.1/gatk
GENOME=/home/bicw/Desktop/work/Branchiostoma/Reference/Branchiostoma_18h27.fa
DataDir=/home/bicw/Desktop/work/Branchiostoma/Sequencing_Data/cleandata
Fai=/home/bicw/Desktop/work/Branchiostoma/Reference/Branchiostoma_18h27.fa.fai
BWA=/home/bicw/bin/bwa_v0.7.15/bwa
SAMBAMBA=/home/bicw/bin/sambamba_v0.6.6/sambamba
QUALIMAP=/home/bicw/bin/qualimap_v.2.2.1/qualimap
SAMTOOLS=/home/bicw/bin/samtools_v1.3.1/samtools
# bwa index ${GENOME}
# samtools faidx ${GENOME}
for sample in BB_ZJ{1..5} BB_XM{1..5}
do
	mkdir ${sample}
        Read1=${DataDir}/${sample}/${sample}_1.clean.fq.gz
        Read2=${DataDir}/${sample}/${sample}_2.clean.fq.gz

	########1. Map the clean reads to the reference genome
        time ${BWA} mem -t 20 -M -k 19 -R "@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPL:ILLUMINA" ${GENOME}  $Read1 $Read2 | samtools view -Sb -t ${Fai} -o ${sample}/${sample}.mem.M.bam && echo "===== bwa mapping done ====="

	########2. Sort BAM files and filter secondary and supplementary alignments
        time ${SAMBAMBA} sort -F "not (secondary_alignment or supplementary)" -t 20 -p -l 9 ${sample}/${sample}.mem.M.bam -o ${sample}/${sample}.mem.M.filtered.bam && echo "===== sambamba sort and filter done"

	########3. Count the number of alignments for each FLAG type
        time ${SAMBAMBA} flagstat ${sample}/${sample}.mem.M.filtered.bam > ${sample}/${sample}.mem.M.filtered.flagstat && echo "===== sambamba flagstat done  ====="

	########4. Report information for the evaluation of the quality of the provided alignment data (a BAM file)
        time ${SUALIMAP} bamqc --java-mem-size=24G -c -nt 20 -outdir ${sample}/${sample} -outformat PDF:HTML -outfile ${sample}/${sample}.mem.M.filtered.pdf -os -bam ${sample}/${sample}.mem.M.filtered.bam && echo "===== qualimap done ====="

	########5. Mark PCR duplicate reads
        time ${GATK4} MarkDuplicates -I ${sample}/${sample}.mem.M.filtered.bam -O ${sample}/${sample}.mem.M.filtered.markdup.bam -M ${sample}/${sample}.mem.M.filtered.markdup_metrics.txt && echo "===== GATK markdup done ====="

	########6. Index the duplication marking BAM files
        time ${SAMTOOLS} index ${sample}/${sample}.mem.M.filtered.markdup.bam && echo "===== samtools index done ====="
done



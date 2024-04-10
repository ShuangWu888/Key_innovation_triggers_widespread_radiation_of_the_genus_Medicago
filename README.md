# Key_innovation_triggers_widespread_radiation_of_the_genus_Medicago

Key innovation triggers widespread radiation of the genus Medicago

## Genotype calling

```
# Mapped clean reads onto the Medicago sativa subsp. caerulea (voucher PI464715) reference genome using BWA
bwa mem -t 40 -M -R '@RG\tID:sample\tPL:illumina\tPU:illumina\tLB:sample\tSM:sample' reference.fa sample_1_clean.fq.gz sample_2_clean.fq.gz | samtools sort -O bam -T tmp/sample -o 1.bwa/sample.sort.bam
# Duplicated reads were removed using the “MarkDuplicates” option in Picard v.2.25.0
java -Xmx10g -jar picard.2.25.0.jar MarkDuplicates INPUT=1.bwa/sample.sort.bam OUTPUT=2.rehead/sample.dedup.bam METRICS_FILE=2.rehead/sample.dup.txt REMOVE_DUPLICATES=true ; samtools index 2.rehead/sample.dedup.bam
```

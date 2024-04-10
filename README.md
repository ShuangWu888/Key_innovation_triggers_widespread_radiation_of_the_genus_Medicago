# Key_innovation_triggers_widespread_radiation_of_the_genus_Medicago

Key innovation triggers widespread radiation of the genus Medicago

## Genotype calling

```
# Mapped these clean reads onto the Medicago sativa subsp. caerulea (voucher PI464715) reference genome33 using BWA
bwa mem -t 40 -M -R '@RG\tID:sample\tPL:illumina\tPU:illumina\tLB:sample\tSM:sample' reference.fa sample_1_clean.fq.gz sample_2_clean.fq.gz | samtools sort -O bam -T tmp/sample -o 1.bwa/sample.sort.bam
```

# Key_innovation_triggers_widespread_radiation_of_the_genus_Medicago

Key innovation triggers widespread radiation of the genus Medicago

## Genotype calling

```
#Mapped these clean reads onto the Medicago sativa subsp. caerulea (voucher PI464715) reference genome33 using BWA
bwa mem -t 40 -M -R '@RG\tID:Ms170\tPL:illumina\tPU:illumina\tLB:Ms170\tSM:Ms170' Msa.fa Ms170_1_clean.fq.gz Ms170_2_clean.fq.gz | samtools sort -O bam -T tmp/Ms170 -o 1.bwa/Ms170.sort.bam
```

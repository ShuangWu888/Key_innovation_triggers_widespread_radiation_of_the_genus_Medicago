# Key_innovation_triggers_widespread_radiation_of_the_genus_Medicago

Key innovation triggers widespread radiation of the genus Medicago

## Genotype calling

```
#mapped these clean reads onto the Medicago sativa subsp. caerulea (voucher PI464715) reference genome33 using BWA
bwa mem -t 40 -M -R '@RG\tID:Ms170\tPL:illumina\tPU:illumina\tLB:Ms170\tSM:Ms170' Msa.fa /data/cold00/user109/usb_2t/42.Medicago.resequencing/cleandata/Ms170_1_clean.fq.gz /data/cold00/user109/usb_2t/42.Medicago.resequencing/cleandata/Ms170_2_clean.fq.gz | /data/00/user/user113/software/samtools-1.10/samtools sort -O bam -T tmp/Ms170 -o 1.bwa/Ms170.sort.bam
```

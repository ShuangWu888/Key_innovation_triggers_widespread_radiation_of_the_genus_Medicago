# Key_innovation_triggers_widespread_radiation_of_the_genus_Medicago

Key innovation triggers widespread radiation of the genus Medicago

## Genotype calling
Mapped clean reads onto the Medicago sativa subsp. caerulea (voucher PI464715) reference genome using BWA
```
bwa mem -t 40 -M -R '@RG\tID:sample\tPL:illumina\tPU:illumina\tLB:sample\tSM:sample' reference.fa sample_1_clean.fq.gz sample_2_clean.fq.gz | samtools sort -O bam -T tmp/sample -o 1.bwa/sample.sort.bam
```
Duplicated reads were removed using the “MarkDuplicates” option in Picard v.2.25.0
```
java -Xmx10g -jar picard.2.25.0.jar MarkDuplicates INPUT=1.bwa/sample.sort.bam OUTPUT=2.rehead/sample.dedup.bam METRICS_FILE=2.rehead/sample.dup.txt REMOVE_DUPLICATES=true ; samtools index 2.rehead/sample.dedup.bam
```
Identified single nucleotide polymorphisms (SNPs)
```
java -Xmx10g -jar gatk-package-4.1.4.1-local.jar HaplotypeCaller -R reference.fa -I sample.dedup.bam -O get/sample.g.vcf.gz --emit-ref-confidence GVCF --native-pair-hmm-threads 30
```
Filtered SNPs
```
vcftools --gzvcf rm_indel_5bp.sativa132.SNP.vcf.gz --minDP 6 --maxDP 50 --minGQ 20 --max-missing 0.8 --recode --recode-INFO-all --stdout | gzip -c > DP_6-50_miss_0.2_sativa132.SNP.vcf.gz
```

## Identification of single copy genes
Coverage statistical analysis for each nuclear gene
```
java -jar bamstats05.jar -B reference.gff.CDS.bed sample.dedup.bam | gzip -c > 02.bed_cov.pl.out/sample.bed_cov.gz
```
Used a custom script to identification of single copy genes
```
perl identify_single_copy_gene.132sativa.pl
```

## Plastomes assembly and annotation
Assembly

[NOVOPlasty3.8.3](https://github.com/ndierckx/NOVOPlasty)

```
perl NOVOPlasty3.8.3.pl -c config.txt
```

[GetOrganelle](https://github.com/Kinggerm/GetOrganelle)

```
get_organelle_from_reads.py -1 sample_1_clean.fq.gz -2 sample_2_clean.fq.gz -t 10 -o sample.plastome -F embplant_pt -R 10
```

Annotation

[GeSeq](https://chlorobox.mpimp-golm.mpg.de/geseq.html)

[Geneious v.10.2.6](https://www.geneious.com/)


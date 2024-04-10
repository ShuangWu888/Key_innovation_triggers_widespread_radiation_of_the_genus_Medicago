# Key_innovation_triggers_widespread_radiation_of_the_genus_Medicago

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

## Phylogenomic analysis
### nuclear genome
#### For the concatenated whole genome SNPs
```
iqtree2 -s iqtree_DP_6-50_miss_0.2_sativa132_SNP.varsites.phy -st DNA -m GTR+ASC -B 1000 --bnni --prefix iqtree_DP_6-50_miss_0.2_sativa132_BS_SNP -T 50
```
#### For the concatenated CDS alignment of 7,990 single copy nuclear genes
##### The concatenated method
[RaxML](https://github.com/stamatak/standard-RAxML)
```
raxmlHPC-PTHREADS -s merge.cds.codon123.fa -n merge.cds.codon123.fa -m GTRCAT -f a -x 12345 -N 100 -p 12345 -T 30
```
##### The coalescent method
Reconstructed phylogenetic trees of each single copy nuclear gene using RAxML
```
raxmlHPC-PTHREADS -s MsaG000053.cds.nonmiss.fa -n MsaG000053.cds.nonmiss.fa -m GTRGAMMAI -f a -x 12345 -N 100 -p 12345 -T 4
```
Estimated a coalescent tree using ASTRAL v.5.6.3
```
java -Xmx200G -jar astral.5.6.3.jar -i 3.arstral.CDS.tree.BS10.tre -o 4.arstral.CDS.tree.BS10.individual.tre 2>4.arstral.CDS.tree.BS10.individual.tre.log
```
### Plastome
```
raxmlHPC-PTHREADS -s z131-70_cp_cds-half_gap.fasta -n z131-70_cp_cds-half_gap.fasta -m GTRCAT -f a -x 12345 -N 100 -p 12345 -T 30
```

## Divergence time estimation
[r8s](https://sourceforge.net/projects/r8s/)

[RelTime method in MEGA X](https://academic.oup.com/mbe/article/35/9/2334/5042667?login=false)

## Ancestral area reconstruction
Reconstructed ancestral distributions
```
Rscript M0_sativa-origin_maxareas6.R
```
Used biogeographic stochastic mapping (BSM) to calculate the type and number of biogeographical events
```
Rscript sativa-origin_BAYAREALIKE+J_BSM_maxareas6.R
```

## Ancestral character reconstruction
[Mesquite](https://www.mesquiteproject.org/)

## Macroevolutionary rate estimation

### Diversification rates estimation

#### Fit a series of time-dependent likelihood diversification models for Medicago using [RPANDA](https://github.com/hmorlon/PANDA/tree/master)
```
Rscript RPANDA-time.R
```
#### Estimated diversification rates using [BAMM](https://github.com/macroevolution/bamm/tree/master)
```
bamm -c BAMM_divcontrol.txt
```

### Traits rates estimation

[BAMM](https://github.com/macroevolution/bamm/tree/master)

```
bamm -c BAMM_traitcontrol.txt
```

### Niche rate estimation

#### Ordinated all 35 environmental variables using phylogenetic principal component analysis (PCA) implemented in the “phyl.pca” function of [phytools](https://github.com/liamrevell/phytools) package
```
Rscript phyl.pca.R
```

#### Used the trait model of [BAMM](https://github.com/macroevolution/bamm/tree/master) on the first axis of the phylogenetic PCA of the niche data
```
bamm -c BAMM_traitcontrol.txt
```

## Key innovation test
[HiSSE analysis in RevBayes](https://revbayes.github.io/tutorials/sse/hisse.html)

```
singularity exec -B /data /data/00/user/user109/software/RevBayes/ rb mcmc_HiSSE.Rev.txt
```

## Reconstruction of ancestral niches for annual and perennial Medicago species
[RevBayes](https://revbayes.github.io/tutorials/)

```
singularity exec -B /data /data/00/user/user109/software/RevBayes/ rb niche.charactor.rev.txt
```

## ILS simulation, hybridization inference, and admixture analysis
### Simplified phylogenetic trees (nuclear genome and plastome)
Selected one individual from each subclade (M1-M18) as a representative and reconstructed phylogenetic trees with Vicia sativa (DRR053677) as an outgroup
#### For nuclear genome
##### Used the concatenation method of RAxML based on the CDS dataset of 7,990 single copy nuclear genes
```
raxmlHPC-PTHREADS -s merge.cds.codon123.fa -n merge.cds.codon123.fa -m GTRCAT -f a -x 12345 -N 100 -p 12345 -T 30
```
##### Used the coalescent method of ASTRAL based on the CDS dataset of 7,990 single copy nuclear genes
Reconstructed phylogenetic trees of each single copy nuclear gene using RAxML
```
raxmlHPC-PTHREADS -s MsaG000053.cds.nonmiss.fa -n MsaG000053.cds.nonmiss.fa -m GTRGAMMAI -f a -x 12345 -N 100 -p 12345 -T 4
```
Estimated a coalescent tree using ASTRAL v.5.6.3
```
java -Xmx200G -jar astral.5.6.3.jar -i 3.arstral.CDS.tree.BS10.tre -o 4.arstral.CDS.tree.BS10.individual.tre 2>4.arstral.CDS.tree.BS10.individual.tre.log
```
#### For plastome
```
raxmlHPC-PTHREADS -s z19-70_cp_cds-half_gap.fa -n z19-70_cp_cds-half_gap.fa -m GTRCAT -f a -x 12345 -N 100 -p 12345 -T 15
```

### ILS simulation
#### Simulated 10,000 gene trees using the multispecies coalescent model of [Phybase](https://github.com/lliu1871/phybase) v.1.5 using the ASTRAL tree
```
Rscript phybase.R
```
#### Calculated gene-tree quartet frequencies of the five incongruent internal nodes in observed and simulated datasets using Twisst
##### Observed dataset
```
python twisst.py -t 2.arstral.CDS.simple.reroot.order.rename.tree -w Observed.node1.weights -g A M5 -g B M4 -g C M9 -g D M6
python twisst.py -t 2.arstral.CDS.simple.reroot.order.rename.tree -w Observed.node2.weights -g A M8 -g B M9 -g C M7 -g D M10
python twisst.py -t 2.arstral.CDS.simple.reroot.order.rename.tree -w Observed.node3.weights -g A M11 -g B M12 -g C M10 -g D M9
python twisst.py -t 2.arstral.CDS.simple.reroot.order.rename.tree -w Observed.node4.weights -g A M13 -g B M12 -g C M15 -g D M18
python twisst.py -t 2.arstral.CDS.simple.reroot.order.rename.tree -w Observed.node5.weights -g A M17 -g B M16 -g C M15 -g D M18
```
##### Simulated dataset
```
python twisst.py -t simulated.reroot.order.tree -w simulated.node1.weights -g A M5 -g B M4 -g C M9 -g D M6
python twisst.py -t simulated.reroot.order.tree -w simulated.node2.weights -g A M8 -g B M9 -g C M7 -g D M10
python twisst.py -t simulated.reroot.order.tree -w simulated.node3.weights -g A M11 -g B M12 -g C M10 -g D M9
python twisst.py -t simulated.reroot.order.tree -w simulated.node4.weights -g A M13 -g B M12 -g C M15 -g D M18
python twisst.py -t simulated.reroot.order.tree -w simulated.node5.weights -g A M17 -g B M16 -g C M15 -g D M18
```

### Hybridization inference
[PhyloNetworks](https://crsl4.github.io/PhyloNetworks.jl/latest/)
```
sh phyloNetworks.sh
```

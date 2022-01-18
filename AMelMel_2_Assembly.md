# Assembly

<!-- MDTOC maxdepth:6 firsth1:1 numbering:0 flatten:0 bullets:1 updateOnSave:1 -->

- [Assembly](#assembly)   
   - [1. Assembly of reads into contigs](#1-assembly-of-reads-into-contigs)   
   - [2. From contigs to chromosomes with a genetic map](#2-from-contigs-to-chromosomes-with-a-genetic-map)   
      - [2.1. Download fastq from SRA](#21-download-fastq-from-sra)   
         - [2.1.1 datasheet (supplementary table 2)](#211-datasheet-supplementary-table-2)   
         - [2.1.2. Download SRA data, validate download and write fastq](#212-download-sra-data-validate-download-and-write-fastq)   
      - [2.2. Map reads to the contigs](#22-map-reads-to-the-contigs)   
         - [Script for mapping:](#script-for-mapping)   
      - [2.3. Individual genotyping](#23-individual-genotyping)   
      - [2.4. Joint genotyping](#24-joint-genotyping)   
      - [2.5. Quality filters on variants](#25-quality-filters-on-variants)   
      - [Filter out any SNP or indel with missing genotype](#filter-out-any-snp-or-indel-with-missing-genotype)   

<!-- /MDTOC -->

## 1. Assembly of reads into contigs

* Raw PacBio reads from the single individual were assembled with Canu 1.3 (Koren *et al.* 2017) using standard parameters
* a first polishing of the assembly was done with quiver (version SMRT_Link v4.0.0) (https://github.com/PacificBiosciences/GenomicConsensus) using standard parameters.
* Error correction was then performed with Illumina reads from the same individual sequenced with an Illumina NovaSeq6000 instrument, producing over 28 000 000 reads (estimated raw sequencing depth = 33.7 X).
* Contigs were then assigned to chromosomes by alignment to the Amel4.5 reference genome using LAST (Frith and Kawaguchi 2015).

## 2. From contigs to chromosomes with a genetic map

Data from Liu et al (2015) was used to order and orient contigs with a genetic map.

In Liu et al. (2015), drones from 3 colonies were used to estimate recombination rates in the Amel4.5.

We will use the same sequencing data, to detect SNPs in our contigs. Minimizing recombination rates between contigs ends will then be used to order and orientate the contigs.

### 2.1. Download fastq from SRA

#### 2.1.1 datasheet (supplementary table 2)

* All the sequences of Liu et al. (2015) are in SRA under project SRP043350
  * Downloaded the corresponding SraRunTable from NCBI

  * Edited the table  as new sheet “Edited” in /Users/avignal/Documents/Stats/2017_LiuNewGeneticMap/SraRunTable.xlsx

    * Remove unnecessary columns
    * Add caste and colony information
    * Column “Duplicate”: "yes" for the 3 drones sequenced twice independantly

#### 2.1.2. Download SRA data, validate download and write fastq
* List: list of all 49 sequences to download
```
SRR1424586
SRR1425460
...
```

* Generate qarray file:

```
for i in `cat ../List`
do echo prefetch $i \&\& vdb-validate $i \&\& fastq-dump --split-files --gzip $i
done > Lance
```

* download sequences 5 at a time:

```
qarray -tc 5 Lance
```


### 2.2. Map reads to the contigs

#### Script for mapping:
* /genphyse/cytogen/seqapipop/ScriptsAV/mapLiu.sh
* Called by
  * /work/project/cytogen/Alain/LiuRecombination/MapColo1.bash
  * /work/project/cytogen/Alain/LiuRecombination/MapColo2.bash
  * /work/project/cytogen/Alain/LiuRecombination/MapColo3.bash
  * /work/project/cytogen/Alain/LiuRecombination/MapQueens.bash
* Colonies and queens done one after the others, due to memory constraints.

### 2.3. Individual genotyping
* 2017-11-06
* /genphyse/cytogen/seqapipop/ScriptsAV/snpsLiuOK.sh
    * Called by:
        * /work/project/cytogen/Alain/LiuRecombination/SnpsQueens.bash
        * /work/project/cytogen/Alain/LiuRecombination/SnpsColony1.bash
        * /work/project/cytogen/Alain/LiuRecombination/SnpsColony2.bash
        * /work/project/cytogen/Alain/LiuRecombination/SnpsColony3.bash

### 2.4. Joint genotyping
* List of samples to genotype:
```{bash, eval=FALSE}
ls ~/LiuRecombination/AlignmentsHaplotypeCaller/*/*/vcfs/*.g.vcf > GenotypeGVCF.list
```

* Then run the joint genotyping:
```{bash, eval=FALSE}
qsub -q workq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 geno-sitesLiuGvcf.bash
```

### 2.5. Quality filters on variants
* makes one vcf with SNPs and another with indels
```{bash, eval=FALSE}
qsub -q workq -l mem=4G -l h_vmem=48G FilterLiuGvcf.bash
```

### Filter out any SNP or indel with missing genotype
* Snps

```
qsub -b y -N filter 'vcftools \
--gzvcf ~/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF_snps_OK.vcf.gz \
--max-missing 1.0 \
--recode \
--out ~/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF_snps_OKNoMissing'
qsub -b y -N compress 'bgzip ~/LiuAllGenotypesGVCF_snps_OKNoMissing.recode.vcf'
qsub -b y -N index 'tabix -p vcf ~/LiuAllGenotypesGVCF_snps_OKNoMissing.recode.vcf.gz'
```

* Indels

```
qsub -b y -N filter 'vcftools \
--gzvcf ~/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF_indels_OK.vcf.gz \
--max-missing 1.0 \
--recode \
--out ~/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF_indels_OKNoMissing'
qsub -b y -N compress 'bgzip ~/LiuAllGenotypesGVCF_indels_OKNoMissing.recode.vcf'
qsub -b y -N index 'tabix -p vcf /LiuAllGenotypesGVCF_indels_OKNoMissing.recode.vcf.gz'
```

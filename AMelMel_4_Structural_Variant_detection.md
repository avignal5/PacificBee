# Structural Variant search





## minimap2



github link : https://github.com/lh3/minimap2



minimap2 is used to align the AMelMel.1 assembly of *Apis mellifera mellifera* to the HAv3.1 assembly, genetically similar to the  *Apis mellifera ligustica* or *Apis mellifera carnica* sub species.

* Output: bam file to use as input for svim-asm.

check for newer versions :

```bash
module avail -t 2>&1 | grep -i minimap2
```

On the genologin server, we tested the latest available minimap2 version, minimap2-2.19, but it produced a redundent error, so we used verion 2.5 instead :

```bash
module load bioinfo/minimap2-2.5
```



### Running minimap2 :

```bash
minimap2 -ax asm5 \
         ~/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
         ~/GCA_003314205.1_INRA_AMelMel_1.0_genomic.fna \
         > minimap_aln_genomes.sam
```

- -a	output in the SAM format (PAF by default)
- -x	preset (always applied before other options; see minimap2.1 for details)
     - map-pb/map-ont - PacBio CLR/Nanopore vs reference mapping
     - map-hifi - PacBio HiFi reads vs reference mapping
     - ava-pb/ava-ont - PacBio/Nanopore read overlap
     - asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5% sequence divergence
     - splice/splice:hq - long-read/Pacbio-CCS spliced alignment
     - sr - genomic short-read mapping

We use asm5 since both genomes are from the same species so have little divergence.



We then use samtools to change the file from sam to bam an prepare the bam file to be used by svim-asm.

Checking version and loading in genologin :

```bash
module avail -t 2>&1 | grep -i samtools
module load bioinfo/samtools-1.12
```

### Running samtools :

```bash
samtools view -S -b minimap_aln_genomes.sam > minimap_aln_genomes.bam
samtools sort minimap_aln_genomes.bam -o minimap_aln_genomes.sorted.bam
samtools index minimap_aln_genomes.sorted.bam
```

* The bam file was sorted and indexed to be used properly by svim-asm.



## svim-asm



github link : https://github.com/eldariont/svim-asm

As recomended, a conda environment was created and used for svim-asm (creat only once) :

```bash
conda create -n svimasm_env --channel bioconda svim-asm
```

Basic commands of conda environment :

```bash
conda env list # To see a list of all environments
conda activate svimasm_env # activate environment to use
conda deactivate # deactivate currently used invironment
```



### Running svim-asm :

```bash
svim-asm haploid . minimap_aln_genomes.sorted.bam /home/qboone/work/AMelMel/genome_abeille/GCF_003254395.2_Amel_HAv3.1_genomic.fna
```

The output will be "variants.vcf".



additional usful tools, vcftools (optional) :

```bash
conda install -c conda-forge scikit-allel
```

https://scikit-allel.readthedocs.io/en/stable/

http://alimanfoo.github.io/2017/06/14/read-vcf.html



For graphtyper2 the vcf file must be compressed and indexed. For that we use tabix : http://genometoolbox.blogspot.com/2014/09/how-to-index-vcf-file.html

http://www.htslib.org/doc/bgzip.html

```
bgzip -c variants.vcf > variants.vcf.gz
tabix -fp vcf variants.vcf.gz
```

Both commands are included in samtools that we loaded earlier.



## graphtyper2



github link : https://github.com/DecodeGenetics/graphtyper

installing graphtyper2 (linux) :

```bash
wget https://github.com/DecodeGenetics/graphtyper/releases/download/v2.7.1/graphtyper
chmod a+x graphtyper
```

Inputs for graphtyper:
* the genome's fasta file
* the compressed and indexed vcf file from svim-asm
* a file containing the paths to the bam files to analyse
* the coordinates of the target region

```bash
./graphtyper genotype_sv genome.fna variants.vcf.gz --sams=BAMLIST --region="chr:start-stop"
```

### Two regions tested

#### Scripts for genotyping:

* A 745 bp InDel on chromosome 2, present in AMelMel1.1 and absent in HAv3.1

```bash
#!/bin/bash

./graphtyper genotype_sv /home/gencel/vignal/save/Genomes/Abeille/HAv3_1_indexes/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
            /home/gencel/vignal/work/DepthChr3Inverstion/allInsDel/variants.vcf.gz \
            --sams=/home/gencel/vignal/work/DepthChr3Inverstion/allBAMsHAv3.list \
            --region="NC_037639.1:12047000-12048000"
```



* A 576 bp InDel on chromosome 10, absent in AMelMel1.1 and present in HAv3.1

```bash
#!/bin/bash

./graphtyper genotype_sv /home/gencel/vignal/save/Genomes/Abeille/HAv3_1_indexes/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
            /home/gencel/vignal/work/DepthChr3Inverstion/allInsDel/variants.vcf.gz \
            --sams=/home/gencel/vignal/work/DepthChr3Inverstion/allBAMsHAv3.list \
            --region="NC_037647.1:568701-773291"
```

#### To extract the genotypes for further analysis:

```bash
bcftools query -f'[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\n]' 012047000-012048000.vcf.gz > genotypes.txt
```

#### Unexpectedly, only the InDel on chromosome 2 could be genotyped. No results for the one on chromosome 10.
#### To see if this could be an effect of the genomes used as reference and query, the whole process was done again with the two genomes swapped for the minimap alignment:

### minipap

```bash
#!/bin/bash

module load bioinfo/minimap2-2.5

minimap2 -ax asm5 \
        ~/GCA_003314205.1_INRA_AMelMel_1.0_genomic.fna \
        ~/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
        > minimap_HAv3_on_AMel.sam
```

### graphtyper

* chromosome 2:


```bash
#!/bin/bash

../graphtyper genotype_sv /genphyse/cytogen/PacificBee/FromNcbi/GCA_003314205.1_INRA_AMelMel_1.0_genomic.fna \
            /home/gencel/vignal/work/DepthChr3Inverstion/allInsDel/HAv3_on_AMel/variants.vcf.gz \
            --sams=/home/gencel/vignal/work/DepthChr3Inverstion/allBAMsAMelMel.list \
            --region="CM010320.1:12211278-12214019"
```

* chromosome 10:

```bash
#!/bin/bash

../graphtyper genotype_sv /genphyse/cytogen/PacificBee/FromNcbi/GCA_003314205.1_INRA_AMelMel_1.0_genomic.fna \
            /home/gencel/vignal/work/DepthChr3Inverstion/allInsDel/HAv3_on_AMel/variants.vcf.gz \
            --sams=/home/gencel/vignal/work/DepthChr3Inverstion/allBAMsAMelMel.list \
            --region="CM010328.1:663884-665903"
```

* This time, both Indels were succesfully genotyped
* Comment:
  - Due to the change in reference? But the remapping by graphtyper is done on a graph, which is supposed to be neutral in this point of view.
  - Due to the use of bam files from alignments to AMelMel instead of HAv3? If only reads pre-mapped to the region are re-aligned locally, this could make a difference.
  - The bam files must correspond tyo the reference, so we can't test by swapping the reference and not the bams.
  - When comparing the two results on chromosome 2, missing genotypes in one instance can be successful in the other. All cases in which genotypes were successful in both instances gave the same results.

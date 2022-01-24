# ﻿Genetic map scripts
Sonia Eynard 12/08/2019

* All scripts are run throuh:

```bash
Pipeline_vcftogeno.sh
```

## 1. Preparation of genotype files
### 1.1 PacBio

* **vcftogeno.sh**: Extract genotypes from the vcf  with
  - **Output**: genotype.txt and genotypemap.txt

* **genomade.py**: Edit genotype file to get markers names as contig position and individual names
  - **Output**: genotype_f.txt and genotypemap_f.txt

* **genoedit.py**: Perform quality control on genotypes: sequencing errors using re-sequenced individuals, markers more than bi-allelic or for which one allele version does not exist in the queen, drone heterozygosity, queen homozygosity.
  - **Output**: summary files (out_genoedit.txt and marker_aftercontrol.txt), genotype file for each colony genotype_colony\*.txt and genotype file for each colony only for markers that are present in all colonies genotype_all_colony\*.txt

* **ContigsOnChromosomes_complete.r**: combining all prior information available for each contig, such as: size, a priori position, orientation, chromosome affiliation
  - **Output**: ContigsOnChromosomes_complete.csv

* **queen_cover.sh**: Extract genotyping coverage for each queen, for each chromosomes
  - **Output**: samCovQueenx.txt and samCovQueenx_c_sub.txt

### 1.2 Liu

* **info_all-colony.py**: Extract colony genotypes from markers in all colonies
  - **Output**: genotype.txt and genotypemap.txt

* **ContigsOnChromosomes_complete.r**: combining all prior information available for each contig, such as: size, a priori position, orientation, chromosome affiliation
  - **Output**: ContigsOnChromosomes_complete.csv

## 2. Contig orientation and ordering
### 2.1 PacBio and Liu (same procedure) for union of 3 colonies (1) and intersect of 3 colonies (2)

* **mk_bin.py**: script to create bins of unique vector of 0 and 1 for each contig. Assuming that the marker order is valid we create bins of multiple successive markers for which no allele change is observe across the individuals of the population, basically creating bins of recombinants.
  - **Output**: out_get_mk1.txt summarise the information extracted from the script and binx.txt 1 file per contig containing the bins and additional information about them (number of markers, list of markers, size in bp for each bin

* **mk_missing_bin.py**: script to create bins of unique vectors along each contig for markers shared by the 3 colonies only and identify missing contigs across the 3 colonies.
  - **Output**: out_get_mk_missingcontig.txt

* **bin_correction.r**: script to correct bins: repeated bins with bin in between smaller than the decided threshold for NCO/double recombinants events merged
  - **Output**: out_bin_description1.txt summary of analysis and binx.txt

* **bin_phasing.r**: script to ‘phase’ bins for each contig and order them following prior information on their chromosome affiliation. When not possible to place and order bins list the contig as missing.
  - **Output**: bin_description_colonyc.txt summary statistics for each colony and chromsomez.txt,
* **chr_phase.r**: in the case of full chromosomes no missing bins exist, however we still want to phase vectors of chromosome.
  - **Output**: all_chromosome.txt

* **missing_bins.r**: script to place and order bins for which we did not have prior information on position or that were not genotyped in all 3 colonies.
  - **Output**: out_missing1.txt summary file, all_chr_new.txt

* **chromosome_info_contig.r**: script sumarizing information on the contigs and chromosomes ones phase, order …
  - **Output**: out_chr1.txt

* **contigonchromosome.r**: script to order contigs on chromosome using prior knowledge e.g; ch1= contig10(+)-contg11(-)-contig12(+)-contig13(+)
  - **Output**: all_chr_new2.txt and combined_info.txt summary file compiling all knowledge

* **geneticmapbuild.r**: script to construct chromosome_map and genetic_map for recombinaison map calculation. Chromosome_map gives chromosome, marker and position information for all markers present in the genotype files of the 3 colonies (markers in common and exclusive). Genetic_map contains position of all the recombination events happening across the 3 colonies and the 43 individuals.
  - **Output**: chromosome_map_colony.txt, genetic_map_colony.txt

* **recombination_map.py**: script from Petit et al. 2017
Output: Mb_map.txt and parent_recombination.txt

* **mbmap.r**: script to plot the recombination rates
  - **Output**: recomb_rate.pdf

* **nb_recomdrone.r**: summary script for recombination events, outputting detail recombination events position, individual, size ...
  - **Output**: info_all.txt, recombintion_segments.txt, intercross_segments.txt

### 2.2 PacBio and Liu (same procedure) for 3 colonies combined in 1(3)

* **combine_colony.r**: script to combine genotypes for 3 colonies into 1.
  - **Output**: genotype_phased.txt

* **colony_phased.py** (PacBio): script to create phase bins for 3 colonies in 1. Using data from runs for intersect of colonies
  - **Output**: bins2.txt

* **order_contig_phased.r** (PacBio): script to order bins for 3 colonies in 1. Using bins2.txt
  - **Output**: all_chr_new2.txt and contigonchromsomes.txt

## 3. Comparison PacBio/Liu recombination profil and summary information (queen coverage, GC content, SNP density ….)

* **plot.r**: script to plot comparison between Liu and PacBio recombination profile
  - **Output**: PacBiovsLiu.pdf

# 4. Details on new assembly
### 4.1 Recombination hotspots

* **recom_hotspot.r**: script to localise recombination hotspot and check for repeated sequence on the genome at this positions
Output: hotspot.pdf and hotspot.txt

### 4.2 BUSCO
### 4.3 mapping quality
* **qual.sh**: script to measure mapping quality for multiple individuals, across genome and within windows on chromosomes
  - **Output**: mapping_qual.out

* **qual_plot.r**: script to plot mapping quality
  - **Output**: ind_mapq_chr.pdf

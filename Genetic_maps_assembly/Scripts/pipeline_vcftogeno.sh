#!/bin/bash
########################################################
#pipeline to extract and format genotype files from vcf#
########################################################
#files needed: 
#vcf + sample name 
#file containing info on contig size 
#files containing info on chromosome (chromosome contigid average_position orientation)
#genotypes from Liu et al. 2015
# modules: python3, R, augustus, busco3, hmmer, ncbi-blast, samtools)
#######################################################
initial_folder=vcf10
nb_colony=3
contig_size=contigSizes
contig_chr1=ContigsOnChromosomes.txt
contigorient=ChrOrientationAlain.csv
contig_chro=ContigsOnChromosomes_complete.csv

########## prep genotypes ########## 
### PacBio ###
parent_folder=PacBio
vcf=LiuAllGenotypesGVCF_snps_OKNoMissing.recode.vcf #vcf
id=SampleNames.txt #summary information of Samples ID, coverage, run, sex and colony
initialgeno=genotype.txt 	 
initialgenomap=genotypemap.txt
finalgeno=genotype_f.txt
finalgenomap=genotypemap_f.txt
name_colony='A,C,E'
name_ind_resequenced='B'
pattern_c='tig' 

cd $parent_folder
vcftogeno.sh $vcf $initialgeno $initialgenomap #get genotypes from vcf
genomade.py $id $initialgeno $initialgenomap $finalgeno $finalgenomap $pattern_c #edit the genotype file 
genoedit.py $finalgeno $nb_colony $name_colony $name_ind_resequenced > out_genoedit.txt #perform quality control on genotypes
rm $initialgeno $initialgenomap
ContigsOnChromosomes_complete.r $nb_colony $contig_size $contig_chr1 $contigorient $pattern_c genotype_colony #output: ContigsOnChromosomes_complete.csv combining all initial information available for each contig, such as: size, a priori position, orientation, chromosome affiliation
qsub -N queen_cover queen_cover.sh #extract information on sequence coverage for queens of each colony
cd ../

### Liu ###
parent_folder=Liu
name_colony_i=(I II III)
pattern_c='chr'

cd $parent_folder
for c in ${name_colony_i[*]}
	do
	echo $c
	sed 's/ *$//' Colony$c.genotype > colony$c.txt #get colony genotype files from complete genotype file Liu et al. 2015
	done
rm Colony*
info_all_colony.py colony #get colony genotype files for markers in all colonies
mv colonyI.txt colony_colony1.txt
mv colonyII.txt colony_colony2.txt
mv colonyIII.txt colony_colony3.txt
ContigsOnChromosomes_complete.r $nb_colony $contig_size $contig_chr1 $contigorient $pattern_c colony_colony
cd ../

########## Contig orientation and ordering genotypes ##########
size_threshold=2000 #threshold choosen to distinguish between CO and NCO/double recombinant (based on literature and test in Liu et al 2015)
### PacBio ###
##### 1 ##### each colony seperately and combined afterwards (union of 3 colonies)
parent_folder=PacBio
pattern_c='tig' 
pattern='group_colony'

cd $parent_folder
for c in $(seq 1 $nb_colony)
	do
	mkdir -p $pattern$c
	done
genotype_bin='genotype'
mk_bin.py $nb_colony $genotype_bin $pattern > out_get_mk1.txt # create bins of unique vectors along each contig
for c in $(seq 1 $nb_colony)
	do
	cd $pattern$c/
	mkdir -p bin_initial
	mv bins*.txt bin_initial
	cd ../
	done
bin_correction.r $nb_colony $size_threshold $pattern > out_bin_description1.txt #correct bins if size below threshold for NCO/double recombinant events by combining bins  
for c in $(seq 1 $nb_colony)
	do
	cd $pattern$c/
	mkdir -p chromosome
	cd ../
	done
if [[ "$pattern_c" == "tig" ]]
	then 
	echo contig #in the case where contigs are available
	bin_phasing.r $nb_colony $contig_chro $pattern > out_phasing1.txt #phase and order bins when possible based on prior knowledge
	missing_bins.r $nb_colony $contig_chro $pattern > out_missing1.txt #complete bins association to chromosomes for bins that are colony specific or that were not affiliated to chromosome
	chromosome_info_contig.r $nb_colony $pattern > out_chr1.txt # summarise information on chromosomes and contigs
else
	echo chromosome #in the case where chromosomes are available
	chr_phase.r $nb_colony $pattern #phase chromosomes
	chromosome_info_chr.r $nb_colony $pattern > out_chr1.txt 
fi
contigonchromosome.r $nb_colony $contig_chro $pattern #combined_info.txt summarise information on contigs, orientation, order and chromosome allocation
mv combined_info.txt combined_info_colony.txt
mv bin_description_colony1.txt bin_description_colony11.txt
mv bin_description_colony2.txt bin_description_colony21.txt
mv bin_description_colony3.txt bin_description_colony31.txt
mv bin_description.pdf bin_description1.pdf
genotype_bin='genotype_colony'
geneticmapbuild.r $nb_colony $pattern combined_info_colony.txt $genotype_bin #create chromosome_map and genetic_map necessary for recombination_map
mv genetic_map.txt genetic_map_colony.txt
mv chromosome_map.txt chromosome_map_colony.txt

recombination_map.py chromosome_map_colony.txt genetic_map_colony.txt 1000000 False 0.23 0.01 #inference of recombination rate along the sequence. The last 3 arguments are necessary to calibrate prior: if False then provide values for alpha and beta for gamma
mv Mb_map.txt Mb_map1_colony.txt
mbmap.r Mb_map1_colony.txt chromosome_map_colony.txt genetic_map_colony.txt #descriptive plots of the recombination map created
recombination_map.py chromosome_map_colony.txt genetic_map_colony.txt 500000 False 0.23 0.01
mv Mb_map.txt Mb_map0.5_colony.txt
mbmap.r Mb_map0.5_colony.txt chromosome_map_colony.txt genetic_map_colony.txt
mv parent_recombination.txt parent_recombination_colony.txt
mv recomb_rate1.pdf recomb_rate1_colony.pdf
mv recomb_rate0.5.pdf recomb_rate0.5_colony.pdf
nb_recomdrone.r $nb_colony $size_threshold $pattern combined_info_colony.txt > out_recombination1.txt #output of the number of recombination per individual, per chromosome, estimation of recombination rate cM/Mb
mv info_all1.txt info_all11.txt
mv info_all2.txt info_all21.txt
mv info_all3.txt info_all31.txt
mv intercross_segments.txt intercross_segments1.txt
mv recombination_segments.txt recombination_segments1.txt
mv segment.pdf segment1.pdf
cd ../

##### 2 ##### intersect of 3 colonies
parent_folder=PacBio
pattern_c='tig'
pattern='group_colony_all'
cd $parent_folder
for c in $(seq 1 $nb_colony)
	do
	mkdir -p $pattern$c
	done
genotype_bin='genotype_all'
mk_bin.py $nb_colony $genotype_bin $pattern > out_get_mk2.txt
for c in $(seq 1 $nb_colony)
	do
	cd $pattern$c/
	mkdir -p bin_initial
	mv bins*.txt bin_initial
	cd ../
	done
genotype_bin='genotype'
mk_missing_bin.py $nb_colony $genotype_bin $pattern > out_get_mk_missingcontig.txt # create bins of unique vectors along each contig for markers shared by the 3 colonies only
bin_correction.r $nb_colony $size_threshold $pattern > out_bin_description2.txt
genotype_bin='genotype_all'
for c in $(seq 1 $nb_colony)
	do
	cd $pattern$c/
	mkdir -p chromosome
	cd ../
	done
if [[ "$pattern_c" == "tig" ]]
	then 
	echo contig
	bin_phasing.r $nb_colony $contig_chro $pattern > out_phasing2.txt
	missing_bins.r $nb_colony $contig_chro $pattern > out_missing2.txt
	chromosome_info_contig.r $nb_colony $pattern > out_chr2.txt
else
	echo chromosome
	chr_phase.r $nb_colony $pattern
	chromosome_info_chr.r $nb_colony $pattern > out_chr2.txt
fi
contigonchromosome.r $nb_colony $contig_chro $pattern
mv combined_info.txt combined_info_all.txt
mv bin_description_colony1.txt bin_description_colony12.txt
mv bin_description_colony2.txt bin_description_colony22.txt
mv bin_description_colony3.txt bin_description_colony32.txt
mv bin_description.pdf bin_description2.pdf
genotype_bin='genotype_all_colony'
geneticmapbuild.r $nb_colony $pattern combined_info_all.txt $genotype_bin
mv genetic_map.txt genetic_map_all.txt
mv chromosome_map.txt chromosome_map_all.txt

recombination_map.py chromosome_map_all.txt genetic_map_all.txt 1000000 False 0.23 0.01 
mv Mb_map.txt Mb_map1_all.txt
mbmap.r Mb_map1_all.txt chromosome_map_all.txt genetic_map_all.txt
recombination_map.py chromosome_map_all.txt genetic_map_all.txt 500000 False 0.23 0.01
mv Mb_map.txt Mb_map0.5_all.txt
mbmap.r Mb_map0.5_all.txt chromosome_map_all.txt genetic_map_all.txt
mv parent_recombination.txt parent_recombination_all.txt
mv recomb_rate1.pdf recomb_rate1_all.pdf
mv recomb_rate0.5.pdf recomb_rate0.5_all.pdf
nb_recomdrone.r $nb_colony $size_threshold $pattern combined_info_all.txt > out_recombination2.txt
mv info_all1.txt info_all12.txt
mv info_all2.txt info_all22.txt
mv info_all3.txt info_all32.txt
mv intercross_segments.txt intercross_segments2.txt
mv recombination_segments.txt recombination_segments2.txt
mv segment.pdf segment2.pdf
cd ../

##### 3 ##### 3 colonies as 1 colony
parent_folder=PacBio
genotype_bin='genotype_all'
geno_to_reorder='genotype_allphased.txt'
pattern='group_colony_all'
cd $parent_folder
combine_colony.r $nb_colony $genotype_bin $pattern combined_info_all.txt # combine phased vectors of the 3 colonies into 1 
mkdir bins2
colony_phased.py $geno_to_reorder # create bins of unique vectors along each contig for markers when 3 colonies combined in 1
cd bins2
mkdir -p bin_initial
mv bins*.txt bin_initial
cd ../
nb_colony=1
pattern='bins2'
bin_correction.r $nb_colony $size_threshold $pattern > out_bin_description3.txt
mkdir bins2/chromosome
order_contig_phased.r $contig_chro # order contigs when 3 colonies combined in 1
contigonchromosome.r $nb_colony $contig_chro $pattern
cd ../

### Liu ###
##### 1 ##### each colony seperately and combined afterwards (union of 3 colonies)
parent_folder=Liu
pattern_c='chr'
pattern='colony_group'
cd $parent_folder
for c in $(seq 1 $nb_colony)
	do
	mkdir -p $pattern$c
	done
genotype_bin='colony'
mk_bin.py $nb_colony $genotype_bin $pattern > out_get_mk1.txt
for c in $(seq 1 $nb_colony)
	do
	cd $pattern$c/
	mkdir -p bin_initial
	mv bins*.txt bin_initial
	cd ../
	done
bin_correction.r $nb_colony $size_threshold $pattern > out_bin_description1.txt
for c in $(seq 1 $nb_colony)
	do
	cd $pattern$c/
	mkdir -p chromosome
	cd ../
	done
if [[ "$pattern_c" == "tig" ]]
	then 
	echo contig
	bin_phasing.r $nb_colony $contig_chro $pattern > out_phasing1.txt
	missing_bins.r $nb_colony $contig_chro $pattern > out_missing1.txt
	chromosome_info_contig.r $nb_colony $pattern > out_chr1.txt
else
	echo chromosome
	chr_phase.r $nb_colony $pattern
	chromosome_info_chr.r $nb_colony $pattern > out_chr1.txt
fi
contigonchromosome.r $nb_colony $contig_chro $pattern
mv combined_info.txt combined_info_colony.txt
awk -F',' '{$38="+" ; print}' OFS=',' combined_info_colony.txt > combined_info_colony2.txt #specific of chr values, no orientation needed, always (+)
awk -F',' '{$39="+" ; print}' OFS=',' combined_info_colony2.txt > combined_info_colony.txt
awk -F',' 'NR==1{$38="orientation_consensus"}1' OFS=',' combined_info_colony.txt > combined_info_colony2.txt
awk -F',' 'NR==1{$39="orientation_final"}1' OFS=',' combined_info_colony2.txt > combined_info_colony.txt
rm combined_info_colony2.txt
mv bin_description_colony1.txt bin_description_colony11.txt
mv bin_description_colony2.txt bin_description_colony21.txt
mv bin_description_colony3.txt bin_description_colony31.txt
mv bin_description.pdf bin_description1.pdf

genotype_bin='colony_colony'
geneticmapbuild.r $nb_colony $pattern combined_info_colony.txt $genotype_bin
mv genetic_map.txt genetic_map_colony.txt
mv chromosome_map.txt chromosome_map_colony.txt

recombination_map.py chromosome_map_colony.txt genetic_map_colony.txt 1000000 False 0.37 0.01 
mv Mb_map.txt Mb_map1_colony.txt
mbmap.r Mb_map1_colony.txt chromosome_map_colony.txt genetic_map_colony.txt
recombination_map.py chromosome_map_colony.txt genetic_map_colony.txt 500000 False 0.37 0.01
mv Mb_map.txt Mb_map0.5_colony.txt
mbmap.r Mb_map0.5_colony.txt chromosome_map_colony.txt genetic_map_colony.txt
mv parent_recombination.txt parent_recombination_colony.txt
mv recomb_rate1.pdf recomb_rate1_colony.pdf
mv recomb_rate0.5.pdf recomb_rate0.5_colony.pdf
nb_recomdrone.r $nb_colony $size_threshold $patternl combined_info_colony.txt > out_recombination1.txt
mv info_all1.txt info_all11.txt
mv info_all2.txt info_all21.txt
mv info_all3.txt info_all31.txt
mv intercross_segments.txt intercross_segments1.txt
mv recombination_segments.txt recombination_segments1.txt
mv segment.pdf segment1.pdf
cd ../

##### 2 ##### intersect of 3 colonies
parent_folder=Liu
pattern_c='chr'
pattern='colony_groupall'
cd $parent_folder
for c in $(seq 1 $nb_colony)
	do
	mkdir -p $pattern$c
	done
genotype_bin='colony_all'
mk_bin.py $nb_colony $genotype_bin $pattern > out_get_mk2.txt
for c in $(seq 1 $nb_colony)
	do
	cd $pattern$c/
	mkdir -p bin_initial
	mv bins*.txt bin_initial
	cd ../
	done
bin_correction.r $nb_colony $size_threshold $pattern > out_bin_description2.txt
for c in $(seq 1 $nb_colony)
	do
	cd $pattern$c/
	mkdir -p chromosome
	cd ../
	done
if [[ "$pattern_c" == "tig" ]]
	then 
	echo contig
	bin_phasing.r $nb_colony $contig_chro $pattern > out_phasing2.txt
	missing_bins.r $nb_colony $contig_chro $pattern > out_missing2.txt
	chromosome_info_contig.r $nb_colony $pattern > out_chr2.txt
else
	echo chromosome
	chr_phase.r $nb_colony $pattern
	chromosome_info_chr.r $nb_colony $pattern > out_chr2.txt
fi
contigonchromosome.r $nb_colony $contig_chro $pattern
mv combined_info.txt combined_info_all.txt
awk -F',' '{$38="+" ; print}' OFS=',' combined_info_all.txt > combined_info_all2.txt
awk -F',' '{$39="+" ; print}' OFS=',' combined_info_all2.txt > combined_info_all.txt
awk -F',' 'NR==1{$38="orientation_consensus"}1' OFS=',' combined_info_all.txt > combined_info_all2.txt
awk -F',' 'NR==1{$39="orientation_final"}1' OFS=',' combined_info_all2.txt > combined_info_all.txt
rm combined_info_all2.txt
mv bin_description_colony1.txt bin_description_colony12.txt
mv bin_description_colony2.txt bin_description_colony22.txt
mv bin_description_colony3.txt bin_description_colony32.txt
mv bin_description.pdf bin_description2.pdf

genotype_bin='colony_all_colony'
geneticmapbuild.r $nb_colony $pattern combined_info_all.txt $genotype_bin
mv genetic_map.txt genetic_map_all.txt
mv chromosome_map.txt chromosome_map_all.txt

recombination_map.py chromosome_map_all.txt genetic_map_all.txt 1000000 False 0.37 0.01
mv Mb_map.txt Mb_map1_all.txt
mbmap.r Mb_map1_all.txt chromosome_map_all.txt genetic_map_all.txt
recombination_map.py chromosome_map_all.txt genetic_map_all.txt 500000 False 0.37 0.01
mv Mb_map.txt Mb_map0.5_all.txt
mbmap.r Mb_map0.5_all.txt chromosome_map_all.txt genetic_map_all.txt
mv parent_recombination.txt parent_recombination_colony.txt
mv recomb_rate1.pdf recomb_rate1_all.pdf
mv recomb_rate0.5.pdf recomb_rate0.5_all.pdf
nb_recomdrone.r $nb_colony $size_threshold $pattern combined_info_all.txt > out_recombination2.txt
mv info_all1.txt info_all12.txt
mv info_all2.txt info_all22.txt
mv info_all3.txt info_all32.txt
mv intercross_segments.txt intercross_segments2.txt
mv recombination_segments.txt recombination_segments2.txt
mv segment.pdf segment2.pdf
cd ../

##### 3 ##### 3 colonies as 1 colony
parent_folder=Liu
genotype_bin='colony_all'
pattern='colony_groupall'
nb_colony=3
cd $parent_folder
combine_colony.r $nb_colony $genotype_bin $pattern combined_info_all.txt
cd ../

########## comparison PacBio Liu ##########
nb_colony=3
kb=2000
parent_folder=PacBio
pattern='group_colony'
parent_folderl=Liu
patternl='colony_group'
plot.r $nb_colony $parent_folder $pattern $parent_folderl $patternl 'Mb_map1_' combined_info_colony.txt 1000000 $kb
mv PacBiovsLiu_all.pdf PacBiovsLiu_all1.pdf
mv PacBiovsLiu_colony.pdf PacBiovsLiu_colony1.pdf
mv PacBiovsLiu_chr*.png PacBiovsLiu_chr*_colony1.png
mv PacBiovsLiu_chr3Paper.png PacBiovsLiu_chr3Paper_colony1.png
plot.r $nb_colony $parent_folder $pattern $parent_folderl $patternl 'Mb_map0.5_' combined_info_colony.txt 500000 $kb
mv PacBiovsLiu_all.pdf PacBiovsLiu_all05.pdf
mv PacBiovsLiu_colony.pdf PacBiovsLiu_colony05.pdf
mv PacBiovsLiu_chr*.png PacBiovsLiu_chr*_colony0.5.png
mv PacBiovsLiu_chr3Paper.png PacBiovsLiu_chr3Paper_colony0.5.png

########## new assembly details ##########
### recombination hotspots ###
recom_hotspot.r # localise recombination hotspot

### BUSCO ###
target='/genphyse/dynagen/BeeStrong/Fasta/GCF_000002195.4_Amel_4.5_genomic.fna'
species='/usr/local/bioinfo/src/BUSCO/datasets/metazoa_odb9'
out='busco_Amel_4.5_metazoa3'
python  /usr/local/bioinfo/src/BUSCO/busco-3.0.2/scripts/run_BUSCO.py --in '$target' --out '$out' --lineage '$species' --mode genome

target='/genphyse/dynagen/BeeStrong/Fasta/GCA_003254395.1_Amel_HAv3_genomic.fna'
species='/usr/local/bioinfo/src/BUSCO/datasets/metazoa_odb9'
out='busco_Amel_HAv3_metazoa3'
python  /usr/local/bioinfo/src/BUSCO/busco-3.0.2/scripts/run_BUSCO.py --in '$target' --out '$out' --lineage '$species' --mode genome

target='/genphyse/dynagen/BeeStrong/Fasta/GCA_003254395.2_Amel_HAv3.1_genomic.fna'
species='/usr/local/bioinfo/src/BUSCO/datasets/metazoa_odb9'
out='busco_Amel_HAv3.1_metazoa3'
python  /usr/local/bioinfo/src/BUSCO/busco-3.0.2/scripts/run_BUSCO.py --in '$target' --out '$out' --lineage '$species' --mode genome

target='/genphyse/dynagen/BeeStrong/Fasta/GCA_003314205.1_INRA_AMelMel_1.0_genomic.fna'
species='/usr/local/bioinfo/src/BUSCO/datasets/metazoa_odb9'
out='busco_AMelMel_metazoa3'
python  /usr/local/bioinfo/src/BUSCO/busco-3.0.2/scripts/run_BUSCO.py --in '$target' --out '$out' --lineage '$species' --mode genome

target='/genphyse/dynagen/BeeStrong/Fasta/AMelMel_1.vcftools.fasta'
species='/usr/local/bioinfo/src/BUSCO/datasets/metazoa_odb9'
out='busco_AMelMel_1_metazoa3'
python  /usr/local/bioinfo/src/BUSCO/busco-3.0.2/scripts/run_BUSCO.py --in '$target' --out '$out' --lineage '$species' --mode genome

target='/genphyse/dynagen/BeeStrong/Fasta/GCF_000002195.4_Amel_4.5_genomic.fna'
species='/usr/local/bioinfo/src/BUSCO/datasets/hymenoptera_odb9'
out='busco_Amel_4.5_hymenoptera3'
python  /usr/local/bioinfo/src/BUSCO/busco-3.0.2/scripts/run_BUSCO.py --in '$target' --out '$out' --lineage '$species' --mode genome

target='/genphyse/dynagen/BeeStrong/Fasta/GCA_003254395.1_Amel_HAv3_genomic.fna'
species='/usr/local/bioinfo/src/BUSCO/datasets/hymenoptera_odb9'
out='busco_Amel_HAv3_hymenoptera3'
python  /usr/local/bioinfo/src/BUSCO/busco-3.0.2/scripts/run_BUSCO.py --in '$target' --out '$out' --lineage '$species' --mode genome

target='/genphyse/dynagen/BeeStrong/Fasta/GCA_003254395.2_Amel_HAv3.1_genomic.fna'
species='/usr/local/bioinfo/src/BUSCO/datasets/hymenoptera_odb9'
out='busco_Amel_HAv3.1_hymenoptera3'
python  /usr/local/bioinfo/src/BUSCO/busco-3.0.2/scripts/run_BUSCO.py --in '$target' --out '$out' --lineage '$species' --mode genome

target='/genphyse/dynagen/BeeStrong/Fasta/GCA_003314205.1_INRA_AMelMel_1.0_genomic.fna'
species='/usr/local/bioinfo/src/BUSCO/datasets/hymenoptera_odb9'
out='busco_AMelMel_hymenoptera3'
python  /usr/local/bioinfo/src/BUSCO/busco-3.0.2/scripts/run_BUSCO.py --in '$target' --out '$out' --lineage '$species' --mode genome

target='/genphyse/dynagen/BeeStrong/Fasta/AMelMel_1.vcftools.fasta'
species='/usr/local/bioinfo/src/BUSCO/datasets/hymenoptera_odb9'
out='busco_AMelMel_1_hymenoptera3'
python  /usr/local/bioinfo/src/BUSCO/busco-3.0.2/scripts/run_BUSCO.py --in '$target' --out '$out' --lineage '$species' --mode genome

### mapping quality ###
sh qual.sh
qual_plot.r

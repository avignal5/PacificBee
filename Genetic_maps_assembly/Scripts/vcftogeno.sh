################################################
### script to get genotypes from VCF
################################################
#!/bin/bash

# creation of a genotype file 
bcftools query -f '%CHROM %POS[ %GT]\n' $1 > geno1.txt
id=`grep '#CHROM' $1`
idsub=`echo "$id" | grep -o 'SRR.*'`
topline='contigID position '
topline+=$idsub
echo 'genotype file printed'
sed '1i\'"$topline" geno1.txt > $2
rm geno1.txt

bcftools query -f '%CHROM %POS %REF %ALT\n' $1 > genotype_map1.txt
sed '1i contigID position ref alt' genotype_map1.txt > $3
rm genotype_map1.txt
echo 'genotype map file printed'

################################################
### script to extract queen genotyping coverage
################################################
#!/bin/bash

nb_colony=3
for c in $(seq 1 $nb_colony)
do
echo $c
for z in $(seq 1 16)
do
echo $z
n=chr
n+=$z
awk -v var=$n '$1==var' samCovQueen${c}.txt > samCovQueen${c}_${z}_sub.txt 
done
done

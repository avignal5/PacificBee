################################################
### script to measure mapping quality for a set of individuals on average and on genome windows
################################################
#!/bin/bash

bam_files_amelmel='/genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnAMelMel/*/mapping/'
arr_amelmel=(${bam_files_amelmel// /})  
for i in ${arr_amelmel[@]}
do 
id=$(echo ${i}|cut -d'/' -f 8)
echo ${id}
echo 'AMelMel'
b_amelmel='/genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnAMelMel/'${id}'/mapping/'
samtools flagstat ${b_amelmel}*.bam
IFS=$'\n'
chr=($(samtools view -H  ${b_amelmel}*.bam |grep 'CM'))
for j in ${chr[@]}
do 
chr_name=$(echo ${j}|cut -f2 | cut -d':' -f2)
echo ${chr_name}
chr_size=$(echo ${j}|cut -f3 | cut -d':' -f2)
echo ${chr_size}
n=1
m=1000000
while [ ${n} -lt ${chr_size} ]
do 
echo ${n} ${m} 
samtools view ${b_amelmel}*.bam ${chr_name}:${n}-${m} | awk '{sum+=$5} END { print "Mean MAPQ =",sum/NR}'
n=$((n+1000000))
m=$((m+1000000))
done
done

echo 'HAV3'
b_hav3='/genphyse/cytogen/seqapipop/Data/Apis-mellifera/seqapipopOnHAV3/'${id}'/mapping/'
samtools flagstat ${b_hav3}*.bam
IFS=$'\n'
chr=($(samtools view -H  ${b_hav3}*.bam |grep 'CM'))
for j in ${chr[@]}
do 
chr_name=$(echo ${j}|cut -f2 | cut -d':' -f2)
echo ${chr_name}
chr_size=$(echo ${j}|cut -f3 | cut -d':' -f2)
echo ${chr_size}
n=1
m=1000000
while [ ${n} -lt ${chr_size} ]
do 
echo ${n} ${m} 
samtools view ${b_hav3}*.bam ${chr_name}:${n}-${m} | awk '{sum+=$5} END { print "Mean MAPQ =",sum/NR}'
n=$((n+1000000))
m=$((m+1000000))
done
done
done







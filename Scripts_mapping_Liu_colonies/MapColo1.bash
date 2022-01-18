#!/bin/bash

#------------------------------------------------------------------
#EDIT FOR ACCESS TO FOLDER WITH THE SEQUENCE FILES
IN="/work/project/cytogen/Alain/LiuRecombination/FromNCBI/Colony1"
#------------------------------------------------------------------

# The following lines will recover an array of SAMPLE names:
# Read in list of files
SAMPLELIST=(${IN}/*fastq.gz)
# Strip out the suffix
tmp1=(${SAMPLELIST[@]/_R1.fastq.gz/})
tmp2=(${tmp1[@]/_R2.fastq.gz/})
# Strip out the prefix
tmp3=(${tmp2[@]/${IN}\//})
# Reduce to unique list
tmp4=( $( printf "%s\n" "${tmp3[@]}" | awk 'x[$0]++ == 0' ) )

#----------------------------------------------------------------
#EDIT FOR WRITING RESULTS
# This is the path to where the resulting files will be dumped
DUMP="/work/project/cytogen/Alain/LiuRecombination/Alignments/Colony1"
#----------------------------------------------------------------

# This is the path to pipeline scripts,
PIPE="/genphyse/cytogen/seqapipop/ScriptsAV"

#----------------------------------------------------------------
# This is the loop to send each sample to the cluster to be mapped
for ID in ${tmp4[@]:0:${#tmp4[@]}}
do
  mkdir -p ${DUMP}/logs/${ID}
  cd ${DUMP}
  qsub -q workq -l mem=4G -l h_vmem=48G -pe parallel_smp 4 \
    -o ${DUMP}/logs/${ID} \
    -e ${DUMP}/logs/${ID} \
    ${PIPE}/mapLiu.bash  -s ${ID} -i ${IN} -o ${DUMP}
done
#-----------------------------------------------------------------


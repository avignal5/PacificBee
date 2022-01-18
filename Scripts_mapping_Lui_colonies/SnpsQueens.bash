#!/bin/bash

#------------------------------------------------------------------
#EDIT FOR ACCESS TO SEQUENCE FILES AND BAM FILES
DATASET=Queens

#------------------------------------------------------------------
# List of fastqs: 
IN=/work/project/cytogen/Alain/LiuRecombination/FromNCBI/${DATASET}
# List of bams
IN2=/work/project/cytogen/Alain/LiuRecombination/Alignments/${DATASET}
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
DUMP=/work/project/cytogen/Alain/LiuRecombination/AlignmentsHaplotypeCaller/${DATASET}
#----------------------------------------------------------------

# This is the path to the pipeline scripts,
PIPE="/genphyse/cytogen/seqapipop/ScriptsAV"

#----------------------------------------------------------------
# This is the loop to send each sample to the cluster to be mapped
for ID in ${tmp4[@]:0:${#tmp4[@]}}
do
  mkdir -p ${DUMP}/logs/${ID}
  cd ${DUMP}
  qsub -q workq -l mem=30G -l h_vmem=55G \
    -o ${DUMP}/logs/${ID} \
    -e ${DUMP}/logs/${ID} \
    ${PIPE}/snpsLiuOK.bash  -s ${ID} -o ${DUMP} -i ${IN2}
done
#-----------------------------------------------------------------


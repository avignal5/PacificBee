#!/bin/bash
#$ -M alain.vignal@inra.fr
#$ -m a

module load bioinfo/Java8

# set -e will cause script to terminate on an error, set +e allows it to continue
set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

# ==============================================================================
# map.sh
# ==============================================================================
echo -e "${BLUE}====================================${NC}"
echo -e "${BLUE}SeqApiPop: genotyping${NC}"
echo -e "${BLUE}====================================${NC}"

function usage()
{
	echo "Usage:"
	echo "map.sh -s <sample name> -f <path to fastq.gz files> -o <output path> [options]"
}
# ==============================================================================
# Clear variables and set constant variables
# ==============================================================================
REF=/home/gencel/vignal/save/Genomes/Abeille/PacificBee/BlackBee_genome.AllAlternative.NoMultipleAllele.fasta
GATK=/usr/local/bioinfo/src/GATK/GenomeAnalysisTK-3.7
OUT=
SAMPLE=

# ==============================================================================
# Read in variables from optins
# ==============================================================================
while getopts ":i:s:o:" opt; do
  case $opt in
    s) SAMPLE=${OPTARG};;
    o) OUT=${OPTARG};;
    i) IN=${OPTARG};;
  esac
done

if [[ -z ${SAMPLE} ]] | [[ -z ${OUT} ]] | [[ -z ${IN} ]]
then
  usage
  exit 1
fi

# ==============================================================================
# Message
# ==============================================================================

echo -e "\n${RED}Variables:\n\nimported succesfully:${RED}"
echo -e "${GREEN}SAMPLE: ${SAMPLE}${NC}"
echo -e "${GREEN}REF: ${REF}${NC}"
echo -e "${GREEN}IN: ${IN}${NC}"
echo -e "${GREEN}OUT: ${OUT}${NC}"
echo -e "${GREEN}GATK: ${GATK}${NC}"

# ==============================================================================
# Create folders for storing files
# ==============================================================================
mkdir -p ${OUT}/${SAMPLE}/logs
mkdir -p ${OUT}/${SAMPLE}/metrics
mkdir -p ${OUT}/${SAMPLE}/vcfs

# ==============================================================================
# Genotype with Haplotype Caller
# ==============================================================================

echo -e "\n${RED}GATK:HaplotypeCaller${NC}:"
logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_HaplotypeCaller.err

java -Xmx30g -jar ${GATK}/GenomeAnalysisTK.jar \
    -T HaplotypeCaller \
    -R ${REF} \
    -I ${IN}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
    --genotyping_mode DISCOVERY \
    -stand_call_conf 10 \
    -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_HC.g.vcf \
    -variant_index_type LINEAR \
    -variant_index_parameter 128000 \
    -ERC GVCF
    2> >(tee "$logfile")


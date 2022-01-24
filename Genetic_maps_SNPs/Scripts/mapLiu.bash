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
echo -e "${BLUE}SeqApiPop: mapping${NC}"
echo -e "${BLUE}====================================${NC}"

function usage()
{
	echo "Usage:"
	echo "mapLiu.sh -s <sample name> -i <path to fastq.gz files> -o <output path> [options]"
}
# ==============================================================================
# Clear variables and set defaults for nThreads, BQSR bootstraps and RAM
# ==============================================================================
BOOTS=2
REF=/home/gencel/vignal/save/Genomes/Abeille/PacificBee/BlackBee_genome.AllAlternative.NoMultipleAllele.fasta
IN=
OUT=
SAMPLE=
BWA=/usr/local/bioinfo/src/bwa/bwa-0.7.15
PICARD=/usr/local/bioinfo/src/picard-tools/picard-tools-2.1.1
GATK=/usr/local/bioinfo/src/GATK/GenomeAnalysisTK-3.7

# ==============================================================================
# Read in variables from optins
# ==============================================================================
while getopts ":i:s:o:" opt; do
  case $opt in
    i) IN=${OPTARG};;
    s) SAMPLE=${OPTARG};;
    o) OUT=${OPTARG};;
  esac
done

if [[ -z ${SAMPLE} ]] | [[ -z ${IN} ]] | [[ -z ${OUT} ]]
then
  usage
  exit 1
fi

echo -e "\n\e[91mScript variables:\e[0m"
echo -e "\e[96mParameters file: \e[92m${PARAM}\e[0m"
echo -e "\e[96mBQSR iterations: \e[92m${BOOTS}\e[0m"


# ==============================================================================
#For BAM header (read in first line of FASTQ)
# ==============================================================================

header=`zcat ${IN}/${SAMPLE}_R1.fastq.gz | head -1 | cut -c 2-`	# removes @ from fastq
tmp1=${header% *}	# All before first space	HWI...1953
tmp2=${tmp1%:*}		# All before last :	HWI...1718
tmp3=${tmp2%:*}		# All before last :	HWI...1101
TECHNO=${tmp3%:*}	# All before last :	HWI...6

INDIVIDUAL=${SAMPLE%%_*}	# All before first _		Gri-w

tmp4=${SAMPLE%_*}	# All before last _		Gri-w_GTGAAA
TAG=${tmp4##*_}	# All after first _		GTGAAA

#FASTQ files
READ1=${IN}/${SAMPLE}_R1.fastq.gz
READ2=${IN}/${SAMPLE}_R2.fastq.gz

# ==============================================================================
#Message
# ==============================================================================

echo -e "\n${RED}Variables:\n\nimported succesfully:${RED}"
echo -e "${GREEN}SAMPLE: ${SAMPLE}${NC}"
echo -e "${GREEN}REF: ${REF}${NC}"
echo -e "${GREEN}IN: ${IN}${NC}"
echo -e "${GREEN}OUT: ${OUT}${NC}"
echo -e "${GREEN}BWA: ${BWA}${NC}"
echo -e "${GREEN}PICARD: ${PICARD}${NC}"
echo -e "${GREEN}GATK: ${GATK}${NC}"

echo -e "\n${RED}For BAM header:${RED}"
echo -e "${GREEN}TECHNO: ${TECHNO}${NC}"
echo -e "${GREEN}SAMPLE: ${INDIVIDUAL}${NC}"
echo -e "${GREEN}TAG: ${TAG}${NC}"

echo -e "\n${RED}FASTQ files:${NC}"
echo -e "Reads1: ${READ1}"
echo -e "Reads2: ${READ2}"

# ==============================================================================
# Create folders for storing files
# ==============================================================================
mkdir -p ${OUT}/${SAMPLE}/logs
mkdir -p ${OUT}/${SAMPLE}/metrics
mkdir -p ${OUT}/${SAMPLE}/vcfs

# ==============================================================================
# Number of steps in pipeline, this will need updating if additional steps are added, to provide an idea of how long is remaining whilst running the pipe
# ==============================================================================
STEPS=7

# ==============================================================================
# Mapping to Reference (BWA)
# ==============================================================================

echo -e "\n${RED}1/${STEPS}: Mapping reads (BWA MEM)${NC}"
##${BWA}/bwa mem -M -R @RG"\t"ID:${SAMPLE}"\t"SM:${INDIVIDUAL}"\t"PL:ILLUMINA"\t"LB:${TAG}"\t"PU:${TECHNO} ${REF} ${READ1} ${READ2} > ${OUT}/${SAMPLE}/${SAMPLE}.sam
${BWA}/bwa mem -M -R @RG"\t"ID:${SAMPLE}"\t"SM:${INDIVIDUAL}"\t"PL:ILLUMINA"\t"LB:${TAG} ${REF} ${READ1} ${READ2} > ${OUT}/${SAMPLE}/${SAMPLE}.sam


# ==============================================================================
# Sort, mark duplicates and index
# ==============================================================================

# Sort SAM into coordinate order and save as BAM
echo -e "\n${RED}2/${STEPS}: Sorting SAM and save as BAM: ${OUT}/${SAMPLE}/${SAMPLE}.bam${NC}"
logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_PicardSortSam.log
java -Xmx4g -jar $PICARD/picard.jar SortSam \
  INPUT=${OUT}/${SAMPLE}/${SAMPLE}.sam \
  OUTPUT=${OUT}/${SAMPLE}/${SAMPLE}.bam \
  SORT_ORDER=coordinate \
  QUIET=T \
  VERBOSITY=ERROR \
  VALIDATION_STRINGENCY=LENIENT \
  2> >(tee "$logfile")
rm ${OUT}/${SAMPLE}/${SAMPLE}.sam

# Mark duplicates 
echo -e "\n${RED}3/${STEPS}: Mark Duplicates${NC}"
logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_PicardMarkDuplicates.log
java -Xmx4g -jar  ${PICARD}/picard.jar MarkDuplicates \
  INPUT=${OUT}/${SAMPLE}/${SAMPLE}.bam \
  OUTPUT=${OUT}/${SAMPLE}/${SAMPLE}_dup.bam \
  METRICS_FILE=${OUT}/${SAMPLE}/metrics/${SAMPLE}_dup.metrics \
  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
  QUIET=T \
  VERBOSITY=ERROR \
  VALIDATION_STRINGENCY=LENIENT \
  2> >(tee "$logfile")
rm ${OUT}/${SAMPLE}/${SAMPLE}.bam

# Index BAM file 
echo -e "\n${RED}4/${STEPS}: BuildBamIndex${NC}"
logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_PicardBuildBamIndex.log
java -Xmx4g -jar ${PICARD}/picard.jar BuildBamIndex \
  INPUT=${OUT}/${SAMPLE}/${SAMPLE}_dup.bam \
  QUIET=T \
  VERBOSITY=ERROR \
  2> >(tee "$logfile")



# ==============================================================================
# Local realignment around INDELs
# ==============================================================================

# Create target interval list
echo -e "\n${RED}5/${STEPS}: GATK:RealignerTargetCreator${NC}"
logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_RealignerTargetCreator.err

java -Xmx4g -jar ${GATK}/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R ${REF} \
    -I ${OUT}/${SAMPLE}/${SAMPLE}_dup.bam \
    -o ${OUT}/${SAMPLE}/logs/${SAMPLE}_dup_intervals.list \
    -l FATAL \
    2> >(tee "$logfile")

# Perform realignment around INDELs 
echo -e "\n${RED}6/${STEPS}: GATK:IndelRealigner${NC}"
logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_IndelRealigner.err
    	
java -Xmx4g -jar ${GATK}/GenomeAnalysisTK.jar \
  -T IndelRealigner \
  -R ${REF} \
  -I ${OUT}/${SAMPLE}/${SAMPLE}_dup.bam \
  -targetIntervals ${OUT}/${SAMPLE}/logs/${SAMPLE}_dup_intervals.list \
  -o ${OUT}/${SAMPLE}/${SAMPLE}_realn.bam \
  -l FATAL \
  2> >(tee "$logfile")	


# Delete surplus files
rm ${OUT}/${SAMPLE}/${SAMPLE}_dup.bam
rm ${OUT}/${SAMPLE}/${SAMPLE}_dup.bai


# ==============================================================================
# BQSR based on boostrapping UnifiedGenotyper and BaseRecalibrator
# ==============================================================================

echo -e "\n${RED}7/${STEPS}: GATK:Base Quality Score Recalibration (BQSR)${NC}"

# Call SNPs > BQSR > Call SNPs > BQSR > Call SNPs > BQSR > Call SNPs
cp ${OUT}/${SAMPLE}/${SAMPLE}_realn.bam ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam

# Boostrap a number of times
for (( n = 1; n <= ${BOOTS}; n++ ))
do
  
  echo -e "Processing iteration ${n}"

  # Sort and Index bootstram bam
  #samtools sort ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap
  samtools sort -o ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam
  samtools index ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam 
  
  echo -e "index OK"

 
    # GATK UnifiedGenotyper on realn BAM file, qual 20
    logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_UnifiedGenotyper.err
    java -Xmx4g -jar ${GATK}/GenomeAnalysisTK.jar \
      -T UnifiedGenotyper \
      -R ${REF} \
      -I ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
      --genotyping_mode DISCOVERY \
      --sample_ploidy 1 \
      --min_base_quality_score 20 \
      -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_${n}.vcf \
      -l FATAL \
      2> >(tee "$logfile")

    # Filter SNPs for high quality
    logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_SNPfilter.err
    java -Xmx4g -jar ${GATK}/GenomeAnalysisTK.jar \
      -T VariantFiltration \
      -R ${REF} \
      -V ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_${n}.vcf \
      --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
      --filterName "HQ_fail" \
      --genotypeFilterExpression "GQ < 90.0" \
      --genotypeFilterName "GQ_fail" \
      -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_filtered_${n}.vcf \
      -l FATAL \
      2> >(tee "$logfile")

    # Extract SNPs from VCF that passed the filters
    logfile=${OUT}/${SAMPLE}/${SAMPLE}_SelectVariants.err
    java -Xmx4g -jar ${GATK}/GenomeAnalysisTK.jar \
      -T SelectVariants \
      -R ${REF} \
      --variant ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_filtered_${n}.vcf \
      -o ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_pass_${n}.vcf \
      -l FATAL \
      --excludeFiltered

    # Analyze patterns of covariation in the sequence dataset
    logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_BaseRecalibrator_pass1.err
    java -Xmx4g -jar ${GATK}/GenomeAnalysisTK.jar \
      -T BaseRecalibrator \
      -R ${REF} \
      -I ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
      -knownSites:VCF ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_pass_${n}.vcf \
      -o ${OUT}/${SAMPLE}/${SAMPLE}_BQSR_${n}.table \
      -l FATAL \
      2> >(tee "$logfile")

    # Do a second pass to analyze covariation post-recalibration
    logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_BaseRecalibrator_pass2.err
    java -Xmx4g -jar ${GATK}/GenomeAnalysisTK.jar \
      -T BaseRecalibrator \
      -R ${REF} \
      -I ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
      -knownSites:VCF ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn_pass_${n}.vcf \
      -BQSR ${OUT}/${SAMPLE}/${SAMPLE}_BQSR_${n}.table \
      -o ${OUT}/${SAMPLE}/${SAMPLE}_post_BQSR_${n}.table \
      -l FATAL \
      2> >(tee "$logfile")

    # Apply the recalibration
    logfile=${OUT}/${SAMPLE}/logs/${SAMPLE}_ApplyRecalibration.err
    java -Xmx4g -jar ${GATK}/GenomeAnalysisTK.jar \
      -T PrintReads \
      -R ${REF} \
      -I ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam \
      -BQSR ${OUT}/${SAMPLE}/${SAMPLE}_post_BQSR_${n}.table \
      -o ${OUT}/${SAMPLE}/${SAMPLE}_tmp.bam \
      --emit_original_quals \
     -l FATAL \
      2> >(tee "$logfile")
        
  cp ${OUT}/${SAMPLE}/${SAMPLE}_tmp.bam ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam

done

# Re-index bam file resulting from the above BQSR bootstrapping
samtools index ${OUT}/${SAMPLE}/${SAMPLE}_bootstrap.bam

# Remove surplus bootstrap files and realn.bam
rm ${OUT}/${SAMPLE}/${SAMPLE}_tmp.ba*
rm ${OUT}/${SAMPLE}/${SAMPLE}*BQSR*
rm ${OUT}/${SAMPLE}/vcfs/${SAMPLE}_GATK_realn*
rm ${OUT}/${SAMPLE}/${SAMPLE}_realn.b*




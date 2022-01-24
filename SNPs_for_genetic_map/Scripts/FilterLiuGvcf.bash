module load bioinfo/Java8

#Select SNP
java -d64 -Xmx32g -jar /usr/local/bioinfo/src/GATK/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R /home/gencel/vignal/save/Genomes/Abeille/PacificBee/BlackBee_genome.AllAlternative.NoMultipleAllele.fasta \
    -V /work/project/cytogen/Alain/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF.vcf \
    -selectType SNP \
    -o /work/project/cytogen/Alain/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF_raw_snps.vcf

#Filter SNP
java -d64 -Xmx32g -jar /usr/local/bioinfo/src/GATK/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R /home/gencel/vignal/save/Genomes/Abeille/PacificBee/BlackBee_genome.AllAlternative.NoMultipleAllele.fasta \
    -V /work/project/cytogen/Alain/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF_raw_snps.vcf \
    --filterExpression "FS > 60.0 || MQ < 40.0 || SOR > 3.0" \
    --filterName "snp_filter" \
    -o /work/project/cytogen/Alain/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF_filtered_snps.vcf

#Extract SNPs from VCF that passed the filters
java -d64 -Xmx32g -jar /usr/local/bioinfo/src/GATK/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R /home/gencel/vignal/save/Genomes/Abeille/PacificBee/BlackBee_genome.AllAlternative.NoMultipleAllele.fasta \
    -V /work/project/cytogen/Alain/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF_filtered_snps.vcf \
    --excludeFiltered \
    -o /work/project/cytogen/Alain/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF_snps_OK.vcf
    
#Clean up and compress
bgzip -f /work/project/cytogen/Alain/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF_snps_OK.vcf
rm /work/project/cytogen/Alain/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF_raw_snps.vcf*
rm /work/project/cytogen/Alain/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF_filtered_snps.vcf*

#Filter out 


#Select InDels
java -d64 -Xmx32g -jar /usr/local/bioinfo/src/GATK/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R /home/gencel/vignal/save/Genomes/Abeille/PacificBee/BlackBee_genome.AllAlternative.NoMultipleAllele.fasta \
    -V /work/project/cytogen/Alain/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF.vcf \
    -selectType INDEL \
    -o /work/project/cytogen/Alain/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF_raw_indels.vcf
    
#Filter InDels
java -d64 -Xmx32g -jar /usr/local/bioinfo/src/GATK/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R /home/gencel/vignal/save/Genomes/Abeille/PacificBee/BlackBee_genome.AllAlternative.NoMultipleAllele.fasta \
    -V /work/project/cytogen/Alain/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF_raw_indels.vcf \
    --filterExpression "FS > 200.0 || SOR > 10.0" \
    --filterName "indel_filter" \
    -o /work/project/cytogen/Alain/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF_filtered_indels.vcf

#Extract InDels from VCF that passed the filters
java -d64 -Xmx32g -jar /usr/local/bioinfo/src/GATK/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R /home/gencel/vignal/save/Genomes/Abeille/PacificBee/BlackBee_genome.AllAlternative.NoMultipleAllele.fasta \
    -V /work/project/cytogen/Alain/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF_filtered_indels.vcf \
    --excludeFiltered \
    -o /work/project/cytogen/Alain/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF_indels_OK.vcf

##Clean up and compress
bgzip -f /work/project/cytogen/Alain/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF_indels_OK.vcf
rm /work/project/cytogen/Alain/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF_raw_indels.vcf*
rm /work/project/cytogen/Alain/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF_filtered_indels.vcf*


module load bioinfo/Java8

java -d64 -Xmx32g -jar /usr/local/bioinfo/src/GATK/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
	-T GenotypeGVCFs \
	-R /home/gencel/vignal/save/Genomes/Abeille/PacificBee/BlackBee_genome.AllAlternative.NoMultipleAllele.fasta \
	-V /work/project/cytogen/Alain/LiuRecombination/AlignmentsHaplotypeCaller/GenotypeGVCF.list \
	-o /work/project/cytogen/Alain/LiuRecombination/AlignmentsHaplotypeCaller/LiuAllGenotypesGVCF.vcf

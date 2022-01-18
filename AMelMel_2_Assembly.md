# Assembly

## Assembly of reads into contigs

* Raw PacBio reads from the single individual were assembled with Canu 1.3 (Koren *et al.* 2017) using standard parameters
* a first polishing of the assembly was done with quiver (version SMRT_Link v4.0.0) (https://github.com/PacificBiosciences/GenomicConsensus) using standard parameters.
* Error correction was then performed with Illumina reads from the same individual sequenced with an Illumina NovaSeq6000 instrument, producing over 28 000 000 reads (estimated raw sequencing depth = 33.7 X). 
* Contigs were then assigned to chromosomes by alignment to the Amel4.5 reference genome using LAST (Frith and Kawaguchi 2015).

## From contigs to chromosomes with a genetic map

Data from Lui et al (2015) was used to order and orient contigs with a genetic map.

### Download fastq from SRA




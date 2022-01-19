#!/usr/bin/env python3
# -*- coding: Utf-8 -*-

#Generate_Agp_Fasta.py

"""
Alain Vignal 18/01/2022
"""

"""
Edit the path to the fasta file containing the contigs.
Edit the name of the file with the contig order.
"""
################ EDITS #######################
path_to_fasta = '/home/gencel/vignal/save/Genomes/Abeille/PacificBee/BlackBee_genome.AllAlternative.NoMultipleAllele.fasta' #fasta file containing the contigs
file_contig_order_orientation = '/home/gencel/vignal/save/Genomes/Abeille/PacificBee/assemblingV2/combined_infoV2MT.txt' #file with the contig order: chromosome(int), order(int), contigName, orientation (+ or -)
nb_ns = 100 #number of N to insert between contigs
################ END EDITS ###################

import argparse
import csv
import re

from collections import defaultdict

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.SeqRecord import SeqRecord

fileOut = open('AMelMelPacBioChromosomes.agp', 'w')
fileOut2 = open('AMelMelPacBioChrUnknown.agp', 'w')
fileOut3 = open('AMelMelPacBioUnlocalized.agp', 'w')

'''
Build dictionary of contig lenths and list of all fasta contig sequence records
this will avoid searching the whole assembly text file, just the list
'''
contigLengths = dict()
genomeRecords = list()
inputFasta = open(path_to_fasta, 'r')
for genomeRecord in SeqIO.parse(inputFasta, "fasta"):
	contigLengths[genomeRecord.id] = str(len(genomeRecord.seq))
	genomeRecords.append(SeqRecord(seq=genomeRecord.seq, id=genomeRecord.id))

"""
Build dictionary of chromosomes, contigs, order and orientation, from imput txt file.
"""
chrBase = defaultdict(dict)
unplacedBase = defaultdict(dict)
unlocalizedBase = defaultdict(dict)

with open(file_contig_order_orientation, newline='') as csvFile:
	data = csv.reader(csvFile, delimiter = "\t")
	for row in data:
		if re.match("\A\d",row[35]):		#Contigs assigned to a chromosome: name starts with digit
			chrBase[int(row[35])][int(row[36])] = (row[1], row[38]) #Assign contig name and orientation to chromosome number and order
		elif re.match("\AUn\Z", row[35]):	#ChrUnknown: name = Un
			unplacedBase[int(row[0])] = row[1]
		elif re.match("\AUn\d",row[35]):	#Unlocalized: name = UnXX
			unlocalizedBase[row[35]][int(row[0])] = (row[1])
#		elif re.match("\AMT\Z",row[35]):	#Mitochondria: name = MT
#			chrBase[int(row[35])][int(row[36])] = (row[1], row[38]) #Assign contig name and orientation to chromosome number and order

"""
Print assembled chromosomes agp files
"""
chrList = sorted(list(chrBase.keys()))
for chromosome in chrList:
	if chromosome == 17:
			chrom = "MT"
	else:
		chrom = "".join(("LG", str(chromosome)))
	orderList = sorted(list(chrBase[chromosome].keys()))
	count = 1
	start = 0
	end = 0
	for order in orderList:
		contigName = chrBase[chromosome][order][0]
		start = end + 1
		end = end + int(contigLengths[contigName])
		contigOrientation = chrBase[chromosome][order][1]
		contigLine = (chrom,str(start),str(end),str(count),"W",contigName,"1",contigLengths[contigName],contigOrientation)
		print ("\t".join(contigLine), file = fileOut)
		count = count + 1
		start = end + 1
		end = end + nb_ns
		if (count / 2) < len(orderList):
			spacingNs = (chrom,str(start),str(end),str(count),"U",str(nb_ns),"scaffold","yes","map")
			print ("\t".join(spacingNs), file = fileOut)
			count = count + 1

"""
Print Unlocalized contigs agp files (known chromosome, but unknown location)
"""
chrList = sorted(list(unlocalizedBase.keys()))
for chromosome in chrList:
	chrom = "".join(("LG_", chromosome))
	orderList = sorted(list(unlocalizedBase[chromosome].keys()))
	count = 1
	start = 0
	end = 0
	for order in orderList:
		contigName = unlocalizedBase[chromosome][order]
		start = end + 1
		end = end + int(contigLengths[contigName])
		contigLine = (chrom,str(start),str(end),str(count),"W",contigName,"1",contigLengths[contigName],"+")
		print ("\t".join(contigLine), file = fileOut3)
		count = count + 1
		start = end + 1
		end = end + nb_ns
		if (count / 2) < len(orderList):
			spacingNs = (chrom,str(start),str(end),str(count),"U",str(nb_ns),"contig","no","NA")
			print ("\t".join(spacingNs), file = fileOut3)
			count = count + 1

"""
Print ChrUnknown agp files (chromosome assignation unknown)
"""
contigList = sorted(list(unplacedBase.keys()))
for contigNb in contigList:
	contigName = "unplaced_" + unplacedBase[contigNb]
	contigLine = (contigName, "1", contigLengths[unplacedBase[contigNb]], "1", "W", unplacedBase[contigNb], "1", contigLengths[unplacedBase[contigNb]], "+")
	print ("\t".join(contigLine), file = fileOut2)

"""
Print sequences files in fasta format
"""
#Assembled chromosomes
chrList = sorted(list(chrBase.keys()))
toInsertNs = 'N' * nb_ns
for chromosome in chrList:
	chrom = "".join(("chr", str(chromosome)))
	orderList = sorted(list(chrBase[chromosome].keys()))
	count = 1
	chrSequence = ""
	for order in orderList:
		contigName = chrBase[chromosome][order][0]
		contigOrientation = chrBase[chromosome][order][1]
		for contig in genomeRecords:
			if contig.id == chrBase[chromosome][order][0] and contigOrientation == '+':
				sequence = str(contig.seq)
			elif contig.id == chrBase[chromosome][order][0] and contigOrientation == '-':
				sequence = str(contig.seq.reverse_complement())
			else:
				continue
		chrSequence = chrSequence + sequence
		count = count + 1
		if (count / 2) < len(orderList):
			chrSequence = chrSequence + toInsertNs
			count = count +1
	chrSequenceComplete = Seq(chrSequence)
	chrSequenceComplete = SeqRecord(chrSequenceComplete, id=chrom, description="")
	fileName = 'AMelMelPacBio_' + chrom + '.fa'
	with open(fileName, 'w') as fout:
		output_fasta = FastaWriter(fout, wrap=60)
		output_fasta.write_header()
		output_fasta.write_record(chrSequenceComplete)
		output_fasta.write_footer()
#ChrUnknown
contigList = sorted(list(unplacedBase.keys()))
with open('chrUnknown.fa', 'w') as fout:
	output_fasta = FastaWriter(fout, wrap=60)
	output_fasta.write_header()
	for contigNb in contigList:
		for contig in genomeRecords:
			if contig.id == unplacedBase[contigNb]:
				contig.description = ""
				output_fasta.write_record(contig)
	output_fasta.write_footer()

#Unlocalized
chrList = sorted(list(unlocalizedBase.keys()))
for chromosome in chrList:
	chrom = "".join(("chr", chromosome))
	orderList = sorted(list(unlocalizedBase[chromosome].keys()))
	fileName = 'AMelMelPacBio_' + chrom + '.fa'
	with open(fileName, 'w') as fout:
		output_fasta = FastaWriter(fout, wrap=60)
		output_fasta.write_header()
		for order in orderList:
			for contig in genomeRecords:
				if contig.id == unlocalizedBase[chromosome][order]:
					contig.description = ""
					output_fasta.write_record(contig)
		output_fasta.write_footer()

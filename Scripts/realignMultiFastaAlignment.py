#!/usr/bin/env python3
# -*- coding: Utf-8 -*-

#For changing the start position in a multiple alignment of circular sequences of tandem repeats consensus
#Input : multi fasta file after alignment
#Output: multi fasta file with all sequences starting at a new position

import argparse
import csv

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.SeqRecord import SeqRecord

def main (input_file, output_file, start):
	fileOut = open(output_file, 'w')
	for record in SeqIO.parse(input_file, "fasta"):
#		id = record.id.split('|')[1] if '|' in record.id else record.id
		sequence = str(record.seq)
		out = Seq(sequence)
		out = SeqRecord(out, id=record.id, description="")
		writer = FastaWriter(fileOut, wrap=60)
		writer.write_header()
		writer.write_record(out[int(start):len(record.seq)] + out[0:int(start)])
		writer.write_footer

def parseArguments():
	'''Parse input arguments.
	IN:		None
	OUT:	parsed arguments
	'''
	parser = argparse.ArgumentParser(description='Input : multi fasta file after alignment with mafft\n \
									Output: multi fasta file with all sequences starting at equivalent positions)')

	parser.add_argument('-i', '--input_file', nargs='?', required=True,
						help='input multi-fasta file with aligned sequences to reorganise.')
	parser.add_argument('-o', '--output_file', nargs='?', required=False, default='fout', #couldnt find how to output to screen by default
						help='output file (default file named fout).')
	parser.add_argument('-s', '--start', nargs='?', required=True,
						help='start of sequences')

	return parser.parse_args()


if __name__ == '__main__':
	args = parseArguments()
	main(args.input_file, args.output_file, args.start)

# end of file

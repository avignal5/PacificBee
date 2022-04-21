#!/usr/bin/env python3
# -*- coding: Utf-8 -*-

#Extracts one sequence from a multifasta file.
#Either the complete sequence or start-end coordinates can be provided

import argparse

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.SeqRecord import SeqRecord

def main(input_file, output_file, name, start, end):
	with open(output_file, "w") as fileOut:
		for record in SeqIO.parse(input_file, "fasta"):
			id = record.id.split('|')[1] if '|' in record.id else record.id
			if id in name:
				writer = FastaWriter(fileOut, wrap=60)
				writer.write_header()
				if end == 1:
					writer.write_record(record[start-1:len(record.seq)]) #Writes a single record. For several, must use write_records
				elif end != 1:
					sequence = str(record.seq)
					out = Seq(sequence)
					out = SeqRecord(out, id=name + "_" + str(start) + "_" + str(end), description="Sequence of " + name + ", from " + str(start) + " to " + str(end))
					writer.write_record(out[start-1:end])
				writer.write_footer

def parseArguments():
	'''Parse input arguments.
	IN:		None
	OUT:	parsed arguments
	'''
	parser = argparse.ArgumentParser(description='Extracting fasta sequences')
	
	parser.add_argument('-i', '--input_file', nargs='?', required=True,
						help='input fasta file.')
	parser.add_argument('-o', '--output_file', nargs='?', required=False, default='out.fa', #Writes in a file named out.fa Couldn't manage to send to screen by default
						help='output file. Default = out.fa')
	parser.add_argument('-n', '--name', nargs='?', required=True,
						help='sequence name')
	parser.add_argument('-s', '--start', nargs='?', required=False, type=int, default=1,
						help='start position. Default = 1')
	parser.add_argument('-e', '--end', nargs='?', required=False, type=int,  default=1,
						help='end position. Default = all the sequence')
						
	return parser.parse_args()


if __name__ == '__main__':
	args = parseArguments()	
	main(args.input_file, args.output_file, args.name, args.start, args.end)
	
# end of file

#!/usr/bin/env python3
# -*- coding: Utf-8 -*-

#PlotAlignments.py

import re
import pylab as pl
from matplotlib import collections  as mc
from matplotlib.colors import ListedColormap
import argparse
import csv

"""
mafToTxt:
Function for converting a maf file to a compact txt file, to be later used for plotting.
Redundant with the use of psl files or tab files produced by maf-convert (see LAST documentstion)
And probably slow (check)
"""
def mafToTxt(input_file,data):
	seq = 0
	fileName = re.split("\.",input_file)
	outTextFile = open(fileName[0] + ".txt", "w")
	header = ("Score", "Mismap", "tName", "tStart", "tEnd", "tAlignLength", "tStrand", "tLength", "qName", "qStart", "qEnd", "qAlignLength", "qStrand", "qLength", "coordinates")
	print("\t".join(header), file=outTextFile)
	for line in data:
		if re.search('^a',line):
			row = re.split("\s+",line)
			score = row[1].split('=')[1]
			mismap = row[2].split('=')[1]
		if (re.search('^s',line)) and (seq == 0):
			row = re.split("\s+",line)
			tName = row[1]
			tStart = int(row[2])
			tAlignLength = int(row[3])
			tEnd = tStart + tAlignLength
			tStrand = row[4]
			tLength = int(row[5])
			seq = 1
		elif (re.search('^s',line)) and (seq == 1):
			row = re.split("\s+",line)
			qName = row[1]
			qAlignLength = int(row[3])
			qStrand = row[4]
			qLength = int(row[5])
			if qStrand == "+":
				qStart = int(row[2])
				qEnd = qStart + qAlignLength
			if qStrand == "-":
				qEnd = qLength - int(row[2])
				qStart = qEnd - qAlignLength
			seq = 0
			lineOut = (str(score),str(mismap),tName,str(tStart),str(tEnd),str(tAlignLength),tStrand,str(tLength),qName,str(qStart),str(qEnd),str(qAlignLength),qStrand,str(qLength))
			print("\t".join(lineOut), file=outTextFile)

"""
Read lines in the input (usually psl) files
"""
def getData(input_file, line, scale_factor):
	if (re.search('.txt$', input_file)):
		if re.search('^[0-9]',line):
			row = re.split("\t",line)
			compData = {'tChr' : row[2],
						'tStart' : int(row[3])/scale_factor,
						'tEnd' : int(row[4])/scale_factor,
						'qChr' : row[8],
						'qStart' : int(row[9])/scale_factor,
						'qEnd' : int(row[10])/scale_factor,
						'strand' : row[12]}
			return compData
	elif (re.search('.psl$', input_file)):
		if re.search('^[0-9]',line):
			row = re.split("\s+",line)
			compData = {'qChr' : row[9],
						'qStart' : int(row[11])/scale_factor,
						'qEnd' : int(row[12])/scale_factor,
						'tChr' : row[13],
						'tStart' : int(row[15])/scale_factor,
						'tEnd' : int(row[16])/scale_factor,
						'strand' : row[8]}
			return compData

"""
Coordinates of alignment blocks, according to:
- plotting option => which genome is on which axis
- if alignment is on the forward or reverse strand
"""
def buildCoordinates(qStart, tStart, qEnd, tEnd, strand, swap_xy):
	if swap_xy == 'False':
		if strand == '+':
			forCoordinates = [(qStart,tStart),(qEnd,tEnd)]
		if strand == '-':
			forCoordinates = [(qEnd,tStart),(qStart,tEnd)]
		return forCoordinates
	elif swap_xy == 'True':
		if strand == '+':
			forCoordinates = [(tStart,qStart),(tEnd,qEnd)]
		if strand == '-':
			forCoordinates = [(tStart,qEnd),(tEnd,qStart)]
		return forCoordinates

"""
To allow forward and reverse matches to be colored differently
"""
def colorStrand(strand):
	if strand == '+':
		return('b') #b => blue
	else:
		return('r') #r => red
"""
To allow for different x and y axes legends than the chromosome name from the psl file.
"""
def axLegends(query, target, swap_xy):
	if swap_xy == 'False':
		pl.xlabel(query)
		pl.ylabel(target)
	elif swap_xy == 'True':
		pl.xlabel(target)
		pl.ylabel(query)

def main (input_file, output_fig_file, target_chromosome, target_start, target_end, query_chromosome, query_start, query_end, \
			vertical_lines, horizontal_lines, swap_xy, line_width, scale_factor, x_axis_legend, y_axis_legend):

	data = open(input_file, "r")
	coordinates = []
	colors = []
#	cmap = ListedColormap(['r','b','g'])

	if (re.search('.maf$', input_file)):	#if a maf file given, will convert to txt file
		mafToTxt(input_file,data)
	#else: #Tenter comme ça, vu que les différents cas sont traités dans def getData(input_file, line):
	else:
		if target_chromosome != "allchromosomes" and query_chromosome != "allchromosomes":
			if target_end == 1 and query_end == 1:
				for line in data:
					compData = getData(input_file, line, scale_factor)
					if compData:
						if compData['tChr'] == target_chromosome and compData['qChr']  == query_chromosome:
							forCoordinates = buildCoordinates(compData['qStart'],compData['tStart'],compData['qEnd'],compData['tEnd'],compData['strand'],swap_xy)
							coordinates.append(forCoordinates)
							colors.append(colorStrand(compData['strand']))
			elif target_end != 1 and query_end == 1:
				for line in data:
					compData = getData(input_file, line, scale_factor)
					if compData:
						if compData['tChr'] == target_chromosome and compData['qChr']  == query_chromosome:
							if (compData['tStart'] <= target_start and compData['tEnd'] >= target_start) \
							or (compData['tStart'] <= target_end and compData['tEnd'] >= target_end) \
							or (compData['tStart'] >= target_start and compData['tEnd'] <= target_end):
								forCoordinates = buildCoordinates(compData['qStart'],compData['tStart'],compData['qEnd'],compData['tEnd'],compData['strand'],swap_xy)
								coordinates.append(forCoordinates)
								colors.append(colorStrand(compData['strand']))
			elif target_end == 1 and query_end != 1:
				for line in data:
					compData = getData(input_file, line, scale_factor)
					if compData:
						if compData['tChr'] == target_chromosome and compData['qChr']  == query_chromosome:
							if (compData['qStart'] <= query_start and compData['qEnd'] >= query_start) \
							or (compData['qStart'] <= query_end and compData['qEnd'] >= query_end) \
							or (compData['qStart'] >= query_start and compData['qEnd'] <= query_end):
								forCoordinates = buildCoordinates(compData['qStart'],compData['tStart'],compData['qEnd'],compData['tEnd'],compData['strand'],swap_xy)
								coordinates.append(forCoordinates)
								colors.append(colorStrand(compData['strand']))
			elif target_end != 1 and query_end != 1:
				for line in data:
					compData = getData(input_file, line, scale_factor)
					if compData:
						if compData['tChr'] == target_chromosome and compData['qChr']  == query_chromosome:
							if ((compData['qStart'] <= query_start and compData['qEnd'] >= query_start) \
							or (compData['qStart'] <= query_end and compData['qEnd'] >= query_end) \
							or (compData['qStart'] >= query_start and compData['qEnd'] <= query_end)) \
							and ((compData['tStart'] <= target_start and compData['tEnd'] >= target_start) \
							or (compData['tStart'] <= target_end and compData['tEnd'] >= target_end) \
							or (compData['tStart'] >= target_start and compData['tEnd'] <= target_end)):
								forCoordinates = buildCoordinates(compData['qStart'],compData['tStart'],compData['qEnd'],compData['tEnd'],compData['strand'],swap_xy)
								coordinates.append(forCoordinates)
								colors.append(colorStrand(compData['strand']))

		outFigFile = open(output_fig_file, "w")
		fig, ax = pl.subplots()

#Plot vertical and horizontal lines
		vertical_color = 'black'
		horizontal_color = 'black'
#Plot vertical lines at breakpoints (or horizontal if swap)
		if not vertical_lines == "None":
			with open(vertical_lines, "r") as csvFile:
				data=csv.reader(csvFile, delimiter = "\t")
				for row in data:
					if row[0] == query_chromosome:
						if swap_xy == 'False':
							ax.axvline(x=float(row[1])/scale_factor,color=vertical_color,alpha=0.5, linewidth=0.5, linestyle='--')
						elif swap_xy == 'True':
							ax.axhline(y=float(row[1])/scale_factor,color=vertical_color,alpha=0.5, linewidth=0.5, linestyle='--')
#Plot horizontal lines at breakpoints (or vertical if swap)
		if not horizontal_lines == "None":
			with open(horizontal_lines, "r") as csvFile:
				data=csv.reader(csvFile, delimiter = "\t")
				for row in data:
					if row[0] == target_chromosome:
						if swap_xy == 'False':
							ax.axhline(y=float(row[1])/scale_factor,color=horizontal_color,alpha=0.5, linewidth=0.5, linestyle='--')
						elif swap_xy == 'True':
							ax.axvline(x=float(row[1])/scale_factor,color=horizontal_color,alpha=0.5, linewidth=0.5, linestyle='--')
		lc = mc.LineCollection(coordinates, linewidths=line_width, colors=colors)
		ax.add_collection(lc)
		ax.set_aspect('equal','box')
		#ax.autoscale()
		ax.margins(0.0)
		if target_end != 1 and query_end == 1:
			ax.set_ylim(target_start, target_end)
		elif target_end == 1 and query_end != 1:
			ax.set_xlim(query_start, query_end)
		elif target_end != 1 and query_end != 1:
			if swap_xy == 'False':
				ax.set_xlim(query_start, query_end)
				ax.set_ylim(target_start, target_end)
			elif swap_xy == 'True':
				ax.set_xlim(target_start, target_end)
				ax.set_ylim(query_start, query_end)

		if x_axis_legend == "None" and y_axis_legend == "None":
			axLegends(query_chromosome,target_chromosome, swap_xy)
		elif x_axis_legend != "None" and y_axis_legend == "None":
			axLegends(x_axis_legend,target_chromosome, swap_xy)
		elif x_axis_legend == "None" and y_axis_legend != "None":
			axLegends(query_chromosome,x_axis_legend, swap_xy)
		elif x_axis_legend != "None" and y_axis_legend != "None":
			axLegends(x_axis_legend,y_axis_legend, swap_xy)

#Create figure as png or pdf

		if (re.search('.png',output_fig_file)):
			pl.savefig(output_fig_file, format="png")
		if (re.search('.pdf',output_fig_file)):
			pl.savefig(output_fig_file, format="pdf")

def parseArguments():
	'''Parse input arguments.
	IN:		None
	OUT:	parsed arguments
	'''
	parser = argparse.ArgumentParser(description='Plot alignments')

	parser.add_argument('-i', '--input_file', nargs='?', required=True,
						help='input psl, maf or txt file.')
	parser.add_argument('-o', '--output_fig_file', nargs='?', required=False, default='alignmentOut.png', #Writes in a file named alignmentOut.png, or .pdf. Couldn't manage to send to screen by default
						help='output file. Default = alignmentOut.png')						             #May work if default is a switch to a print function with no file specified.
	parser.add_argument('-c', '--target_chromosome', nargs='?', required=False, default="allchromosomes",
						help='target chromosome name')
	parser.add_argument('-s', '--target_start', nargs='?', required=False, type=int, default=0,
						help='target chromosome start position. Default = 0')
	parser.add_argument('-e', '--target_end', nargs='?', required=False, type=int,  default=1,
						help='target end position. Default = all the sequence')
	parser.add_argument('-q', '--query_chromosome', nargs='?', required=False, default="allchromosomes",
						help='query chromosome name. Default = all the sequence')
	parser.add_argument('-t', '--query_start', nargs='?', required=False, type=int,  default=0, #was default=1 in older versions.
						help='query start position. Default = all the sequence')
	parser.add_argument('-f', '--query_end', nargs='?', required=False, type=int,  default=1,
						help='query end position. Default = all the sequence')
	parser.add_argument('-v', '--vertical_lines', nargs='?', required=False, default="None",
						help='File for vertical lines. Default = None')
	parser.add_argument('-z', '--horizontal_lines', nargs='?', required=False, default="None",
						help='File for horizontal lines. Default = None')
	parser.add_argument('-w', '--swap_xy', nargs='?', required=False, default="False",
						help='Swap X and Y. True or False. Default = False')
	parser.add_argument('-l', '--line_width', nargs='?', required=False, type=float,  default=1.0,
						help='Line thickness. Default = 1.0')
	parser.add_argument('-r', '--scale_factor', nargs='?', required=False, type=int,  default=1000000,
						help='Scale of DNA sequence. Default = 1000000, to output in Mb')
	parser.add_argument('-x', '--x_axis_legend', nargs='?', required=False, type=str,  default="None",
						help='Chromosome name legend. Default = the names from the fasta file')
	parser.add_argument('-y', '--y_axis_legend', nargs='?', required=False, type=str,  default="None",
						help='Chromosome name legend. Default = the names from the fasta file. Default = 1.0')

	return parser.parse_args()

if __name__ == '__main__':
	args = parseArguments()
	main(args.input_file, args.output_fig_file, args.target_chromosome, args.target_start, args.target_end, args.query_chromosome, args.query_start, args.query_end, \
			args.vertical_lines, args.horizontal_lines, args.swap_xy, args.line_width, args.scale_factor, args.x_axis_legend, args.y_axis_legend)

# end of file

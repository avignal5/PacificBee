################################################
### script to edit genotype file by changing marker id for contig_position and individual id by sample names
################################################
#!/usr/bin/env python3

import pandas as pd
import numpy as np
import csv
import re
import sys

def read_in(data):
	'''
	Function to read in data file 
	'''
	geno=open(data,'r')
	geno=geno.readlines()
	geno=[x.strip().split() for x in geno]
	genotype=pd.DataFrame(geno[1:],columns=geno[0])
	return genotype
	
def rename_markers(data,pattern):
	'''
	Function to rename markers as contig_position 
	'''
	datac=data['contigID'].str.lstrip(pattern).str.lstrip('0')
	m_id=datac+'_'+data['position']
	for i in range(0,len(m_id)):
		if m_id[i][0]=='_':
			m_id[i]='0'+m_id[i]
	data.drop(data.columns[[0,1]],axis=1,inplace=True)
	data.insert(0,'marker_id',m_id)
	return data
	
def column_rename_geno(geno,idinfo,idname):
	'''
	Function to rename the columns of the genotype file by the individual ids from the individual information file
	'''
	head=list(geno)
	id_new=[]
	for i in range(0,len(idinfo)):
		k=idinfo[idinfo==head[i+1]].index.tolist()[0]
		id_new.append(idname[k])
	col_new='marker_id'+' '+' '.join(id_new)
	geno.columns=[col_new.split()]
	return geno

###########################################################################
id=read_in(sys.argv[1])
geno_map_i=read_in(sys.argv[3])
geno_map_f=rename_markers(geno_map_i,sys.argv[6])
geno_map_f.to_csv(sys.argv[5],sep=' ',index=False)
print('genotype_map written')
geno_i=read_in(sys.argv[2])
geno_f1=rename_markers(geno_i,sys.argv[6])
geno_f=column_rename_geno(geno_f1,id.Run_s,id.Sample_Name_s)
geno_f.to_csv(sys.argv[4],sep=' ',index=False)
print('genotype written')

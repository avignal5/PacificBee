################################################
### script to create bins of unique vectors along contigs
################################################
#!/usr/bin/env python3

import pandas as pd
import numpy as np
import csv
import sys
import os

def read_in(data):
	'''
	Function to read in the genotype data 
	'''
	geno=open(data,'r')
	geno=geno.readlines()
	geno=[x.strip().split() for x in geno]
	genotype=pd.DataFrame(geno[1:],columns=geno[0])
	return genotype

def nttonb(data):	
	'''
	Function to convert genotype from nucleotide to 0/1
	'''
	geno=[]
	for i,row in data.iterrows():
		seen=set()
		uniq=[y for y in row if y not in seen and not seen.add(y)][1:]
		v=row.tolist()
		v_sub=''.join(v[1:])
		v_order=v_sub.replace(min(uniq),'0').replace(max(uniq),'1')
		v_reverse=''.join('1' if x=='0' else '0' for x in v_order)
		cont=v[0].split('_')[0]
		marker=v[0].split('_')[1]
		geno.append([v[0],cont,marker,min(uniq),max(uniq),v_order,v_reverse])	
	g=pd.DataFrame(geno)
	g.columns=['marker_id','contig','marker','0','1','vector','vector_']
	return g
	
def bin_creation(data):
	'''
	Function to identify unique genotype vectors and create bins of markers with this genotype
	'''
	v=data.vector.iloc[0]
	v_=data.vector_.iloc[0]
	if(v.startswith('0')==True):
		bin0=[v]
		bin1=[v_]
	else:
		bin0=[v_]
		bin1=[v]
	listm={}
	listm[0]=[]
	i=0
	for m,row in data.iterrows():
		if (row.vector==bin0[i] or row.vector==bin1[i]):
			listm[i].append(row.marker)
		else:
			i=i+1
			if(row.vector.startswith('0')==True):
				bin0.append(row.vector)
				bin1.append(row.vector_)
			else:
				bin0.append(row.vector_)
				bin1.append(row.vector)
			listm[i]=[row.marker]
	length_list=list({key: len(value) for key, value in listm.items()}.values())
	listm_dat=[]
	physical_length=[]
	for l in range(0,len(listm)):
		listm_dat.append(','.join(str(x) for x in listm[l]))
		physical_length.append(max(listm[l])-min(listm[l]))
	physical_length=[1 if x==0 else x for x in physical_length]
	bin_dat={'bin':range(1,len(bin0)+1),'vector0':bin0,'vector1':bin1,'nb_markers':length_list,'group':[str(data.contig.iloc[0])]*len(bin0),'physical_length':physical_length,'list_markers':listm_dat}
	bin_data=pd.DataFrame(bin_dat)
	bin_data=bin_data[['bin','vector0','vector1','nb_markers','group','physical_length','list_markers']]
	return bin_data	
	
##############################################################################################################
nb_colony=int(sys.argv[1])
pattern1=sys.argv[2]
pattern2=sys.argv[3]

for j in range(0,nb_colony):
	geno_name=('%s_colony%d.txt'%(pattern1,(j+1)))
	print('Colony %d'%(j+1))
	geno_c=read_in(geno_name)
	contig_c=set(geno_c.iloc[:,0].str.split('_').str[0])
	print('Number of markers in colony%d:'%(j+1),len(geno_c))
	print('Number of drones in colony%d:'%(j+1),len(geno_c.iloc[1])-1)
	print('number contigs in colony after controls:',len(contig_c))
	geno_nb=nttonb(geno_c)
	geno_nb['contig']=geno_nb['contig'].str.replace('Un','2')
	geno_nb.contig=geno_nb['contig'].str.extract('(\d+)').astype(int)
	geno_nb.marker=geno_nb['marker'].astype(int)
	geno_nb=geno_nb.sort_values(by=['contig','marker'])
	seen=set()
	uniq_cont=[y for y in geno_nb.contig if y not in seen and not seen.add(y)]
	for c in uniq_cont:
		geno_contig=geno_nb[geno_nb.contig==c]
		geno_contig[['contig','marker']] = geno_contig[['contig','marker']].apply(pd.to_numeric)
		if [int(x) for x in geno_contig.marker]!=sorted([int(x) for x in geno_contig.marker]):
			geno_contig.sort_values('marker')
		bins=bin_creation(geno_contig)
		print('Number of unique markers in colony%d contig%d:'%((j+1),int(c)),len(bins))
		bins.to_csv('%s%d/bins_%d.txt' %(pattern2,(j+1),int(c)))
		print('output written')




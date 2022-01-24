################################################
### script to create bins of unique vectors along contigs when 3 colonies combined in 1
################################################
#!/usr/bin/env python3


import pandas as pd
import numpy as np
import csv
import sys
import os
import itertools

##############################################################################################################

geno=pd.read_csv(sys.argv[1],sep=',',header=None)
geno.columns=['chr','contig','marker','vector']
contig_c=set(geno.contig)
print('Number of markers:',len(geno))
print('Number of drones:',len(geno.vector[1]))
print('number contigs in colony after controls:',len(contig_c))
seen=set()
uniq_cont=[y for y in geno.contig if y not in seen and not seen.add(y)]
for c in uniq_cont:
	geno_contig=geno[geno.contig==c]
	geno_contig[['contig','marker']] = geno_contig[['contig','marker']].apply(pd.to_numeric)
	if [int(x) for x in geno_contig.marker]!=sorted([int(x) for x in geno_contig.marker]):
		geno_contig.sort_values('marker')
	geno_contig.reset_index(inplace=True)
	i=0
	j=0
	listm=[geno_contig.marker[i]]
	vect=[geno_contig.vector[i]]
	while(j<(len(geno_contig)-1)):
		vi=geno_contig.vector[i]
		vj=geno_contig.vector[j]
		if(vi==vj):
			j=j+1
			vi=geno_contig.vector[i]
			vj=geno_contig.vector[j]
			listm.append(geno_contig.marker[j])
		else:
			listm.append('break')
			vect.append(geno_contig.vector[j])
			i=j
	listm2=[list(x[1]) for x in itertools.groupby(listm, lambda x: x=='break') if not x[0]] 
	listm_dat=[]
	physical_length=[]
	nb_m=[]
	for l in range(0,len(listm2)):
		listm_dat.append(','.join(str(x) for x in listm2[l]))
		physical_length.append(max(listm2[l])-min(listm2[l]))
		nb_m.append(len(listm2[l]))
	physical_length=[1 if x==0 else x for x in physical_length]
	bin_dat={'chr':geno_contig.chr[0],'contig':c,'bin':range(1,len(vect)+1),'vector':vect,'nb_markers':nb_m,'physical_length':physical_length,'list_markers':listm_dat}
	bin_data=pd.DataFrame(bin_dat)
	print('Number of unique markers in contig%d:'%(int(c)),len(listm2))
	bin_data.to_csv('bins2/bins_%d.txt' %(int(c)))
	print('output written')




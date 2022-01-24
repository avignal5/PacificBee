################################################
### script to extract genotypes from markers in all colonies
################################################
#!/usr/bin/env python3

import pandas as pd
import numpy as np
import csv
import sys
import os
from functools import reduce
import glob

def subset(pattern):
	'''
	Function to extract subset of genotypes present in all colonies 
	'''
	pos={}
	filename=glob.glob('%s*.txt'%pattern)
	for i in range(0,len(filename)):
		col=pd.read_csv(filename[i],sep=' ')
		col.rename(columns={col.columns[0]:'position'},inplace=True)
		pos[i]=col.position
	pos_inter=reduce(set.intersection, map(set, pos.values()))

	for i in range(0,len(filename)):
		col=pd.read_csv(filename[i],sep=' ')
		col.rename(columns={col.columns[0]:'position'},inplace=True)
		col_sub=col[col.position.isin(pos_inter)]
		col_sub.to_csv('%s_all_colony%d.txt' %(pattern,(i+1)),sep=' ',index=False)

###########################################################################################
pattern=sys.argv[1]
subset(pattern)

################################################
### script to perform quality controls on genotype files
################################################
#!/usr/bin/env python3

import pandas as pd
import numpy
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
		
def select_rows(df,search_strings):
    unq,IDs = np.unique(df,return_inverse=True)
    unqIDs = np.searchsorted(unq,search_strings)
    return df[((IDs.reshape(df.shape) == unqIDs[:,None,None]).any(-1)).all(0)]
    
def reseq_check(data,name1,name2):
	'''
	Function to check for sequencing errors by comparing re-sequenced individuals
	'''
	c2=data.filter(regex=name2)
	id_reseq=' '.join(list(c2))
	c1=data[id_reseq.replace(name2,name1).split()]
	m_seq_ok={}
	m_seq_nonok={}
	for j in range(0,len(list(c2))):
		m_s_o=[]
		m_s_no=[]
		for i in range(0,len(c2)):
			if c1.iloc[i,j]==c2.iloc[i,j]:
				if c1.iloc[i,1]!='0/1':
					m_s_o.append(i)
				else:
					m_s_no.append(i) 
			else:
				m_s_no.append(i)
		m_seq_ok[j]=m_s_o
		m_seq_nonok[j]=m_s_no 
	from functools import reduce
	m_seq_ok_tot=reduce(lambda x,y:x&y,map(set,m_seq_ok.values()))
	m_seq_nonok_tot=c2.drop(m_seq_ok_tot).index.tolist()
	return m_seq_ok_tot,m_seq_nonok_tot

def multiallele(data):
	'''
	Function to keep only bi-allelic markers 
	'''
	g_multi_all_colony=[]
	for index, row in geno_f.iterrows():
		if sum(x.isdigit() for x in set(str(set(row[1:]))))>2:
			g_multi_all_colony.append(index)
	return g_multi_all_colony
	
def het_queen_check(data,name1,namequeen):
	'''
	Function to check for queen genotype (kept=heterozygote)
	'''
	geno=data.filter(regex=name1)
	c_queen=geno.filter(regex=namequeen)
	queen_name=name1+namequeen
	m_hetQ=c_queen[c_queen[queen_name].str[0]!=c_queen[queen_name].str[2]].index.tolist() #to keep
	m_homQ=c_queen.drop(m_hetQ).index.tolist() #removed
	v1Q=c_queen[queen_name].str[0]
	v2Q=c_queen[queen_name].str[2]
	return m_hetQ,m_homQ,v1Q,v2Q

def hom_drone_check(data,name1,namequeen,v1Q,v2Q):
	'''
	Function to check for drone genotype (kept=homozygote (haploid))
	'''
	geno=data.filter(regex=name1)
	c_drone=geno.drop(geno.filter(regex=namequeen),axis=1)
	m_error=[]
	m_hom=[]
	m_het=[]
	m_multi=[]
	for i in range (0,len(c_drone)):
		v1q=v1Q[i]
		v2q=v2Q[i]
		a=set(''.join(list(set(c_drone.iloc[i,:]))).replace('/',''))
		if ({v1q,v2q}!=a):
			m_error.append(i)
		if len(a)==1:
			m_hom.append(i)
		elif len(a)>2:
			m_multi.append(i)
		elif len(a)==2:
			if ((list(a)[0]+'/'+list(a)[1] in list(c_drone.iloc[i,:])) or (list(a)[1]+'/'+list(a)[0] in list(c_drone.iloc[i,:]))):
				m_het.append(i)
	m_homk=c_drone.drop(m_error+m_het+m_hom+m_multi).index.tolist()
	return m_homk,m_het,m_hom,m_multi,m_error

###########################################################################
nb_colony=int(sys.argv[2])
geno_f=read_in(sys.argv[1])

contig=geno_f.iloc[:,0].str.split('_').str[0]
print('number markers before controls:',len(geno_f))
print('number contigs before controls:',len(set(contig)))
marker_id=geno_f.iloc[:,0]

m_seq_ok_tot,m_seq_nonok_tot=reseq_check(geno_f,list(sys.argv[3])[0],sys.argv[4])
g_multi_all_colony=multiallele(geno_f)
print('number of multi allelic markers across colonies:',len(g_multi_all_colony))

namequeen='Q'
colony_list=sys.argv[3].split(',')
if len(colony_list)!=nb_colony: 
	print('error number colonies')
keep_remove={}
for i in range(0,nb_colony):
	m_hetQ,m_homQ,v1Q,v2Q=het_queen_check(geno_f,colony_list[i],'Q')
	m_homk,m_het,m_hom,m_multi,m_error=hom_drone_check(geno_f,colony_list[i],'Q',v1Q,v2Q)
	print('number non polymorphic markers in colony%d :'%(i+1),(len(m_hom)))
	print('number homozygous markers in queen%d :'%(i+1),len(m_homQ))
	print('number heterozygous markers in drones from colony%d :'%(i+1),len(m_het))
	print('number multiallelic(>3) markers in colony%d :'%(i+1),(len(m_multi)))
	print('number sequencing errors in drone colony%d :'%(i+1),(len(m_error)))
	to_keep_final=set(m_seq_ok_tot)&set(m_hetQ)&set(m_homk)
	print('number markers kept in colony%d :'%(i+1),len(to_keep_final))
	geno=geno_f.filter(regex=colony_list[i])
	print('number of drones in colony%d :'%(i+1),len(list(geno))-1)
	geno=geno.drop(geno.filter(regex=namequeen),axis=1)
	gn=geno.iloc[list(to_keep_final),:]
	gn.insert(0,'marker_id',marker_id)
	gn.to_csv('genotype_colony%d.txt'%(i+1),header=True,index=False,sep=' ')
	kept_remove=[None]*len(marker_id)	
	kr=numpy.array(kept_remove)
	kr[list(to_keep_final)]="keep"
	kr[list(m_error)]="sequencing_error"
	kr[list(m_multi)]="multi_allelic"
	kr[list(m_seq_nonok_tot)]="sequencing_error_resequenced"
	kr[list(m_het)]="drone_heterozygote"
	kr[list(m_homQ)]="queen_homozygote"
	kr[list(m_hom)]="non_polymorphic"
	keep_remove[i]=kr.tolist()
#	kr1=numpy.array(kept_remove)
#	kr1[list(to_keep_final)]="keep"
#	kr2=numpy.array(kept_remove)
#	kr2[list(m_hom)]="non_polymorphic"
#	kr3=numpy.array(kept_remove)
#	kr3[list(m_homQ)]="queen_homozygote"
#	kr4=numpy.array(kept_remove)
#	kr4[list(m_het)]="drone_heterozygote"
#	kr5=numpy.array(kept_remove)
#	kr5[list(m_seq_nonok_tot)]="sequencing_error_resequenced"
#	kr6=numpy.array(kept_remove)
#	kr6[list(m_multi)]="multi_allelic"
#	kr7=numpy.array(kept_remove)
#	kr7[list(m_error)]="sequencing_error"
#	keep_remove[i]=list(zip(kr1,kr2,kr3,kr4,kr5,kr6,kr7))
	print('Done colony %d'%(i+1))

k_df=pd.DataFrame.from_dict(keep_remove)
k_df.insert(0,'marker_id',marker_id)
k_df.to_csv('marker_aftercontrol.txt',header=True,index=False,sep=' ')

m_k=[]
for i in range(0,len(k_df)):
	if k_df.iloc[i,1]==k_df.iloc[i,2]==k_df.iloc[i,3]=='keep':
#	if k_df.iloc[i,1]==k_df.iloc[i,2]==k_df.iloc[i,3]=="('keep',None,None,None,None,None,None)":
		m_k.append(i)
geno=geno_f.iloc[list(m_k),:]
marker_id_all=geno.marker_id
print('number of markers in all colonies:',len(geno))

namequeen='Q'
for i in range(0,nb_colony):
	g=geno.filter(regex=colony_list[i])
	g=g.drop(g.filter(regex=namequeen),axis=1)
	g.insert(0,'marker_id',marker_id_all)
	g.to_csv('genotype_all_colony%d.txt'%(i+1),header=True,index=False,sep=' ')

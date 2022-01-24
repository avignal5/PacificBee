#!/usr/bin/env R
args<-commandArgs(TRUE)
nb_colony<-as.numeric(args[1])
pattern<-args[5]
library('plyr')
library('colorout')
library('reshape2')

chr_s<-read.table(args[2])
colnames(chr_s)<-c('contig_id','size')
chr_c<-read.table(args[3],header=F,sep='\t')
colnames(chr_c)<-c('chr_c','contig_id','pos','orientation_c')
chr_a<-read.table(args[4],header=T,sep=';')
colnames(chr_a)<-c('contig_id','chr_lg','orientation_a','cont_start','cont_end','cont_length','chr_start','chr_end','chr_length','remarques')
b<-colsplit(chr_a$chr_lg,'LG',c('LG_id','chr'))
chr_a$chr_a<-b$chr
chr<-merge(chr_c,chr_s,by='contig_id',all=TRUE)
chr<-merge(chr,chr_a,by='contig_id',all=TRUE)
for(i in 1:nrow(chr)){
	if(is.na(chr$chr_a[i]) & grepl('Un',chr$chr_c[i])==TRUE){chr$chr[i]<-'Un'
	}else if(chr$chr_c[i]==chr$chr_a[i]){chr$chr[i]<-chr$chr_a[i]
	}else if (chr$chr_c[i]!=chr$chr_a[i] & is.numeric(chr$chr_a[i])==TRUE & is.numeric(chr$chr_c[i])==FALSE){chr$chr[i]<-chr$chr_a[i]
	}else if (chr$chr_c[i]!=chr$chr_a[i] & is.numeric(chr$chr_a[i])==FALSE & is.numeric(chr$chr_c[i])==TRUE){chr$chr[i]<-chr$chr_c[i]}
	if(!is.na(chr$remarques[i])&chr$remarques[i]=='LG8 or LG10'){chr$chr[i]<-'Un'}}

geno<-args[6]
if(pattern=='tig'){
	chr$cont<-as.numeric(gsub('tig','',chr$contig_id))
	for(c in 1:nb_colony){
	col<-read.table(paste0(geno,c,'.txt'),header=T)
	cl<-unlist(strsplit(as.character(col$marker_id),'_'))
	cont_col<-unique(cl[c(T,F)])
	for(i in 1:nrow(chr)){if(chr$cont[i]%in%cont_col){chr[i,paste0('geno',c)]<-'y'}else{chr[i,paste0('geno',c)]<-'n'}}
	}
	chr_geno<-vector()
	for(i in 1:nrow(chr)){if(paste0(chr[i,grepl('geno',colnames(chr))],collapse='')=='yyy'){chr_geno<-c(chr_geno,'y')}else{chr_geno<-c(chr_geno,'n')}}
	chr$geno<-chr_geno	
}else if (pattern=='chr'){
	chr_a2<-unique(chr$chr_c)
	chr_a2[grepl('Un',chr_a2)]<-'Un'
	chr<-data.frame(contig_id=unique(chr$chr_c),chr_c=unique(chr$chr_c),pos=0,orientation_c=0,size=0,chr_lg=0,orientation_a=0,cont_start=0,cont_end=0,cont_length=0,chr_start=0,chr_end=0,chr_length=0,remarques=0,chr_a=chr_a2,chr=chr_a2)
	chr$cont<-gsub('Group','',unique(chr$chr_c))
	for(c in 1:nb_colony){
	col<-read.table(paste0(geno,c,'.txt'),header=T)
	cl<-unlist(strsplit(as.character(col$Position),'_'))
	cont_col<-unique(cl[c(T,F)])
	cont_col<-gsub('LG','',cont_col)
	for(i in 1:nrow(chr)){if(chr$cont[i]%in%cont_col){chr[i,paste0('geno',c)]<-'y'}else{chr[i,paste0('geno',c)]<-'n'}}
	}
	chr_geno<-vector()
	for(i in 1:nrow(chr)){if(paste0(chr[i,grepl('geno',colnames(chr))],collapse='')=='yyy'){chr_geno<-c(chr_geno,'y')}else{chr_geno<-c(chr_geno,'n')}}
	chr$geno<-chr_geno
}
write.table(chr,'ContigsOnChromosomes_complete.csv',sep=',',col.names=T,row.names=F,quote=F)

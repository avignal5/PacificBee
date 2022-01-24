################################################
### script to plot mapping qualities
################################################
#!/usr/bin/env R
library(data.table)
library(ggplot2)
library(reshape2)
library(tidyr)

pdf('ind_mapq_chr.pdf',width=20,height=10)
dat<-fread('mapping_qual.out',sep='|',data.table=F,header=F)

n_AMelMel<-vector()
n_HAV3<-vector()

for(i in 1:nrow(dat)){
	if(dat[i,]=='AMelMel'){n_AMelMel<-c(n_AMelMel,i)
	}else if(dat[i,]=='HAV3'){n_HAV3<-c(n_HAV3,i)}
}

ind<-dat[n_AMelMel-1,]
average_mapq_AMelMel<-data.frame('ind'=ind,'info'=dat[n_AMelMel+5,])
average_mapq_HAV3<-data.frame('ind'=ind[!(dat[n_HAV3+1,]%in%ind)],'info'=dat[n_HAV3+5,][!(dat[n_HAV3+1,]%in%ind)])
D_ind<-merge(average_mapq_AMelMel,average_mapq_HAV3,by='ind',all=T)
colnames(D_ind)<-c('ind','info_AMelMel','info_HAV3')
D_ind$nb_map_AMelMel<-NA
D_ind$mapq_AMelMel<-NA
D_ind$nb_map_HAV3<-NA
D_ind$mapq_HAV3<-NA
for(i in 1:nrow(D_ind)){
	D_ind$nb_map_AMelMel[i]<-unlist(strsplit(as.character(D_ind$info_AMelMel[i]),' '))[1]	
	D_ind$mapq_AMelMel[i]<-gsub('\\(','',unlist(strsplit(as.character(D_ind$info_AMelMel[i]),' '))[5])
	D_ind$nb_map_HAV3[i]<-unlist(strsplit(as.character(D_ind$info_HAV3[i]),' '))[1]	
	D_ind$mapq_HAV3[i]<-gsub('\\(','',unlist(strsplit(as.character(D_ind$info_HAV3[i]),' '))[5])}
D_ind$mapq_AMelMel<-gsub('%','',D_ind$mapq_AMelMel)
D_ind$mapq_HAV3<-gsub('%','',D_ind$mapq_HAV3)
D_ind$info_AMelMel<-NULL
D_ind$info_HAV3<-NULL
a<-colsplit(D_ind$ind,'_',c('ind','x','y'))
D_ind$ind<-a$ind
D_g<-D_ind[,c('ind','mapq_AMelMel','mapq_HAV3')]
D_g<-gather(D_g,version,mapq,mapq_AMelMel:mapq_HAV3)

plot(D_ind$nb_map_AMelMel,D_ind$nb_map_HAV3,xlim=c(0,60000000),ylim=c(0,60000000),pch=16,xlab='number of mapped AMelMel',ylab='number of mapped HAV3')
abline(0,1,col='red')
  
ggplot(D_g,aes(x=ind,y=as.numeric(mapq),col=version))+
	geom_point()+
	geom_line(aes(group=ind),col='grey',linetype='dashed')+
	theme(axis.text.x=element_text(angle=90))+
	xlab('individual')+
	ylab('average mapping quality')

cm_AMelMel<-unique(dat[grep('CM01',dat[1:nrow(dat),]),])
cm_HAV3<-unique(dat[grep('CM00',dat[1:nrow(dat),]),])
cm<-cbind('cm_AMelMel'=c(cm_AMelMel,NA),cm_HAV3)

M_region<-list()
H_region<-list()
for(x in 1:length(ind)){
print(ind[x])
ind_dat<-dat[n_AMelMel[x]:n_HAV3[x],]
cA<-vector()
for(j in 1:length(cm_AMelMel)){cA[j]<-grep(cm_AMelMel[j],ind_dat)}
mapq_region<-list()
for(j in 1:length(cA)){
if(j==length(cA)){
d<-ind_dat[cA[j]:length(ind_dat)]
n_mean<-grep('Mean',d)
n_pos<-n_mean-1
mapq_region[[j]]<-data.frame('pos'=rep(NA,length(n_pos)),'mean'=rep(NA,length(n_mean)),'pos_alt'=rep(NA,length(n_pos)))
mapq_region[[j]][,'pos']=paste0(cm_AMelMel[j],'_',d[n_pos])
mapq_region[[j]][,'mean']=gsub("[^0-9.]", "",d[n_mean])
mapq_region[[j]][,'pos_alt']=paste0(cm_HAV3[j],'_',d[n_pos])
}else{
d<-ind_dat[cA[j]:cA[j+1]]
n_mean<-grep('Mean',d)
n_pos<-n_mean-1
mapq_region[[j]]<-data.frame('pos'=rep(NA,length(n_pos)),'mean'=rep(NA,length(n_mean)),'pos_alt'=rep(NA,length(n_pos)))
mapq_region[[j]][,'pos']=paste0(cm_AMelMel[j],'_',d[n_pos])
mapq_region[[j]][,'mean']=gsub("[^0-9.]", "",d[n_mean])
mapq_region[[j]][,'pos_alt']=paste0(cm_HAV3[j],'_',d[n_pos])
}}
M_region[[x]]<-do.call(rbind,mapq_region)
print(nrow(M_region[[x]]))

if(x==length(ind)){ind_dat<-dat[n_HAV3[x]:nrow(dat),]
}else{ind_dat<-dat[n_HAV3[x]:n_AMelMel[x+1],]
}
if(length(ind_dat)>3){
cH<-vector()
for(j in 1:length(cm_HAV3)){cH[j]<-grep(cm_HAV3[j],ind_dat)}
mapq_region<-list()
for(j in 1:length(cH)){
if(j==length(cH)){
d<-ind_dat[cH[j]:length(ind_dat)]
n_mean<-grep('Mean',d)
n_pos<-n_mean-1
mapq_region[[j]]<-data.frame('pos'=rep(NA,length(n_pos)),'mean'=rep(NA,length(n_mean)),'pos_alt'=rep(NA,length(n_pos)))
mapq_region[[j]][,'pos_alt']=paste0(cm_HAV3[j],'_',d[n_pos])
mapq_region[[j]][,'mean']=gsub("[^0-9.]", "",d[n_mean])
mapq_region[[j]][,'pos']=paste0(cm_AMelMel[j],'_',d[n_pos])
}else{d<-ind_dat[cH[j]:cH[j+1]]
n_mean<-grep('Mean',d)
n_pos<-n_mean-1
mapq_region[[j]]<-data.frame('pos'=rep(NA,length(n_pos)),'mean'=rep(NA,length(n_mean)),'pos_alt'=rep(NA,length(n_pos)))
mapq_region[[j]][,'pos_alt']=paste0(cm_HAV3[j],'_',d[n_pos])
mapq_region[[j]][,'mean']=gsub("[^0-9.]", "",d[n_mean])
mapq_region[[j]][,'pos']=paste0(cm_AMelMel[j],'_',d[n_pos])
}}
H_region[[x]]<-do.call(rbind,mapq_region)
print(nrow(H_region[[x]]))
}else{H_region[[x]]<-data.frame('pos'='not_done','mean'='not_done','pos_alt'='not_done')}
}

D<-list()
for(x in 1:length(ind)){D[[x]]<-merge(M_region[[x]],H_region[[x]],by=c('pos','pos_alt'),all=T)
colnames(D[[x]])<-c('pos_AMelMel','pos_HAV3','mapq_AMelMel','mapq_HAV3')
D[[x]]$ind<-ind[x]
}
D_all<-do.call(rbind,D)
a<-colsplit(D_all$pos_AMelMel,'_',c('chr','pos'))
a2<-colsplit(a$pos,' ',c('start','stop'))
D_all$chr<-a$chr
D_all$pos_start<-as.numeric(a2$start)

for(x in 1:length(ind)){
	D_plot<-D_all[D_all$ind==ind[x],]
	D_g<-gather(D_plot,version,mapq,mapq_AMelMel:mapq_HAV3)	
print(ggplot(D_g,aes(x=pos_start,y=as.numeric(mapq),col=version))+
	geom_point()+
	geom_line(aes(group=pos_start),col='grey',linetype='dashed')+
	theme(axis.text.x=element_text(angle=90))+
	facet_wrap(~chr,scale='free')+
	ggtitle(unlist(strsplit(unique(D_g$ind),'_'))[1])+
	xlab('position on chromosome')+
	ylab('average mapping quality'))
}
dev.off()

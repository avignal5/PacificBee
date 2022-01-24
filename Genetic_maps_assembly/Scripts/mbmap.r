################################################
### script to plot recombination rates along chromosomes
################################################
#!/usr/bin/env R
args<-commandArgs(TRUE)
library('ggplot2')
library('gridExtra')
library('cowplot')
library('data.table')
library('reshape2')
library('colorout')

file_name<-args[1]
size_bin<-gsub('map','',unlist(strsplit(file_name,"_|.txt"))[2])
mbmap<-read.table(file_name,header=T,sep=' ')
chrmap<-read.table(args[2],header=T,sep=',')
gmap<-read.table(args[3],header=T,sep=',')
nchr<-max(mbmap$chr)

pdf(paste0('recomb_rate',size_bin,'.pdf'),width=20,height=20)
par(mfrow=c(4,4))
for(z in 1:nchr){
mbmapz<-mbmap[mbmap$chr==z,]
gmapz<-gmap[gmap$chr==z,]
a<-table((gmapz$pos_d+gmapz$pos_g)/2)
plot(mbmapz$left,mbmapz$m_cj,type='s',ylim=c(0,max(mbmap$q95_cj)),main=paste0('chr',z),xlim=c(0,max(gmapz$pos_d,mbmapz$right)))
segments(as.numeric(names(a)),rep(0,length(a)),as.numeric(names(a)),as.vector(a))
#points(mbmapz$right,mbmapz$q5_cj,type='s',col='blue',lty=3)
#points(mbmapz$right,mbmapz$q95_cj,type='s',col='blue',lty=3)
}

max(mbmap$m_cj)
min(mbmap$m_cj)

dat_bin<-list()
for(z in 1:nchr){
mbmapz<-mbmap[mbmap$chr==z,]
gmapz<-gmap[gmap$chr==z,]
chrmapz<-chrmap[chrmap$chr==z,]
nbbin<-nrow(mbmapz)
nb_cross<-vector()
nb_snp<-vector()
for(i in 1:nbbin){
	ci<-mbmapz[i,'left']
	cj<-mbmapz[i,'right']
	gmapi<-gmapz[gmapz$pos_g>ci & gmapz$pos_d<cj,]
	nb_cross<-c(nb_cross,nrow(gmapi))
	chrmapi<-chrmapz[chrmapz$pos>ci & chrmapz$pos<cj,]
	nb_snp<-c(nb_snp,nrow(chrmapi))}
dat_bin[[z]]<-data.frame(chr=rep(z,nbbin),pos_left=mbmapz$left,pos_right=mbmapz$right,nb_cross,nb_snp)	
}	
dat_bin<-do.call(rbind,dat_bin)
par(mfrow=c(1,1))
hist(dat_bin$nb_cross,breaks=50)
hist(dat_bin$nb_snp,breaks=50)
dev.off()

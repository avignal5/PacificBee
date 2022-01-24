################################################
### script to summarise recombination events
################################################
#!/usr/bin/env R
args<-commandArgs(TRUE)
nb_colony<-as.numeric(args[1])
kb<-as.numeric(args[2])
pattern<-args[3]
data<-read.table(args[4],sep=',',header=TRUE)
library('ggplot2')
library('gridExtra')
library('cowplot')
library('data.table')
library('reshape2')
library('colorout')
b<-as.numeric(as.character(data$chr_nb_consensus))
nchr<-max(b[!is.na(b)])

options(warn=2)
reco<-list()
info_all<-list()
reco_mean<-matrix(nrow=nchr+1,ncol=nb_colony)
ndata<-list()
nintercross<-list()
for(c in 1:nb_colony){
g<-subset(data,data[paste0('geno',c)]=='y')
print(paste0('colony',c,' size placed on chromosome: ',sum(g$size)))
print(paste0('colony',c,' size placed and oriented on chromosome: ',sum(subset(g$size,g[,paste0('orientation',c)]!=0))))
print(paste0('colony',c,' size placed on chromosome from genotyped: ',sum(g[paste0('size',c)])))
print(paste0('colony',c,' size placed and oriented on chromosome from genotyped: ',sum(subset(g[,paste0('size',c)],g[,paste0('orientation',c)]!=0))))
chr<-read.table(paste(pattern,c,'/chromosome/all_chr_new2.txt',sep=''),header=T,sep=',',colClasses=c('character','character','character','character','character','character','character','character','character','character','character','character','character','character','character','character'))
info<-list()
recomb<-data.frame()
for (i in 1:length(unique(chr$chromosome))){
datai<-subset(chr,chr$chromosome==i)
info_j<-list()
start<-0
delta<-50000
for(j in 1:length(unique(datai$contig))){
dataj<-subset(datai,datai$contig==unique(datai$contig)[j])
sizec<-data$size[data$cont==unique(datai$contig)[j]]
if(sizec==0){sizec<-max(as.numeric(dataj$pos_m),as.numeric(dataj$pos_M))}
if(nrow(dataj)==3 & dataj[1,'vector1']==dataj[3,'vector1']){orientation<-0}else if(as.numeric(dataj[1,'bin'])>as.numeric(dataj[nrow(dataj),'bin'])){orientation<-'-'}else if (as.numeric(dataj[1,'bin'])<as.numeric(dataj[nrow(dataj),'bin'])){orientation<-'+'}else if (nrow(dataj==1)){orientation<-0}
dataj$orientation<-rep(orientation,nrow(dataj))
if (orientation=='-'){
	p<-sizec
	for(m in 1:nrow(dataj)){
		p<-c(p,as.numeric(dataj[m,'pos_M']),as.numeric(dataj[m,'pos_m']))
		}
	pos<-vector()	
	for(q in 1:length(p)){
		if(q==1){pos<-p[q]-p[q+1]}else{pos<-c(pos,pos[q-1]+(p[q]-p[q+1]))}
		}
pos<-pos[1:length(pos)-1]
posStart<-pos[seq(1,length(pos),2)]
posStop<-pos[seq(2,length(pos),2)]
info_j[[j]]<-data.frame(dataj$marker1,posStart+start,posStop+start,posStart,posStop,dataj$chromosome,dataj$contig,dataj$orientation,do.call(rbind,strsplit(dataj[,'vector1'],'')))
start<-info_j[[j]][nrow(info_j[[j]]),3]+delta
}else if (orientation=='+'){
info_j[[j]]<-data.frame(dataj$marker1,as.numeric(dataj$pos_m)+start,as.numeric(dataj$pos_M)+start,as.numeric(dataj$pos_m),as.numeric(dataj$pos_M),dataj$chromosome,dataj$contig,dataj$orientation,do.call(rbind,strsplit(dataj[,'vector1'],'')))
start<-info_j[[j]][nrow(info_j[[j]]),3]+delta
} else if (orientation=='0' & as.numeric(dataj$pos_m[1])>as.numeric(dataj$pos_m[nrow(dataj)])){
dataj<-dataj[rev(rownames(dataj)), ]
info_j[[j]]<-data.frame(dataj$marker1,as.numeric(dataj$pos_m)+start,as.numeric(dataj$pos_M)+start,as.numeric(dataj$pos_m),as.numeric(dataj$pos_M),dataj$chromosome,dataj$contig,dataj$orientation,do.call(rbind,strsplit(dataj[,'vector1'],'')))
start<-info_j[[j]][nrow(info_j[[j]]),3]+delta
}else {info_j[[j]]<-data.frame(dataj$marker1,as.numeric(dataj$pos_m)+start,as.numeric(dataj$pos_M)+start,as.numeric(dataj$pos_m),as.numeric(dataj$pos_M),dataj$chromosome,dataj$contig,dataj$orientation,do.call(rbind,strsplit(dataj[,'vector1'],'')))
start<-info_j[[j]][nrow(info_j[[j]]),3]+delta}
colnames(info_j[[j]])<-c('marker','posStart','posStop','posStart_i','posStop_i','chr','contig','orientation',paste0('ind',seq(1:nchar(dataj$vector1[1]))))
for(k in 1:ncol(info_j[[j]])){info_j[[j]][,k]<-as.character(info_j[[j]][,k])}
}
info[[i]]<-do.call(rbind,info_j)
if(i>1){info[[i]]$posPlot<-as.numeric(info[[i]]$posStart)+(max(as.numeric(info[[i-1]]$posPlot))+1000000)}else{info[[i]]$posPlot<-as.numeric(info[[i]]$posStart)}
indgeno<-info[[i]][grepl('ind',colnames(info[[i]]))]
nrecom<-vector()
for(d in 1:ncol(indgeno)){
	nr<-0
	for(k in 1:(nrow(indgeno)-1)){
		if(indgeno[k,d]!=indgeno[k+1,d]){nr<-nr+1}
		}
	nrecom<-c(nrecom,nr)
	}
recomb<-rbind(recomb,nrecom)
reco_mean[i,c]<-mean(as.numeric(recomb[i,]))
}
reco[[c]]<-data.frame(recomb)
rownames(reco[[c]])<-paste0('chr',seq(1:nrow(reco[[c]])))
colnames(reco[[c]])<-paste0('ind',seq(1:ncol(reco[[c]])))
info_all[[c]]<-do.call(rbind,info)
info_all[[c]][is.na(info_all[[c]])]<-0
reco_mean[i+1,c]<-sum(reco_mean[1:i,c])
rownames(reco_mean)<-c(paste0('chr',seq(1:nrow(reco[[c]]))),'sum')
colnames(reco_mean)<-paste0('col',seq(1:nb_colony))
write.table(info_all[[c]],paste0('info_all',c,'.txt'),sep=',',col.names=T,row.names=F,quote=F)

indg<-ncol(info_all[[c]][grepl('ind',colnames(info_all[[c]]))])
ndatai<-list()
intercrossi<-list()
for(d in 1:indg){
ndat<-list()
intercross<-list()
for (z in 1:length(unique(info_all[[c]]$chr))){
b<-info_all[[c]][info_all[[c]]$chr==z,c('marker','posStart_i','posStop_i','chr','contig',paste0('ind',d))]
hap<-b[1,]
i=1
while(i<(nrow(b)-1)){
if(b[i,paste0('ind',d)]==b[(i+1),paste0('ind',d)] & b[i,'contig']==b[(i+1),'contig']){i=i+1}else{
hap<-rbind(hap,b[i,])
i=i+1
hap<-rbind(hap,b[i,])
}
}
hap<-rbind(hap,b[nrow(b),])
j=1
l<-vector()
while(j<nrow(hap)){
	if(hap$contig[j+1]==hap$contig[j]){
	l<-c(l,as.numeric(hap$posStop_i[j+1])-as.numeric(hap$posStart_i[j]))
	j=j+2}
	else{
	l<-c(l,((as.numeric(hap$posStop_i[j])-as.numeric(hap$posStart_i[j]))+(as.numeric(hap$posStop_i[j+1])-as.numeric(hap$posStart_i[j+1]))))
	j=j+2}
	}
ndat[[z]]<-cbind(colony=rep(c,nrow(hap)/2),chr=as.character(hap[1:(nrow(hap)/2),'chr']),contig=as.character(hap$contig[seq(nrow(hap)) %% 2 == 1]),ind=rep(paste0('ind',d),nrow(hap)/2),version=hap[,paste0('ind',d)][seq(nrow(hap)) %% 2 == 1],length=l)	
s=1	
k<-vector()	
ck<-vector()
while(s<(nrow(hap)-1)){
	if(hap$contig[s+2]==hap$contig[s] & hap[s+1,ncol(hap)]==hap[s,ncol(hap)]){
	k<-c(k,as.numeric(hap$posStart_i[s+2])-as.numeric(hap$posStart_i[s]))
	ck<-c(ck,as.character(hap$contig[s]))
	s=s+2}else{s=s+1}
	}
intercross[[z]]<-cbind(colony=rep(c,length(k)),chr=rep(hap$chr[1],length(k)),contig=ck,ind=rep(paste0('ind',d),length(k)),dist_inter_cross=k)
}
ndatai[[d]]<-do.call(rbind,ndat)
intercrossi[[d]]<-do.call(rbind,intercross)
}
ndata[[c]]<-do.call(rbind,ndatai)
nintercross[[c]]<-do.call(rbind,intercrossi)
}
reco
reco_mean
seg_recom<-do.call(rbind,ndata)
seg_intercross<-do.call(rbind,nintercross)
write.table(seg_recom,'recombination_segments.txt',sep=',',col.names=T,row.names=F,quote=F)
write.table(seg_intercross,'intercross_segments.txt',sep=',',col.names=T,row.names=F,quote=F)
seg_recom<-as.data.frame(seg_recom)
seg_recom$length<-as.numeric(as.character(seg_recom$length))
seg_intercross<-as.data.frame(seg_intercross)
seg_intercross$dist_inter_cross<-as.numeric(as.character(seg_intercross$dist_inter_cross))
pdf('segment.pdf',width=10,height=10)
par(mfrow=c(2,2))
hist(seg_recom$length,breaks=20)
hist(log10(seg_recom$length),breaks=20)
hist(seg_intercross$dist_inter_cross,breaks=20)
hist(log10(seg_intercross$dist_inter_cross),breaks=20)
dev.off()


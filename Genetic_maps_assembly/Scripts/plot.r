################################################
### script to plot number of recombination and summary information Liu vs AMelMel
################################################
#!/usr/bin/env R
args<-commandArgs(TRUE)
library('ggplot2')
library('gridExtra')
library('cowplot')
library('data.table')
library('reshape2')
library("seqinr")
slidingwindow<-function(windowsize, inputseq)
  {starts <- seq(1, length(inputseq), by = windowsize)
	starts[2:length(starts)]<-starts[2:length(starts)]-1
	starts<-c(starts,length(datseq))
	n<-length(starts)
   chunkGCs <- numeric(n)
   for (i in 1:n) {chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)]
       if(any(is.na(chunk))){chunk<-chunk[!is.na(chunk)]}
       chunkGC <- GC(chunk)
       chunkGCs[i] <- chunkGC}
datGC<-data.frame(starts,chunkGCs)}
nb_colony<-as.numeric(args[1])
#nb_colony<-3
parent_folder<-args[2]
#parent_folder<-'PacBio'
pattern<-args[3]
#pattern<-'group_colony'
parent_folderl<-args[4]
#parent_folderl<-'Liu'
patternl<-args[5]
#patternl<-'colony_group'
window1<-as.numeric(args[8])
#window1<-1000000
window2<-as.numeric(args[8])
#window2<-50000
data<-read.table(paste0(parent_folder,'/',args[7]),sep=',',header=TRUE)
#data<-read.table(paste0(parent_folder,'/','combined_info_colony.txt'),sep=',',header=TRUE)
b<-as.numeric(as.character(data$chr_nb_consensus))
nchr<-max(b[!is.na(b)])
agpFile = read.table("AMelMelPacBioChromosomes.agp",header=FALSE, sep = "\t")
blankPlot<-ggplot()+geom_blank(aes(1,1))+theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(), axis.text.y = element_blank(),axis.ticks = element_blank(),axis.line = element_blank())
options(scipen=999)

##### snp present in intersect of 3 colonies
type<-'all'
file_name<-paste0(args[6],type,'.txt')
#file_name<-paste0('Mb_map1_',type,'.txt')
size_bin<-gsub('map','',unlist(strsplit(file_name,"_|.txt"))[2])
mbmap<-read.table(paste0(parent_folder,'/',file_name),header=T,sep=' ')
gmap<-read.table(paste0(parent_folder,'/','genetic_map_',type,'.txt'),header=T,sep=',')
chrmap<-read.table(paste0(parent_folder,'/','chromosome_map_',type,'.txt'),header=T,sep=',')
mbmapl<-read.table(paste0(parent_folderl,'/',file_name),header=T,sep=' ')
gmapl<-read.table(paste0(parent_folderl,'/','genetic_map_',type,'.txt'),header=T,sep=',')
chrmapl<-read.table(paste0(parent_folderl,'/','chromosome_map_',type,'.txt'),header=T,sep=',')

chr<-read.table(paste0(parent_folder,'/genotype_allphased.txt'),header=F,sep=',',colClasses=c('numeric','numeric','numeric','character'))
colnames(chr)<-c('chromosome','contig','snp','vector')
datx<-list()
datd<-list()
#datr<-list()
cont_pos<-list()
for(z in 1:nchr){
cont<-data[data$chr_nb_consensus==z,'cont']
orient<-data[data$chr_nb_consensus==z,'orientation_final']
data_x<-list()
for(x in 1:length(cont)){
	datac<-chr[chr$contig==cont[x],]
	if(nrow(datac)==0){x=x+1}else{
	bp_max<-data[data$chr_nb_consensus==z & data$cont==cont[x],'size']
	if(orient[x]=='+'){if(datac$snp[1]>datac$snp[nrow(datac)]){datac<-datac[nrow(datac):1,]}
	}else if (orient[x]=='-'){if(datac$snp[1]<datac$snp[nrow(datac)]){datac<-datac[nrow(datac):1,]}}
	rec<-vector()
	j<-1
	while(j<nrow(datac)){
	vj<-datac$vector[j]
	vk<-datac$vector[j+1]
	r<-0
	for (k in 1:nchar(vj)){if(substr(vj,k,k)!=substr(vk,k,k)){r<-r+1}}
	rec<-c(rec,r)
	j=j+1}
	datac$rec<-c(rec,0)
	datac<-datac[,c('chromosome','contig','snp','rec')]
	if(orient[x]=='-'){datac$snp<-bp_max-datac$snp}	
	if(datac[1,'snp']!=1){d_top<-data.frame('chromosome'=z,'contig'=cont[x],'snp'=1,'rec'=0)
		datac<-rbind(d_top,datac)}
	if(datac[nrow(datac),'snp']!=bp_max){d_bot<-data.frame('chromosome'=z,'contig'=cont[x],'snp'=bp_max,'rec'=datac[nrow(datac),'rec'])
		datac<-rbind(datac,d_bot)}
	data_x[[x]]<-datac
}}
datx[[z]]<-do.call(rbind,data_x)		
for(i in 1:nrow(datx[[z]])){n<-which(cont==datx[[z]][i,'contig'])
	size_add<-sum(data[data$chr_nb_consensus==z & data$cont%in%cont[1:n-1],'size'])
	if(datx[[z]]$contig[i]!=cont[1]){datx[[z]]$snp[i]<-datx[[z]]$snp[i]+size_add}}	
value='snp'
inputseq=datx[[z]]
starts<-c(seq(1,max(inputseq$snp),by=window1),max(inputseq$snp))
n<-length(starts)
snps<-numeric(n)
for(i in 1:(n-1)){
	if(i==(n-1)){chunk<-subset(inputseq[,paste0(value)],inputseq[,paste0(value)]>=starts[i] & inputseq[,paste0(value)]<=starts[i+1])
	}else{chunk<-subset(inputseq[,paste0(value)],inputseq[,paste0(value)]>=starts[i] & inputseq[,paste0(value)]<starts[i+1])} 
	if(any(is.na(chunk))){chunk<-chunk[!is.na(chunk)]}
    sum_<-length(chunk)
    snps[i]<-sum_
     }
datd[[z]]<-data.frame(pos=starts,snp=snps)
#value='rec'
#inputseq=datx[[z]]
#starts<-c(seq(1,max(inputseq$snp),by=window2),max(inputseq$snp))
#n<-length(starts)
#recs<-numeric(n)
#for(i in 1:(n-1)){
#	if(i==(n-1)){chunk<-subset(inputseq[,paste0(value)],inputseq[,'snp']>=starts[i] & inputseq[,'snp']<=starts[i+1])
#	}else{chunk<-subset(inputseq[,paste0(value)],inputseq[,'snp']>=starts[i] & inputseq[,'snp']<starts[i+1])} 
#	if(any(is.na(chunk))){chunk<-chunk[!is.na(chunk)]}
#    sum_<-sum(chunk)
#    recs[i]<-sum_
#     }
#datr[[z]]<-data.frame(pos=starts,rec=recs)
c_pos<-vector()
for(i in 1:length(unique(datx[[z]]$contig))){
	s<-datx[[z]][datx[[z]]$contig==unique(datx[[z]]$contig)[i],'snp']
	c_pos<-c(c_pos,s[1])}
cont_pos[[z]]<-c_pos
}

chrLiu<-read.table(paste0(parent_folderl,'/colony_allphased.txt'),header=F,sep=',',colClasses=c('numeric','numeric','numeric','numeric','character'))
colnames(chrLiu)<-c('chromosome','contig','snp','order','vector')
datxl<-list()
datdl<-list()
#datrl<-list()
for(z in 1:nchr){
	s<-chrLiu[chrLiu$chromosome==z,]
	rec<-vector()
	j<-1
	while(j<nrow(s)){
	vj<-s$vector[j]
	vk<-s$vector[j+1]
	r<-0
	for (k in 1:nchar(vj)){
	if(substr(vj,k,k)!=substr(vk,k,k)){r<-r+1}
	}
	rec<-c(rec,r)
	j=j+1}
datxl[[z]]<-data.frame(snp=s$snp,rec=c(rec,0))
value='snp'
inputseq=datxl[[z]]
starts<-c(seq(1,max(inputseq$snp),by=window1),max(inputseq$snp))
n<-length(starts)
snps<-numeric(n)
for(i in 1:(n-1)){
	if(i==(n-1)){chunk<-subset(inputseq[,paste0(value)],inputseq[,paste0(value)]>=starts[i] & inputseq[,paste0(value)]<=starts[i+1])
	}else{chunk<-subset(inputseq[,paste0(value)],inputseq[,paste0(value)]>=starts[i] & inputseq[,paste0(value)]<starts[i+1])} 
	if(any(is.na(chunk))){chunk<-chunk[!is.na(chunk)]}
    sum_<-length(chunk)
    snps[i]<-sum_
     }
datdl[[z]]<-data.frame(pos=starts,snp=snps)
#value='rec'
#inputseq=datxl[[z]]
#starts<-c(seq(1,max(inputseq$snp),by=window2),max(inputseq$snp))
#n<-length(starts)
#recs<-numeric(n)
#for(i in 1:(n-1)){
#	if(i==(n-1)){chunk<-subset(inputseq[,paste0(value)],inputseq[,'snp']>=starts[i] & inputseq[,'snp']<=starts[i+1])
#	}else{chunk<-subset(inputseq[,paste0(value)],inputseq[,'snp']>=starts[i] & inputseq[,'snp']<starts[i+1])} 
#	if(any(is.na(chunk))){chunk<-chunk[!is.na(chunk)]}
#    sum_<-sum(chunk)
#    recs[i]<-sum_
#     }
#datrl[[z]]<-data.frame(pos=starts,rec=recs)
}

pdf('PacBiovsLiu_all.pdf',width=25,height=25)
for(z in 1:nchr){
		print(z)
png(paste0('PacBiovsLiu_chr',z,'.png'),width=18,height=18,units='in',res=600)
n<-which(datxl[[z]]$rec>=5)
if(length(n)>0){lineslr<-vector()
linesll<-vector()	
for(i in 1:length(n)){
	s<-datxl[[z]][n[i]:(n[i]+1),]
	lineslr<-c(lineslr,s$snp[1])
	linesll<-c(linesll,s$snp[2])
	}}else{lineslr=-100
	linesll=-100}
lineL<-data.frame(lineslr,linesll,grp=seq(1:length(lineslr)))
n<-which(datx[[z]]$rec>=5)
if(length(n)>0){linesr<-vector()
linesl<-vector()	
for(i in 1:length(n)){s<-datx[[z]][n[i]:(n[i]+1),]
	linesr<-c(linesr,s$snp[1])
	linesl<-c(linesl,s$snp[2])
	}}else{linesr=-100
	linesl=-100}
line<-data.frame(linesr,linesl,grp=seq(1:length(linesr)))
chrT= read.table(paste("Chr",z,"_split1_chained.csv",sep=""),header = TRUE)
chromosome = paste("chr",z,sep="")
contigLimits1 = agpFile[agpFile$V1 == chromosome,]$V2
contigLimits2 = agpFile[agpFile$V1 == chromosome,]$V3
a<-seq(1:(nrow(chrT)/3))
grp<-vector()
for(j in 1:length(a)){grp<-c(grp,rep(a[j],3))}
chrT$grp<-grp
size1<-max(chrT$query[!is.na(chrT$query)])
size2<-max(chrT$target[!is.na(chrT$query)])
p1<-ggplot(chrT,aes(query,target,group=grp))+
	geom_point(shape='.')+
	geom_line(size=1)+
	xlab("AMelMel (position in bp)")+
	xlim(0,size1)+
	ylim(0,size2)+
	ylab("Amel4.5 (position in bp)")+
	geom_hline(yintercept=lineL$linesll,size=0.2,color='red')+
	geom_hline(yintercept=lineL$lineslr,size=0.2,color='red')+
	geom_vline(xintercept=line$linesl,size=0.2,color='red')+
	geom_vline(xintercept=line$linesr,size=0.2,color='red')+
	geom_vline(xintercept=cont_pos[[z]],color='black',linetype=3,size=0.8)+
	ggtitle(paste0("Chromosome ",z))+
	theme(plot.title = element_text(size = 20, face = "bold"))+
	theme(axis.text.y=element_text(angle=90,hjust=1))+
	theme(plot.margin=unit(c(2,2.5,1,1),'cm'))+
	geom_rect(data=lineL,inherit.aes=F,aes(xmin=-Inf,xmax=+Inf,ymin=lineslr,ymax=linesll,group=grp),fill='red',alpha=0.2)+
	geom_rect(data=line,inherit.aes=F,aes(ymin=-Inf,ymax=+Inf,xmin=linesr,xmax=linesl,group=grp),fill='red',alpha=0.2)	
p2<-ggplot()+
	geom_step(data=datx[[z]],aes(x=snp,y=rec))+
	xlab('')+
	ylab('Recombination')+
	xlim(0,size1)+
	geom_vline(xintercept=cont_pos[[z]],color='black',linetype=3,size=0.8)+
	geom_vline(xintercept=line$linesl,color='red',size=0.2)+
	geom_vline(xintercept=line$linesr,color='red',size=0.2)+
	geom_rect(data=line,inherit.aes=F,aes(ymin=-Inf,ymax=+Inf,xmin=linesr,xmax=linesl,group=grp),fill='red',alpha=0.2)+
	theme(plot.margin=unit(c(0,2.4,0,1),'cm'))+
	theme(legend.position='none')+
	ylim(0,20)+
	theme(axis.title.y.right=element_text(color='blue'))
p3<-ggplot()+
	geom_step(data=datd[[z]],aes(x=pos,y=snp))+
	xlab('')+
	xlim(0,size1)+
	geom_vline(xintercept=cont_pos[[z]],color='black',linetype=3,size=0.8)+
	geom_vline(xintercept=line$linesl,color='red',size=0.2)+
	geom_vline(xintercept=line$linesr,color='red',size=0.2)+
	geom_rect(data=line,inherit.aes=F,aes(ymin=-Inf,ymax=+Inf,xmin=linesr,xmax=linesl,group=grp),fill='red',alpha=0.2)+
	theme(plot.margin=unit(c(0,2.4,0,0.5),'cm'))+
	theme(legend.position='none')+
	ylab('SNP density')
mbmapz<-mbmap[mbmap$chr==z,]
mbmapz<-rbind(mbmapz,c(z,mbmapz[nrow(mbmapz),'right'],size1,mbmapz[nrow(mbmapz),'m_cj'],mbmapz[nrow(mbmapz),'s_cj'],mbmapz[nrow(mbmapz),'q5_cj'],mbmapz[nrow(mbmapz),'q95_cj']))
p4<-ggplot(mbmapz)+
	geom_step(aes(left,m_cj))+
	geom_hline(aes(yintercept=23),color='blue',linetype='dashed')+
	xlab('')+
	xlim(0,size1)+
	ylim(0,max(mbmap$m_cj))+
	geom_vline(xintercept=cont_pos[[z]],size=0.8,linetype=3)+
	geom_vline(xintercept=line$linesl,color='red',size=0.2)+
	geom_vline(xintercept=line$linesr,color='red',size=0.2)+
	geom_rect(data=line,inherit.aes=F,aes(ymin=-Inf,ymax=+Inf,xmin=linesr,xmax=linesl,group=grp),fill='red',alpha=0.2)+
	theme(plot.margin=unit(c(0,2.4,0,1),'cm'))+
	ylab('Recombination rate cm/Mb')
dat<-read.fasta(file=paste0("Fasta/AMelMelPacBio_chr",z,".fa"))
datseq <- dat[[1]]
datGC<-slidingwindow(window1, datseq)
colnames(datGC)<-c('pos','GC')
cov_mean<-list()
for(c in 1:nb_colony){
cov<-read.table(paste0('PacBio/samCovQueen',c,'_',z,'_sub.txt'))
colnames(cov)<-c('chr','pos','cov')
pos_mean<-vector()
a=1
b=window1
c_mean<-vector()
while(b<nrow(cov)){
c_mean<-c(c_mean,mean(cov$cov[a:b]))
pos_mean<-c(pos_mean,a)
a=b
b=b+window1
}
c_mean<-c(c_mean,mean(cov$cov[tail(pos_mean,1):a]),mean(cov$cov[a:nrow(cov)]))
pos_mean<-c(pos_mean,a,nrow(cov))
cov_mean[[c]]<-c_mean
}
cov2<-data.frame(pos_mean,do.call(cbind,cov_mean))
colnames(cov2)<-c('pos','c1','c2','c3')
dat_GCcov<-merge(datGC,cov2,by='pos')
p5<-ggplot(dat_GCcov)+
	geom_step(aes(x=pos,y=GC*100))+ 
	xlab('')+
	xlim(0,size1)+
	ylim(0,100)+
	geom_vline(xintercept=cont_pos[[z]],size=0.8,linetype=3)+
	geom_vline(xintercept=line$linesl,color='red',size=0.2)+
	geom_vline(xintercept=line$linesr,color='red',size=0.2)+
	geom_rect(data=line,inherit.aes=F,aes(ymin=-Inf,ymax=+Inf,xmin=linesr,xmax=linesl,group=grp),fill='red',alpha=0.2)+
	theme(plot.margin=unit(c(0,0.7,0,0.7),'cm'))+
	geom_step(aes(x=pos,y=c1),linetype='longdash',color='blue')+ 
	geom_step(aes(x=pos,y=c2),linetype='dotdash',color='blue')+ 
	geom_step(aes(x=pos,y=c3),linetype='dotted',color='blue')+ 
	scale_y_continuous(name='GC %',sec.axis=(sec_axis(~.*1,name='Queen coverage')))+	
	theme(axis.title.y.right=element_text(color='blue'))
p6<-ggplot()+
	geom_step(data=datxl[[z]],aes(x=snp,y=rec))+
	xlab('')+
	ylab('Recombination')+
	ylim(0,20)+
	xlim(0,size2)+
	geom_vline(xintercept=lineL$linesll,color='red',size=0.2)+
	geom_vline(xintercept=lineL$lineslr,color='red',size=0.2)+
	geom_rect(data=lineL,inherit.aes=F,aes(ymin=-Inf,ymax=+Inf,xmin=lineslr,xmax=linesll,group=grp),fill='red',alpha=0.2)+
	theme(plot.margin=unit(c(3,0,1,0),'cm'))+
	theme(legend.position='none')+
	theme(axis.text.y=element_text(angle=90))+
	coord_flip()
p7<-ggplot()+
	geom_step(data=datdl[[z]],aes(x=pos,y=snp))+
	xlab('')+
	ylab('SNP density')+
	xlim(0,size2)+
	geom_vline(xintercept=lineL$linesll,color='red',size=0.2)+
	geom_vline(xintercept=lineL$lineslr,color='red',size=0.2)+
	geom_rect(data=lineL,inherit.aes=F,aes(ymin=-Inf,ymax=+Inf,xmin=lineslr,xmax=linesll,group=grp),fill='red',alpha=0.2)+
	theme(plot.margin=unit(c(3,0,1,0),'cm'))+
	theme(legend.position='none')+
	theme(axis.text.y=element_text(angle=90))+
	coord_flip()
mbmaplz<-mbmapl[mbmapl$chr==z,]
p8<-ggplot(mbmaplz)+
	geom_step(aes(left,m_cj))+
	geom_hline(aes(yintercept=37),color='blue',linetype='dashed')+
	geom_vline(xintercept=lineL$linesll,color='red',size=0.2)+
	geom_vline(xintercept=lineL$lineslr,color='red',size=0.2)+
	geom_rect(data=lineL,inherit.aes=F,aes(ymin=-Inf,ymax=+Inf,xmin=lineslr,xmax=linesll,group=grp),fill='red',alpha=0.2)+
	xlab('')+
	xlim(0,size2)+
	ylim(0,160)+
	theme(plot.margin=unit(c(3,0,1,0),'cm'))+
	ylab('Recombination rate cm/Mb')+
	theme(axis.text.y=element_text(angle=90))+
	coord_flip()
dat<-read.fasta(file=paste0("Fasta/ame_ref_Amel_4.5_chrLG",z,".fa"))
datseq <- dat[[1]]
datGC<-slidingwindow(window1, datseq)
p9<-ggplot(datGC)+
	geom_step(aes(x=starts,y=chunkGCs*100))+ 
	xlab('')+
	xlim(0,size2)+
	ylim(0,100)+
	geom_vline(xintercept=lineL$linesll,color='red',size=0.2)+
	geom_vline(xintercept=lineL$lineslr,color='red',size=0.2)+
	geom_rect(data=lineL,inherit.aes=F,aes(ymin=-Inf,ymax=+Inf,xmin=lineslr,xmax=linesll,group=grp),fill='red',alpha=0.2)+
	theme(plot.margin=unit(c(3,1,1,0),'cm'))+
	ylab('GC %')+
	theme(axis.text.y=element_text(angle=90))+
	coord_flip()
grid.arrange(p9,p8,p7,p6,p1,
			blankPlot,blankPlot,blankPlot,blankPlot,p2,
			blankPlot,blankPlot,blankPlot,blankPlot,p3,
			blankPlot,blankPlot,blankPlot,blankPlot,p4,
			blankPlot,blankPlot,blankPlot,blankPlot,p5,
			ncol=5,nrow=5,widths=c(1,1,1,1.5,3),heights=c(3,1.5,1,1,1))	
dev.off()
}
dev.off()

#figures for paper#
png('PacBiovsLiu_chr3Paper.png',width=12,height=12,units='in',res=600)
z=3
n<-which(datxl[[z]]$rec>=5)
if(length(n)>0){lineslr<-vector()
linesll<-vector()	
for(i in 1:length(n)){
	s<-datxl[[z]][n[i]:(n[i]+1),]
	lineslr<-c(lineslr,s$snp[1])
	linesll<-c(linesll,s$snp[2])
	}}else{lineslr=-100
	linesll=-100}
lineL<-data.frame(lineslr,linesll,grp=seq(1:length(lineslr)))
n<-which(datx[[z]]$rec>=5)
if(length(n)>0){linesr<-vector()
linesl<-vector()	
for(i in 1:length(n)){s<-datx[[z]][n[i]:(n[i]+1),]
	linesr<-c(linesr,s$snp[1])
	linesl<-c(linesl,s$snp[2])
	}}else{linesr=-100
	linesl=-100}
line<-data.frame(linesr,linesl,grp=seq(1:length(linesr)))
chrT= read.table(paste("Chr",z,"_split1_chained.csv",sep=""),header = TRUE)
chromosome = paste("chr",z,sep="")
contigLimits1 = agpFile[agpFile$V1 == chromosome,]$V2
contigLimits2 = agpFile[agpFile$V1 == chromosome,]$V3
a<-seq(1:(nrow(chrT)/3))
grp<-vector()
for(j in 1:length(a)){grp<-c(grp,rep(a[j],3))}
chrT$grp<-grp
size1<-max(chrT$query[!is.na(chrT$query)])
size2<-max(chrT$target[!is.na(chrT$query)])
p1<-ggplot(chrT,aes(query,target,group=grp))+
	geom_point(shape='.')+
	geom_line(size=0.3)+
	xlab("AMelMel (position in bp)")+
	xlim(0,size1)+
	ylim(0,size2)+
	ylab("Amel4.5 (position in bp)")+
	geom_vline(xintercept=cont_pos[[z]],size=0.8,linetype=3)+
	geom_hline(yintercept=lineL$linesll,size=0.2,color='red')+
	geom_hline(yintercept=lineL$lineslr,size=0.2,color='red')+
	geom_vline(xintercept=line$linesl,size=0.2,color='red')+
	geom_vline(xintercept=line$linesr,size=0.2,color='red')+
	ggtitle(paste0("Chromosome ",z))+
	theme(plot.title = element_text(size = 20, face = "bold"))+
	theme(axis.text.y=element_text(angle=90,hjust=1))+
	theme(plot.margin=unit(c(2,2.5,1,0),'cm'))+
	geom_rect(data=lineL,inherit.aes=F,aes(xmin=-Inf,xmax=+Inf,ymin=lineslr,ymax=linesll,group=grp),fill='red',alpha=0.2)+
	geom_rect(data=line,inherit.aes=F,aes(ymin=-Inf,ymax=+Inf,xmin=linesr,xmax=linesl,group=grp),fill='red',alpha=0.2)	
p2<-ggplot()+
	geom_step(data=datx[[z]],aes(x=snp,y=rec))+
	xlab('')+
	ylab('Recombination')+
	xlim(0,size1)+
	ylim(0,20)+
	geom_vline(xintercept=cont_pos[[z]],color='black',linetype=3,size=0.8)+
	geom_vline(xintercept=line$linesl,color='red',size=0.2)+
	geom_vline(xintercept=line$linesr,color='red',size=0.2)+
	geom_rect(data=line,inherit.aes=F,aes(ymin=-Inf,ymax=+Inf,xmin=linesr,xmax=linesl,group=grp),fill='red',alpha=0.2)+
	theme(plot.margin=unit(c(0,2.4,0,0),'cm'))+
	theme(axis.text.y = element_text(angle = 90, hjust = 1))+
	theme(legend.position='none')
p3<-ggplot()+
	geom_step(data=datd[[z]],aes(x=pos,y=snp),color='black')+
	xlab('')+
	xlim(0,size1)+
	ylim(0,max(datd[[z]]$snp,datdl[[z]]$snp))+
	geom_vline(xintercept=cont_pos[[z]],color='black',linetype=3,size=0.8)+
	geom_vline(xintercept=line$linesl,color='red',size=0.2)+
	geom_vline(xintercept=line$linesr,color='red',size=0.2)+
	geom_rect(data=line,inherit.aes=F,aes(ymin=-Inf,ymax=+Inf,xmin=linesr,xmax=linesl,group=grp),fill='red',alpha=0.2)+
	theme(plot.margin=unit(c(0,2.4,0,0),'cm'))+
	theme(legend.position='none')+
	theme(axis.text.y = element_text(angle = 90, hjust = 1))+
	ylab('SNP density')
p4<-ggplot()+
	geom_step(data=datxl[[z]],aes(x=snp,y=rec))+
	xlab('')+
	ylab('Recombination')+
	ylim(0,20)+
	xlim(0,size2)+
	theme(axis.text.y=element_text(angle=90))+
	geom_vline(xintercept=lineL$linesll,color='red',size=0.2)+
	geom_vline(xintercept=lineL$lineslr,color='red',size=0.2)+
	geom_rect(data=lineL,inherit.aes=F,aes(ymin=-Inf,ymax=+Inf,xmin=lineslr,xmax=linesll,group=grp),fill='red',alpha=0.2)+
	theme(plot.margin=unit(c(2.6,1,1,0),'cm'))+
	theme(legend.position='none')+
	theme(axis.text.y = element_text(angle = 90, hjust = 1))+
	coord_flip()
p5<-ggplot()+
	geom_step(data=datdl[[z]],aes(x=pos,y=snp),color='black')+
	xlab('')+
	xlim(0,size2)+
	ylim(0,max(datd[[z]]$snp,datdl[[z]]$snp))+
	theme(axis.text.y=element_text(angle=90))+
	geom_vline(xintercept=lineL$linesll,color='red',size=0.2)+
	geom_vline(xintercept=lineL$lineslr,color='red',size=0.2)+
	geom_rect(data=lineL,inherit.aes=F,aes(ymin=-Inf,ymax=+Inf,xmin=lineslr,xmax=linesll,group=grp),fill='red',alpha=0.2)+
	theme(plot.margin=unit(c(2.6,0,1,0),'cm'))+
	theme(legend.position='none')+
	ylab('SNP density')+
	theme(axis.text.y = element_text(angle = 90, hjust = 1))+
	coord_flip()
grid.arrange(p5,p4,p1,
			blankPlot,blankPlot,p2,
			blankPlot,blankPlot,p3,
			ncol=3,nrow=3,widths=c(1,1.5,3),heights=c(3,1.5,1))	
dev.off()

##### snp present in union of 3 colonies
type<-'colony'
file_name<-paste0(args[6],type,'.txt')
#file_name<-paste0('Mb_map1_',type,'.txt')
size_bin<-gsub('map','',unlist(strsplit(file_name,"_|.txt"))[2])
mbmap<-read.table(paste0(parent_folder,'/',file_name),header=T,sep=' ')
gmap<-read.table(paste0(parent_folder,'/','genetic_map_',type,'.txt'),header=T,sep=',')
chrmap<-read.table(paste0(parent_folder,'/','chromosome_map_',type,'.txt'),header=T,sep=',')
mbmapl<-read.table(paste0(parent_folderl,'/',file_name),header=T,sep=' ')
gmapl<-read.table(paste0(parent_folderl,'/','genetic_map_',type,'.txt'),header=T,sep=',')
chrmapl<-read.table(paste0(parent_folderl,'/','chromosome_map_',type,'.txt'),header=T,sep=',')

datx<-list()
datd<-list()
cont_pos<-list()
for(z in 1:nchr){
	datc<-list()
	datdc<-list()
	cont_posc<-list()
for(c in 1:nb_colony){
	chrn<-read.table(paste0(parent_folder,'/info_all',c,'1.txt'),header=T,sep=',')
	nind<-length(which(grepl('ind',colnames(chrn))==T))
	snp<-vector()
	chr<-vector()
	contig<-vector()
	vector<-vector()
	for(i in 1:nrow(chrn)){snp<-c(snp,chrn$posStart[i],chrn$posStop[i])
		chr<-c(chr,chrn$chr[i],chrn$chr[i])
		contig<-c(contig,chrn$contig[i],chrn$contig[i])
		vector<-c(vector,paste0(chrn[i,9:(9+nind-1)],collapse=''),paste0(chrn[i,9:(9+nind-1)],collapse=''))}
	chr<-data.frame(chromosome=chr,contig,snp,vector)
	chr$vector<-as.character(chr$vector)
	s<-subset(chr,chr$chromosome==z)
	rec<-vector()
	j<-1
	while(j<nrow(s)){
	vj<-s$vector[j]
	vk<-s$vector[j+1]
	r<-0
	for (k in 1:nchar(vj)){
	if(substr(vj,k,k)!=substr(vk,k,k)){r<-r+1}
	}
	rec<-c(rec,r)
	j=j+1}
datc[[c]]<-data.frame(chromosome=s$chromosome,contig=s$contig,snp=s$snp,rec=c(rec,0))
datc[[c]]$col<-c
c_pos<-vector()
for(i in 1:length(unique(datc[[c]]$contig))){
	s<-datc[[c]][datc[[c]]$contig==unique(datc[[c]]$contig)[i],'snp']
	c_pos<-c(c_pos,s[1])}
cont_posc[[c]]<-data.frame(contig=unique(datc[[c]]$contig),c_pos,col=c)
}
datx[[z]]<-do.call(rbind,datc)
datx[[z]]<-datx[[z]][order(datx[[z]]$snp),]
value='snp'
inputseq=datx[[z]]
starts<-c(seq(1,max(inputseq$snp),by=window1),max(inputseq$snp))
n<-length(starts)
snps<-numeric(n)
for(i in 1:(n-1)){
	if(i==(n-1)){chunk<-subset(inputseq[,paste0(value)],inputseq[,paste0(value)]>=starts[i] & inputseq[,paste0(value)]<=starts[i+1])
	}else{chunk<-subset(inputseq[,paste0(value)],inputseq[,paste0(value)]>=starts[i] & inputseq[,paste0(value)]<starts[i+1])} 
	if(any(is.na(chunk))){chunk<-chunk[!is.na(chunk)]}
    sum_<-length(chunk)
    snps[i]<-sum_
     }
datd[[z]]<-data.frame(pos=starts,snp=snps)
cont_pos[[z]]<-Reduce(function(x,y) merge(x,y,by='contig',all=T),cont_posc)
for(i in 1:nrow(cont_pos[[z]])){
	contig_start<-cont_pos[[z]][i,grepl('c_pos',colnames(cont_pos[[z]]))==T]
	cont_pos[[z]]$pos[i]<-min(contig_start[!is.na(contig_start)])
	}
}

datxl<-list()
datdl<-list()
for(z in 1:nchr){
	datc<-list()
	datdc<-list()
for(c in 1:nb_colony){
	chrn<-read.table(paste0(parent_folderl,'/info_all',c,'1.txt'),header=T,sep=',')
	nind<-length(which(grepl('ind',colnames(chrn))==T))
	snp<-vector()
	chr<-vector()
	vector<-vector()
	for(i in 1:nrow(chrn)){snp<-c(snp,chrn$posStart[i],chrn$posStop[i])
		chr<-c(chr,chrn$chr[i],chrn$chr[i])
		vector<-c(vector,paste0(chrn[i,9:(9+nind-1)],collapse=''),paste0(chrn[i,9:(9+nind-1)],collapse=''))}
	chr<-data.frame(chromosome=chr,snp,vector)
	chr$vector<-as.character(chr$vector)
	s<-subset(chr,chr$chromosome==z)
	rec<-vector()
	j<-1
	while(j<nrow(s)){
	vj<-s$vector[j]
	vk<-s$vector[j+1]
	r<-0
	for (k in 1:nchar(vj)){
	if(substr(vj,k,k)!=substr(vk,k,k)){r<-r+1}
	}
	rec<-c(rec,r)
	j=j+1}
datc[[c]]<-data.frame(chromosome=s$chromosome,snp=s$snp,rec=c(rec,0))
datc[[c]]$col<-c}
datxl[[z]]<-do.call(rbind,datc)
datxl[[z]]<-datxl[[z]][order(datxl[[z]]$snp),]
value='snp'
inputseq=datxl[[z]]
starts<-c(seq(1,max(inputseq$snp),by=window1),max(inputseq$snp))
n<-length(starts)
snps<-numeric(n)
for(i in 1:(n-1)){
	if(i==(n-1)){chunk<-subset(inputseq[,paste0(value)],inputseq[,paste0(value)]>=starts[i] & inputseq[,paste0(value)]<=starts[i+1])
	}else{chunk<-subset(inputseq[,paste0(value)],inputseq[,paste0(value)]>=starts[i] & inputseq[,paste0(value)]<starts[i+1])} 
	if(any(is.na(chunk))){chunk<-chunk[!is.na(chunk)]}
    sum_<-length(chunk)
    snps[i]<-sum_
     }
datdl[[z]]<-data.frame(pos=starts,snp=snps)
}

pdf('PacBiovsLiu_colony.pdf',width=20,height=20)
for(z in 1:nchr){
	print(z)
linesl<-list()
for(c in 1:nb_colony){
	d<-datxl[[z]][datxl[[z]]$col==c,]
	n<-which(d$rec>=5)
if(length(n)>0){l<-vector()
for(i in 1:length(n)){s<-d[n[i]:(n[i]+1),]
	l<-c(l,mean(s$snp))
	linesl[[c]]<-l}}else{linesl[[c]]<-(-100)}}
lines<-list()
for(c in 1:nb_colony){
	d<-datx[[z]][datx[[z]]$col==c,]
	n<-which(d$rec>=5)
if(length(n)>0){l<-vector()
for(i in 1:length(n)){s<-d[n[i]:(n[i]+1),]
	l<-c(l,mean(s$snp))
	lines[[c]]<-l}}else{lines[[c]]<-(-100)}}
lineslc<-unlist(linesl)
linesc<-unlist(lines)
chrT= read.table(paste("Chr",z,"_split1_chained.csv",sep=""),header = TRUE)
chromosome = paste("chr",z,sep="")
contigLimits1 = agpFile[agpFile$V1 == chromosome,]$V2
contigLimits2 = agpFile[agpFile$V1 == chromosome,]$V3
a<-seq(1:(nrow(chrT)/3))
grp<-vector()
for(j in 1:length(a)){grp<-c(grp,rep(a[j],3))}
chrT$grp<-grp
size1<-max(chrT$query[!is.na(chrT$query)])
size2<-max(chrT$target[!is.na(chrT$query)])
p1<-ggplot(chrT,aes(query,target,group=grp))+
	geom_point(shape='.')+
	geom_line(size=0.3)+
	xlab("AMelMel (position in bp)")+
	xlim(0,size1)+
	ylim(0,size2)+
	ylab("Amel4.5 (position in bp)")+
	geom_vline(xintercept=cont_pos[[z]]$pos,size=0.8,linetype=3)+
	geom_hline(yintercept=lineslc,size=0.3,color='red')+
	geom_vline(xintercept=linesc,size=0.3,color='red')+
	ggtitle(paste0("Chromosome ",z))+
	theme(plot.title = element_text(size = 20, face = "bold"))+
	theme(axis.text.y=element_text(angle=90,hjust=1))+
	theme(plot.margin=unit(c(2,2.5,1,0),'cm'))
p21<-ggplot()+
	geom_step(data=datx[[z]][datx[[z]]$col=='1',],aes(x=snp,y=rec),colour='darkblue')+
	geom_step(data=datd[[z]],aes(x=pos,y=snp),colour='blue')+
	xlab('')+
	xlim(0,size1)+
	geom_vline(xintercept=cont_pos[[z]]$pos,linetype=3,size=0.8)+
	geom_vline(xintercept=lines[[1]],size=0.3,color='red')+
	theme(plot.margin=unit(c(1,1,1,0),'cm'))+
	theme(legend.position='none')+
	scale_y_continuous(name='Recombination',sec.axis=(sec_axis(~.,name='SNP density')))+	
	theme(axis.title.y.right=element_text(color='blue'))
p22<-ggplot()+
	geom_step(data=datx[[z]][datx[[z]]$col=='2',],aes(x=snp,y=rec),colour='magenta')+
	xlab('')+
	xlim(0,size1)+
	ylim(0,10)+
	geom_vline(xintercept=cont_pos[[z]]$pos,linetype=3,size=0.8)+
	geom_vline(xintercept=lines[[2]],size=0.3,color='red')+
	ylab('Recombination')+	
	theme(plot.margin=unit(c(0,2.4,0.8,0),'cm'))+
	theme(legend.position='none')
p23<-ggplot()+
	geom_step(data=datx[[z]][datx[[z]]$col=='3',],aes(x=snp,y=rec),colour='darkgreen')+
	xlab('')+
	xlim(0,size1)+
	ylim(0,10)+
	geom_vline(xintercept=cont_pos[[z]]$pos,linetype=3,size=0.8)+
	geom_vline(xintercept=lines[[3]],size=0.3,color='red')+
	ylab('Recombination')+	
	theme(plot.margin=unit(c(0,2.4,0.8,0),'cm'))+
	theme(legend.position='none')
mbmapz<-mbmap[mbmap$chr==z,]
mbmapz<-rbind(mbmapz,c(z,mbmapz[nrow(mbmapz),'right'],size1,mbmapz[nrow(mbmapz),'m_cj'],mbmapz[nrow(mbmapz),'s_cj'],mbmapz[nrow(mbmapz),'q5_cj'],mbmapz[nrow(mbmapz),'q95_cj']))
p3<-ggplot(mbmapz)+
	geom_step(aes(left,m_cj))+
	geom_hline(aes(yintercept=22),color='blue',linetype='dashed')+
	xlab('')+
	xlim(0,size1)+
	ylim(0,160)+
	geom_vline(xintercept=cont_pos[[z]]$pos,size=0.8,linetype=3)+
	theme(plot.margin=unit(c(0,2,0.8,0),'cm'))+
	ylab('Recombination rate cm/Mb')
dat<-read.fasta(file=paste0("Fasta/AMelMelPacBio_chr",z,".fa"))
datseq <- dat[[1]]
datGC<-slidingwindow(window1, datseq)
colnames(datGC)<-c('pos','GC')
cov_mean<-list()
for(c in 1:nb_colony){
cov<-read.table(paste0('PacBio/samCovQueen',c,'_',z,'_sub.txt'))
colnames(cov)<-c('chr','pos','cov')
pos_mean<-vector()
a=1
b=window1
c_mean<-vector()
while(b<nrow(cov)){
c_mean<-c(c_mean,mean(cov$cov[a:b]))
pos_mean<-c(pos_mean,a)
a=b
b=b+window1
}
c_mean<-c(c_mean,mean(cov$cov[tail(pos_mean,1):a]),mean(cov$cov[a:nrow(cov)]))
pos_mean<-c(pos_mean,a,nrow(cov))
cov_mean[[c]]<-c_mean
}
cov2<-data.frame(pos_mean,do.call(cbind,cov_mean))
colnames(cov2)<-c('pos','c1','c2','c3')
dat_GCcov<-merge(datGC,cov2,by='pos')
p4<-ggplot(dat_GCcov)+
	geom_step(aes(x=pos,y=GC*100))+ 
	xlab('')+
	xlim(0,size1)+
	ylim(0,100)+
	geom_vline(xintercept=cont_pos[[z]]$pos,size=0.8,linetype=3)+
	theme(plot.margin=unit(c(0,1,0.8,0),'cm'))+
	geom_step(aes(x=pos,y=c1),linetype='longdash',color='blue')+ 
	geom_step(aes(x=pos,y=c2),linetype='dotdash',color='blue')+ 
	geom_step(aes(x=pos,y=c3),linetype='dotted',color='blue')+ 
	scale_y_continuous(name='GC %',sec.axis=(sec_axis(~.*1,name='Queen coverage')))+	
	theme(axis.title.y.right=element_text(color='blue'))
p51<-ggplot()+
	geom_step(data=datxl[[z]][datxl[[z]]$col=='1',],aes(x=snp,y=rec),colour='darkblue')+
	geom_step(data=datdl[[z]],aes(x=pos,y=snp),color='blue')+
	xlab('')+
	xlim(0,size2)+
	geom_vline(xintercept=linesl[[1]],color='red',size=0.3)+
	theme(axis.text.y=element_text(angle=90,hjust=1))+
	theme(plot.margin=unit(c(1.2,1,1,0),'cm'))+
	theme(legend.position='none')+
	scale_y_continuous(name='Recombination',sec.axis=(sec_axis(~.,name='SNP density')))+	
	theme(axis.title.x.top=element_text(color='blue'))+
	coord_flip()
p52<-ggplot()+
	geom_step(data=datxl[[z]][datxl[[z]]$col=='2',],aes(x=snp,y=rec),colour='magenta')+
	xlab('')+
	xlim(0,size2)+
	ylim(0,10)+
	geom_vline(xintercept=linesl[[2]],color='red',size=0.3)+
	theme(axis.text.y=element_text(angle=90,hjust=1))+
	theme(plot.margin=unit(c(2.5,1,1,0),'cm'))+
	theme(legend.position='none')+
	ylab('Recombination')+	
	coord_flip()
p53<-ggplot()+
	geom_step(data=datxl[[z]][datxl[[z]]$col=='3',],aes(x=snp,y=rec),colour='darkgreen')+
	xlab('')+
	xlim(0,size2)+
	ylim(0,10)+
	geom_vline(xintercept=linesl[[3]],color='red',size=0.3)+
	theme(axis.text.y=element_text(angle=90,hjust=1))+
	theme(plot.margin=unit(c(2.5,1,1,0),'cm'))+
	ylab('Recombination')+
	theme(legend.position='none')+
	coord_flip()
mbmaplz<-mbmapl[mbmapl$chr==z,]
p6<-ggplot(mbmaplz)+
	geom_step(aes(left,m_cj))+
	geom_hline(aes(yintercept=37),color='blue',linetype='dashed')+
	geom_vline(xintercept=linesl[[1]],color='red',size=0.3)+
	geom_vline(xintercept=linesl[[2]],color='red',size=0.3)+
	geom_vline(xintercept=linesl[[3]],color='red',size=0.3)+
	theme(axis.text.y=element_text(angle=90,hjust=1))+
	xlab('')+
	xlim(0,size2)+
	ylim(0,160)+
	theme(plot.margin=unit(c(2.5,1,1,0),'cm'))+
	ylab('Recombination rate cm/Mb')+
	coord_flip()
dat<-read.fasta(file=paste0("Fasta/ame_ref_Amel_4.5_chrLG",z,".fa"))
datseq <- dat[[1]]
datGC<-slidingwindow(window1, datseq)
p7<-ggplot(datGC)+
	geom_step(aes(x=starts,y=chunkGCs*100))+ 
	xlab('')+
	xlim(0,size2)+
	ylim(0,100)+
	theme(axis.text.y=element_text(angle=90,hjust=1))+
	geom_vline(xintercept=linesl[[1]],color='red',size=0.3)+
	geom_vline(xintercept=linesl[[2]],color='red',size=0.3)+
	geom_vline(xintercept=linesl[[3]],color='red',size=0.3)+
	theme(plot.margin=unit(c(2.5,1,1,0),'cm'))+
	ylab('GC %')+
	coord_flip()
grid.arrange(p7,p6,p53,p52,p51,p1,
			blankPlot,blankPlot,blankPlot,blankPlot,blankPlot,p21,
			blankPlot,blankPlot,blankPlot,blankPlot,blankPlot,p22,
			blankPlot,blankPlot,blankPlot,blankPlot,blankPlot,p23,
			blankPlot,blankPlot,blankPlot,blankPlot,blankPlot,p3,
			blankPlot,blankPlot,blankPlot,blankPlot,blankPlot,p4,
			ncol=6,nrow=6,widths=c(1.5,1.5,1,1,1.5,3),heights=c(3,1.5,1,1,1.5,1.5))	
}
dev.off()


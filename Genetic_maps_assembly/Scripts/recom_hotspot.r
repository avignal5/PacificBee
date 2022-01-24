################################################
### script to localise recombination hotspot
################################################
#!/usr/bin/env R
library('data.table')
library("seqinr")
library("colorout")
library("stringr")
library('ggplot2')
library('gridExtra')
options(scipen=999)

cont_info<-read.table('PacBio/combined_info_all.txt',sep=',',header=T)
rep<-read.table('repeat_pattern.txt',header=T,sep=' ')
genoP<-fread('PacBio/genotype_allphased.txt',sep=',',header=F,data.table=F,colClasses=c('numeric','numeric','numeric','character'))
colnames(genoP)<-c('chromosome','contig','snp','vector')
genoL<-fread('Liu/colony_allphased.txt',sep=',',header=F,data.table=F,colClasses=c('numeric','numeric','numeric','numeric','character'))
colnames(genoL)<-c('chromosome','contig','snp','ordre','vector')

rec<-list()
for(z in 1:max(genoP$chromosome)){	
s<-genoP[genoP$chromosome==z,]
reco<-vector()
j<-1
while(j<nrow(s)){
vj<-s$vector[j]
vk<-s$vector[j+1]
r<-0
for (k in 1:nchar(vj)){if(substr(vj,k,k)!=substr(vk,k,k)){r<-r+1}}
reco<-c(reco,r)
j=j+1}
rec[[z]]<-c(0,reco)}
recomb<-unlist(rec)
genoP$recom<-recomb

rec<-list()
for(z in 1:max(genoL$chromosome)){	
s<-genoL[genoL$chromosome==z,]
reco<-vector()
j<-1
while(j<nrow(s)){
vj<-s$vector[j]
vk<-s$vector[j+1]
r<-0
for (k in 1:nchar(vj)){if(substr(vj,k,k)!=substr(vk,k,k)){r<-r+1}}
reco<-c(reco,r)
j=j+1}
rec[[z]]<-c(0,reco)}
recomb<-unlist(rec)
genoL$recom<-recomb

n<-which(genoP$recom==max(genoP$recom))
repn<-paste0(rep('n',100),collapse='')

seq_cont<-list()
repl<-list()
fasta<-list()
posl<-list()
pdf('hotspot.pdf',width=10,height=10)
for(i in 1:length(n)){
	chri<-genoP[n[i],'chromosome']	
	conti<-genoP[n[i],'contig']
	cont_infoi<-cont_info[cont_info$chr_nb_consensus==chri,]
	cont_ordre<-which(cont_infoi$cont==conti)
	cont_size<-cont_infoi[cont_infoi$cont==conti,'size']
	cont_orient<-cont_infoi[cont_infoi$cont==conti,'orientation_final']
	genoi<-genoP[genoP$contig==conti,]
	fasta[[i]]<-read.fasta(file=paste0("Fasta/AMelMelPacBio_chr",chri,".fa"),as.string=T,set.attributes=F)
#	fastaP<-read.fasta(file=paste0("Fasta/AMelMelPacBio_chr",chri,".fa"),as.string=T,set.attributes=F)
#	fastaL<-read.fasta(file=paste0("Fasta/ame_ref_Amel_4.5_chrLG",chri,".fa"),as.string=T,set.attributes=F)	
	pos=gregexpr(repn,fasta[[i]]) 
	posi<-c(1,pos[[1]],nchar(fasta[[i]]))
	if(length(posi)-1!=nrow(cont_infoi)){
		if(cont_ordre==nrow(cont_infoi)){
				pos_contm=posi[length(posi)-1]+100
				pos_contM=posi[length(posi)]}
	}else{pos_contm=posi[cont_ordre]+100
				pos_contM=posi[cont_ordre+1]
					}
	if(cont_orient=='-'){genoi$seq<-cont_size-genoi$snp+pos_contm
		rec<-list()
		reco<-vector()	
		j<-1
		while(j<nrow(genoi)){
		vj<-genoi$vector[j]
		vk<-genoi$vector[j+1]
		r<-0
		for (k in 1:nchar(vj)){if(substr(vj,k,k)!=substr(vk,k,k)){r<-r+1}}
		reco<-c(reco,r)
		j=j+1}
		genoi$recom<-c(0,reco)
	}else{genoi$seq<-genoi$snp+pos_contm
		rec<-list()
		reco<-vector()	
		j<-1
		while(j<nrow(genoi)){
		vj<-genoi$vector[j]
		vk<-genoi$vector[j+1]
		r<-0
		for (k in 1:nchar(vj)){if(substr(vj,k,k)!=substr(vk,k,k)){r<-r+1}}
		reco<-c(reco,r)
		j=j+1}
		genoi$recom<-c(0,reco)}
	geno_plot<-genoi[,c('recom','seq')]
	x<-which(geno_plot$recom!=0)
	geno_plot$recom[x-1]<-geno_plot$recom[x]
	repi<-subset(rep,rep$chr==chri & rep$start>pos_contm & rep$stop<=pos_contM)
	repz<-list()
	for(j in 1:nrow(repi)){repz[[j]]<-data.frame(pos=c(repi$start[j],repi$stop[j]),nb_copies=rep(repi$nb_copies[j],2))}
	repx<-do.call(rbind,repz)
	l<-geno_plot$seq[which(geno_plot$recom==max(geno_plot$recom))]
	posl[[i]]<-l
	line<-data.frame(min=l[1],max=l[2])
	p1<-ggplot()+
	geom_step(data=geno_plot,aes(x=seq,y=recom))+
	xlab('')+
	ylim(0,5)+
	ylab('Recombination')+
	geom_vline(xintercept=l,color='red',size=0.1)+
	geom_rect(data=line,inherit.aes=F,aes(ymin=-Inf,ymax=+Inf,xmin=min,xmax=max),fill='red',alpha=0.2)+
	theme(legend.position='none')+
	theme(plot.margin=unit(c(0.2,0.5,0.2,1.1),'cm'))
	p2<-ggplot()+
	geom_step(data=repx,aes(x=pos,y=nb_copies))+
	xlab('')+
	ylab('# copies')+
	geom_vline(xintercept=l,color='red',size=0.1)+
	geom_rect(data=line,inherit.aes=F,aes(ymin=-Inf,ymax=+Inf,xmin=min,xmax=max),fill='red',alpha=0.2)+	
	theme(legend.position='none')+
	theme(plot.margin=unit(c(0.2,0.5,0.2,0.5),'cm'))
grid.arrange(p1,p2,
			ncol=1,nrow=2,widths=2,heights=c(1,1))	
seq_cont[[i]]<-substr(fasta[[i]],l[1],l[2])
repl[[i]]<-subset(rep,rep$chr==chri & rep$start>l[1] & rep$stop<=l[2])
write.table(seq_cont[[i]],paste0('Fasta/hotsplot',i,'.txt'),col.names=F,quote=F,row.names=F)
}
dev.off()


################################################
### script to build chromosome_map and genetic_map necessary to calculate recombinaison map
################################################
#!/usr/bin/env R
args<-commandArgs(TRUE)
nb_colony<-as.numeric(args[1])
pattern<-args[2] 
pattern_geno<-args[4]
library('data.table')
library('reshape2')
first.changes <- function(d) {
  p <- cumsum(rle(d)$lengths) + 1
  p[-length(p)]
}

data<-read.table(args[3],sep=',',header=TRUE)
b<-as.numeric(as.character(data$chr_nb_consensus))
nchr<-max(b[!is.na(b)])

options(warn=2)
final_c<-list()
datc<-list()
for(c in 1:nb_colony){
geno<-fread(paste0(pattern_geno,c,'.txt'),sep=' ',header=T,data.table=F)
geno_sub<-colsplit(geno[,1],'_',c('contig','marker'))
geno_sub$contig<-gsub('LG','',geno_sub$contig)
chr<-read.table(paste(pattern,c,'/chromosome/all_chr_new2.txt',sep=''),header=TRUE,sep=',',colClasses=c('character','character','character','character','character','character','character','character','character','character','character','character','character','character','character','character'))
chr_i<-data.frame(chr=chr$chromosome,contig=chr$contig,pos_m=chr$pos_m,pos_M=chr$pos_M,colsplit(chr$vector1,'',paste0('ind',seq(1:nchar(chr$vector1[1])))))
nind<-ncol(chr_i)-4
final_z<-list()
datz<-list()
for(z in 1:nchr){
print(paste0('chromosome: ',z))
cont<-data[data$chr_nb_consensus==z,'cont']
orient<-data[data$chr_nb_consensus==z,'orientation_final']
size<-vector()
for(a in 1:length(cont)){
	sg<-data[data$chr_nb_consensus==z & data$cont==cont[a],'size']
	if(length(geno_sub$marker[geno_sub$contig==cont[a]])>0){
	sc<-max(geno_sub$marker[geno_sub$contig==cont[a]])
	if(sg>sc){size<-c(size,sg)}else{size<-c(size,sc)}}
	else{size<-c(size,sg)}}
for(p in 1:ncol(chr_i)){chr_i[,p]<-as.numeric(as.character(chr_i[,p]))}
final_x<-list()
datx<-list()
for(x in 1:length(cont)){
	X=x
	g3<-subset(geno_sub,geno_sub$contig==cont[X])
	if(nrow(g3)>0){
	g3$chr<-z
	g3$colony<-c
	if(orient[x]=='+'){g3<-g3[order(g3$marker),]}else if (orient[x]=='-'){g3<-g3[order(-g3$marker),]}
	si<-sum(size[0:(x-1)])
	if(orient[x]=='+'){g4<-data.frame(m=seq(1,size[x]),o=seq((si+1),(si+length(seq(1,size[x])))))}else if (orient[x]=='-'){g4<-data.frame(m=rev(seq(1,size[x])),o=seq((si+1),(si+length(seq(1,size[x])))))}			
	g3$seq<-g4[g4$m%in%g3$marker,'o']
	final_x[[x]]<-g3}else{x=x+1}
	chr_n<-subset(chr_i,chr_i$chr==z&chr_i$contig==cont[X])
	if(nrow(chr_n)>1){
	ind<-vector()
	marker_g<-vector()
	marker_d<-vector()
	pos_g<-vector()
	pos_d<-vector()
	for(i in 1:nind){
	w<-first.changes(chr_n[,paste0('ind',i)])
	if(orient[x]=='+'){
	mM<-chr_n[w,'pos_m']
	mm<-chr_n[(w-1),'pos_M']
	pM<-g3[g3$marker%in%mM,'seq']
	pm<-g3[g3$marker%in%mm,'seq']
	}else if (orient[x]=='-'){
	mM<-chr_n[w,'pos_M']
	mm<-chr_n[(w-1),'pos_m']
	pM<-g3[g3$marker%in%mM,'seq']
	pm<-g3[g3$marker%in%mm,'seq']
	}
	ind<-c(ind,rep(i,length(w)))
	marker_g<-c(marker_g,mm)
	marker_d<-c(marker_d,mM)
	pos_g<-c(pos_g,pm)
	pos_d<-c(pos_d,pM)}
	datx[[x]]<-data.frame(col=c,chr=z,contig=cont[x],ind,marker_g,marker_d,pos_g,pos_d)}else{x=x+1}
}
datz[[z]]<-do.call(rbind,datx)
final_z[[z]]<-do.call(rbind,final_x)	
}
datc[[c]]<-do.call(rbind,datz)
final_c[[c]]<-do.call(rbind,final_z)
}
dat_final<-do.call(rbind,datc)
final_final<-do.call(rbind,final_c)
forder<-final_final[order(final_final$chr,final_final$seq),]
forder$contig<-NULL
forder$colony<-NULL
final_write<-forder[!duplicated(forder),]
final_write<-final_write[,c(2,1,3)]
colnames(final_write)<-c('chromosome','marker','position')
write.table(final_write,'chromosome_map.txt',sep=',',col.names=T,row.names=F,quote=F)
write.table(dat_final,'genetic_map.txt',sep=',',col.names=T,row.names=F,quote=F)

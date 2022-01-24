################################################
### script to summarise chromosome and contig information after ordering of bins
################################################
#!/usr/bin/env R
args<-commandArgs(TRUE)
nb_colony<-as.numeric(args[1])
pattern<-args[2]

cover<-list()
marker<-list()
for(c in 1:nb_colony){
chr<-read.table(paste0(pattern,c,'/chromosome/all_chromosome.txt'),sep=',',header=T,colClasses=c('character','character','character','character','numeric','numeric','numeric','numeric','character','numeric','numeric','numeric','numeric','numeric','numeric'),stringsAsFactor=F)
chr[is.na(chr)]<-0
nb_marker<-sum(as.numeric(chr$nb_markers_all))
coverage<-sum(as.numeric(chr$size_all))
print(paste0('contig placed ',length(unique(chr$contig))))
nb_m<-vector()
cov<-vector()
for(i in 1:max(as.numeric(chr$chromosome))){
	chri<-subset(chr,chr$chromosome==i)
	nb_m<-c(nb_m,sum(chri$nb_markers_all))
	cov<-c(cov,sum(chri$size_all))
	}
	
cover[[c]]<-c(coverage,cov)
marker[[c]]<-c(nb_marker,nb_m)
}

coverm<-list()
markerm<-list()
cont<-vector()
for(c in 1:nb_colony){
chr<-read.table(paste0(pattern,c,'/chromosome/all_chr_new.txt'),sep=',',header=T,colClasses=c('character','character','character','character','character','character','character','character','character','character','character','character','character','character','character','character'))
chr[is.na(chr)]<-0
cont<-c(cont,length(unique(as.numeric(chr$contig))))
nb_marker<-sum(as.numeric(chr$nb_markers_all))
coverage<-sum(as.numeric(chr$size_all))
chrs<-subset(chr,chr$chromosome%in%seq(1,16))
print(paste0('contig with disagreeing position ',length(unique(chr$contig))-length(unique(chrs$contig))))
print(paste0('contig possible to place ',length(unique(chrs$contig))))
nb_m<-vector()
cov<-vector()
for(i in 1:max(as.numeric(chrs$chromosome))){
	chri<-subset(chrs,chrs$chromosome==i)
	nb_m<-c(nb_m,sum(as.numeric(chri$nb_markers_all)))
	cov<-c(cov,sum(as.numeric(chri$size_all)))
	}
coverm[[c]]<-c(coverage,cov)
markerm[[c]]<-c(nb_marker,nb_m)
}

cover1<-do.call(cbind,cover)
cover2<-do.call(cbind,coverm)
marker1<-do.call(cbind,marker)
marker2<-do.call(cbind,markerm)

coverage<-data.frame(cover1,cover2)
markers<-data.frame(marker1,marker2)

coverage
markers

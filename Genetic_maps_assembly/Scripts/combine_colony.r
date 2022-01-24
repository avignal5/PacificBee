################################################
### script to combine phased vectors of 3 colonies
################################################
#!/usr/bin/env R
args<-commandArgs(TRUE)
nb_colony<-as.numeric(args[1])
genoname<-args[2]
pattern<-args[3]
info_cont<-read.table(args[4],sep=',',header=T)
info_cont<-info_cont[info_cont$geno=='y',]
info_co<-info_cont[order(info_cont$chr_a,info_cont$pos),]

options(warn=2)
lcol<-list()
for(c in 1:nb_colony){
lmarkers<-read.table(paste0(genoname,'_colony',c,'.txt'),header=T,sep=' ')[,1]
lmarkers<-as.character(lmarkers)
d<-matrix(ncol=2,nrow=length(lmarkers))
for(i in 1:length(lmarkers)){d[i,1]<-unlist(strsplit(lmarkers[i],'_'))[1]
	d[i,2]<-unlist(strsplit(lmarkers[i],'_'))[2]}
d<-as.data.frame(d)
if(length(unique(d[,1]))==16){colnames(d)<-c('chr','marker')
d$chr<-gsub('LG','',d$chr)
d$contig<-d$chr}else{colnames(d)<-c('contig','marker')}
chr<-read.table(paste(pattern,c,'/chromosome/all_chr_new2.txt',sep=''),header=T,sep=',',colClasses=c('character','character','character','character','numeric','numeric','numeric','numeric','character','numeric','numeric','numeric','numeric','numeric','numeric'))
lz<-list()
for(z in 1:max(chr$chromosome)){
	cz<-chr[chr$chromosome==z,]
	lf<-list()
	for(f in 1:length(unique(cz$contig))){
		cc<-cz[cz$contig==unique(cz$contig)[f],]
		lc<-d[d$contig==unique(cz$contig)[f],]
		lc$marker<-as.numeric(as.character(lc$marker))
		cc$pos_m<-as.numeric(as.character(cc$pos_m))
		cc$pos_M<-as.numeric(as.character(cc$pos_M))
		if(nrow(lc)>0){
		for(l in 1:nrow(lc)){
			lc$vector[l]<-subset(cc$vector1,lc$marker[l]<=cc$pos_M & lc$marker[l]>=cc$pos_m)
			lc$order_marker[l]<-l
			lc$chr<-z
			}}
		lf[[f]]<-lc
		}
	lz[[z]]<-do.call(rbind,lf)
	}
lcol[[c]]<-do.call(rbind,lz)
lcol[[c]]<-lcol[[c]][order(lcol[[c]]$contig),]
}
data<-Reduce(function(x, y) merge(x, y, by=c('chr','contig','marker','order_marker')), lcol)
data$vect<-apply(data[,grepl('vect',names(data))],1,paste,collapse='')
data$vector.x<-NULL
data$vector.y<-NULL
data$vector<-NULL

if(length(unique(data$chr))==length(unique(data$contig))){
data$chr<-as.numeric(as.character(data$chr))
data$contig<-as.numeric(as.character(data$contig))
data$marker<-as.numeric(as.character(data$marker))
geno_final<-data[order(data$chr,data$contig,data$marker),]
}else{
geno<-list()
for(z in 1:max(chr$chromosome)){
	cont<-info_cont$cont[info_cont$chr_nb_consensus%in%c(z,paste0('Un',z))]
	orient<-info_cont$orientation_final[info_cont$chr_nb_consensus%in%c(z,paste0('Un',z))]
	ordre<-info_cont$ordre_consensus[info_cont$chr_nb_consensus%in%c(z,paste0('Un',z))]
	info<-data.frame(cont,orient,ordre)
	n<-which(grepl('0_1_',info$ordre)==T)
	if(length(n)>0){info<-rbind(info[n,],info[-n,])	
	info$ordre<-seq(1:nrow(info))}
	info$ordre<-as.numeric(as.character(info$ordre))
	lcont<-unique(data$contig[data$chr==z])
	gf<-list()
	for(f in 1:nrow(info)){
		g<-subset(data,data$contig==info$cont[f])
		g$marker<-as.numeric(as.character(g$marker))
		g<-g[order(g$marker),]
		if(info$orient[f]=='-'){g<-g[rev(rownames(g)), ]}
		gf[[f]]<-g
		}
	geno[[z]]<-do.call(rbind,gf)
	}
geno_final<-do.call(rbind,geno)
geno_final$order_marker<-NULL}

write.table(geno_final,paste0(genoname,'phased.txt'),sep=',',col.names=F,row.names=F,quote=F)



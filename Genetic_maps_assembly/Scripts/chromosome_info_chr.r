################################################
### script to summarise chromosome information after ordering of bins
################################################
#!/usr/bin/env Rargs<-commandArgs(TRUE)
nb_colony<-as.numeric(args[1])
pattern<-args[2]

cover<-list()
marker<-list()
for(c in 1:nb_colony){
chr<-read.table(paste0(pattern,c,'/chromosome/all_chromosome.txt'),sep=';',header=T,colClasses=c('character','character','character','character','numeric','numeric','numeric','character','numeric','numeric'),stringsAsFactor=F)
chr[is.na(chr)]<-0
chr$X<-NULL
chr$chromosome<-NULL
colnames(chr)<-c('marker1','vector1','vector2','nb_markers_all','chromosome','size_all','list','recombination')
for(i in 1:nrow(chr)){
	if(grepl('_',chr$marker1[i])==T){chr$marker2[i]<-unlist(strsplit(chr$marker1[i],'_'))[1]}else{chr$marker2[i]<-paste0(chr$marker1[i],'_')}
	chr$contig[i]<-chr$chromosome[i]
	chr$bin[i]<-unlist(strsplit(chr$marker1[i],'_'))[1]
	if(grepl('/',chr$list[i])==T){
		l1<-unlist(strsplit(chr$list[i],','))
		l<-vector()
		for(n in 1:length(l1)){l2<-strsplit(l1[n],'/')
			l<-c(l,unlist(l2)[1])}
			l<-as.numeric(l)
	}else{l<-as.numeric(unlist(strsplit(chr$list[i],',')))}
	chr$pos_m[i]<-min(l)
	chr$pos_M[i]<-max(l)
	dist<-vector()
	for(j in 1:length(l)-1){dist<-c(dist,l[j+1]-l[j])}
	chr$dist_mean[i]<-mean(dist)
	chr$dist_median[i]<-median(dist)
	chr$dist_var[i]<-var(dist)
	chr$density_marker[i]<-(max(l)-min(l))/length(l)
	}
chr$list<-NULL
chr<-chr[c("marker1", "vector1", "marker2","vector2","contig","bin","nb_markers_all","size_all","dist_mean","dist_median","dist_var","density_marker","recombination","chromosome","pos_m","pos_M")]
nb_marker<-sum(as.numeric(chr$nb_markers_all))
coverage<-sum(as.numeric(chr$size_all))
nb_m<-vector()
cov<-vector()
for(i in 1:max(as.numeric(chr$chromosome))){
	chri<-subset(chr,chr$chromosome==i)
	nb_m<-c(nb_m,sum(chri$nb_markers_all))
	cov<-c(cov,sum(chri$size_all))
	}
cover[[c]]<-c(coverage,cov)
marker[[c]]<-c(nb_marker,nb_m)
chr[is.na(chr)]<-0
write.table(chr,paste0(pattern,c,'/chromosome/all_chr_new.txt'),sep=',',col.names=T,row.names=F,quote=F)
}


cover1<-do.call(cbind,cover)
marker1<-do.call(cbind,marker)

cover1
marker1

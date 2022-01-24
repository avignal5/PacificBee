################################################
### script to phase bins for each contig and order them along genome (chromosome by chromosome and all chromosome together)
################################################
#!/usr/bin/env R
args<-commandArgs(TRUE)
nb_colony<-as.numeric(args[1])
pattern<- args[3]

dat<-list()
datdist2<-list()
for(c in 1:nb_colony){
filenames<-list.files(paste(pattern,c,sep=''),pattern='*.txt',full.names=T)
nb_m<-list()
p<-list()
s<-list()
group_id<-list()
V<-list()
V_<-list()
M<-list()
dMean<-list()
dMed<-list()
dVar<-list()
datdist1<-list()
bins2<-list()
for (i in 1:length(filenames)){
	bin<-read.table(filenames[i],header=T,sep=';',colClasses=c('numeric','numeric','character','character','numeric','character','numeric','character'))
	nb_m[[i]]<-bin$nb_markers
	pos<-0
	size<-0
	M_nb<-0
	dmean<-0
	dmed<-0
	dvar<-0
	bins0<-vector()
	datdist<-list()
	for (j in 1:nrow(bin)){
		bin1=gsub('[[]','',bin[j,8])
		bin2=gsub('[]]','',bin1)
		b0=unlist(strsplit(bin2,','))
		b<-vector()
		for(k in 1:length(b0)){
			b<-c(b,as.numeric(unlist(strsplit(b0[k],'/'))[1]))
			}	
		pos[j]=paste0(min(b),'_',max(b))
		M_nb[j]<-length(b)
		d<-vector()
		for(m in 1:length(b-1)){d<-c(d,b[m+1]-b[m])}
		datdist[[j]]<-data.frame(contig=rep(bin$group[j],length(d)),bin=rep(bin$bin[j],length(d)),distance=d)
		a<-datdist[[j]]$distance[!is.na(datdist[[j]]$distance)]
		dmean[j]<-mean(a)
		dmed[j]<-median(a)
		dvar[j]<-var(a)
		bins0<-c(bins0,bin[j,2])
		}
	bins2[[i]]<-bins0
	datdist1[[i]]<-do.call(rbind,datdist)
	p[[i]]<-unlist(pos)
	s[[i]]<-bin$physical_length
	group_id[[i]]<-bin$group
	V[[i]]<-bin$vector0
	V_[[i]]<-bin$vector1
	M[[i]]<-unlist(M_nb)
	dMean[[i]]<-unlist(dmean)
	dMed[[i]]<-unlist(dmed)
	dVar[[i]]<-unlist(dvar)	
	}
DMean<-unlist(dMean)
DMed<-unlist(dMed)
DVar<-unlist(dVar)
group<-unlist(group_id)
nb_markers_all<-unlist(nb_m)
size_all<-unlist(s)
position_all<-unlist(p)
v<-unlist(V)
v_<-unlist(V_)
m_id<-unlist(M)
bin3<-unlist(bins2)
datdist2[[c]]<-do.call(rbind,datdist1)
datdist2[[c]]$contig<-as.numeric(as.character(datdist2[[c]]$contig))
datdist2[[c]]$bin<-as.numeric(as.character(datdist2[[c]]$bin))
datdist2[[c]]$distance<-as.numeric(as.character(datdist2[[c]]$distance))
dat[[c]]<-data.frame(group,bin=bin3,marker_id_carta=m_id,nb_markers_all,size_all,position_all,vector=v,vector_=v_,dist_mean=DMean,dist_median=DMed,dist_var=DVar)
write.table(dat[[c]],paste('bin_description_colony',c,'.txt',sep=''),col.names=T,row.names=F,quote=F,sep=',')
}

datacontig<-read.table(args[2],sep=',',header=T)

if(identical(datacontig$chr,as.character(datacontig$cont))==FALSE & max(as.numeric(as.vector(datacontig$cont[grepl('Un',datacontig$cont)==FALSE])))!=max(as.numeric(as.vector(datacontig$chr[grepl('Un',datacontig$chr)==FALSE])))){
datadistance<-do.call(rbind,datdist2)
datadistance$col<-c(rep('1',nrow(datdist2[[1]])),rep('2',nrow(datdist2[[2]])),rep('3',nrow(datdist2[[3]])))
col<-c(rep(1,nrow(dat[[1]])),rep(2,nrow(dat[[2]])),rep(3,nrow(dat[[3]])))

databins<-rbind(dat[[1]],dat[[2]],dat[[3]])
databins$col<-col
colnames(databins)[1]<-'contig'
databins$density_markers<-databins$nb_markers_all/databins$size_all
databins$density_markers[!is.finite(databins$density_markers)]<-0
databins$dist_mean[is.na(databins$dist_mean)]<-0
databins$dist_median[is.na(databins$dist_median)]<-0
databins$dist_var[is.na(databins$dist_var)]<-0

#datacontig<-read.table(args[2],sep=',',header=T)
idx<-datacontig[grepl('geno',colnames(datacontig))]
idx<-idx[rowSums(is.na(idx)) != ncol(idx), ]
datacontig<-subset(datacontig,rownames(datacontig)%in%rownames(idx))

nb_markers1<-vector()
nb_bin1<-vector()
size_markers1<-vector()
nb_markers2<-vector()
nb_bin2<-vector()
size_markers2<-vector()
nb_markers3<-vector()
size_markers3<-vector()
nb_bin3<-vector()
for(i in 1:nrow(datacontig)){
	s1<-subset(databins,databins$contig==datacontig$cont[i] & databins$col==1)
	s2<-subset(databins,databins$contig==datacontig$cont[i] & databins$col==2)
	s3<-subset(databins,databins$contig==datacontig$cont[i] & databins$col==3)
	nb_markers1<-c(nb_markers1,sum(s1$nb_markers_all))
	size_markers1<-c(size_markers1,sum(s1$size_all))
	nb_bin1<-c(nb_bin1,nrow(s1))
	nb_markers2<-c(nb_markers2,sum(s2$nb_markers_all))
	size_markers2<-c(size_markers2,sum(s2$size_all))
	nb_bin2<-c(nb_bin2,nrow(s2))	
	nb_markers3<-c(nb_markers3,sum(s3$nb_markers_all))
	size_markers3<-c(size_markers3,sum(s3$size_all))
	nb_bin3<-c(nb_bin3,nrow(s3))
	}
datacontig$nb_marker1<-nb_markers1
datacontig$size_markers1<-size_markers1
datacontig$nb_bin1<-nb_bin1
datacontig$nb_marker2<-nb_markers2
datacontig$size_markers2<-size_markers2
datacontig$nb_bin2<-nb_bin2
datacontig$nb_marker3<-nb_markers3
datacontig$size_markers3<-size_markers3
datacontig$nb_bin3<-nb_bin3

l1<-lm(datacontig$nb_bin1~datacontig$size_markers1)
tr1<-(l1$coef[2]/nchar(as.character(dat[[1]]$vector[1])))*10^6*10^2
l2<-lm(datacontig$nb_bin2~datacontig$size_markers2)
tr2<-(l2$coef[2]/nchar(as.character(dat[[2]]$vector[1])))*10^6*10^2
l3<-lm(datacontig$nb_bin3~datacontig$size_markers3)
tr3<-(l3$coef[2]/nchar(as.character(dat[[3]]$vector[1])))*10^6*10^2

print('colony 1')
print(paste('number of bins:',nrow(databins[databins$col==1,])),sep=' ')
print(paste('number of contigs:',length(unique(databins[databins$col==1,]$contig))),sep=' ')
print(paste('number of markers:',sum(databins[databins$col==1,]$nb_markers_all)))
print('colony 2')
print(paste('number of bins:',nrow(databins[databins$col==2,])),sep=' ')
print(paste('number of contigs:',length(unique(databins[databins$col==2,]$contig))),sep=' ')
print(paste('number of markers:',sum(databins[databins$col==2,]$nb_markers_all)))
print('colony 3')
print(paste('number of bins:',nrow(databins[databins$col==3,])),sep=' ')
print(paste('number of contigs:',length(unique(databins[databins$col==3,]$contig))),sep=' ')
print(paste('number of markers:',sum(databins[databins$col==3,]$nb_markers_all)))

pdf('bin_description.pdf',width=10,height=10)
par(mfrow=c(1,3))
plot(datacontig$size_markers1,datacontig$nb_bin1)
abline(l1,col='orange')
legend('topleft',legend=paste(round(tr1,digit=2),'cM/Mb',sep=' '))
plot(datacontig$size_markers2,datacontig$nb_bin2)
abline(l2,col='orange')
legend('topleft',tr2,legend=paste(round(tr2,digit=2),'cM/Mb',sep=' '))
plot(datacontig$size_markers3,datacontig$nb_bin3)
abline(l3,col='orange')
legend('topleft',tr3,legend=paste(round(tr3,digit=2),'cM/Mb',sep=' '))
library(ggplot2)
datadistance[is.na(datadistance)]<-0
datadistance$contig<-as.numeric(datadistance$contig)
qplot(contig,distance,data=datadistance[datadistance$col==1,],color=bin)
qplot(contig,distance,data=datadistance[datadistance$col==2,],color=bin)
qplot(contig,distance,data=datadistance[datadistance$col==3,],color=bin)
par(mfrow=c(3,3))
plot(databins[databins$col==1,]$size_all,databins[databins$col==1,]$nb_markers_all,ylim=c(0,2500))
plot(databins[databins$col==2,]$size_all,databins[databins$col==2,]$nb_markers_all,ylim=c(0,2500))
plot(databins[databins$col==3,]$size_all,databins[databins$col==3,]$nb_markers_all,ylim=c(0,2500))
plot(databins[databins$size_all<10000&databins$col==1,]$size_all,databins[databins$size_all<10000&databins$col==1,]$nb_markers_all,ylim=c(0,50))
plot(databins[databins$size_all<10000&databins$col==2,]$size_all,databins[databins$size_all<10000&databins$col==2,]$nb_markers_all,ylim=c(0,50))
plot(databins[databins$size_all<10000&databins$col==3,]$size_all,databins[databins$size_all<10000&databins$col==3,]$nb_markers_all,ylim=c(0,50))
plot(databins[databins$size_all<2000&databins$col==1,]$size_all,databins[databins$size_all<2000&databins$col==1,]$nb_markers_all,ylim=c(0,50))
plot(databins[databins$size_all<2000&databins$col==2,]$size_all,databins[databins$size_all<2000&databins$col==2,]$nb_markers_all,ylim=c(0,50))
plot(databins[databins$size_all<2000&databins$col==3,]$size_all,databins[databins$size_all<2000&databins$col==3,]$nb_markers_all,ylim=c(0,50))
dev.off()

hd<-list()
for(c in 1:nb_colony){
	print(paste0('colony ',c))
for(z in 1:16){
chr<-subset(datacontig$cont,datacontig$chr_a==z & datacontig$geno=='y')
nb_cont<-length(chr)
print(paste0('chromosome ',z))
print(paste0('number of contig on chromosome: ',length(chr)))
vect_name<-list()
vect_name_<-list()
vect_list<-list()
vect_list_<-list()
data_vect<-list()
dn<-list()
du<-list()
dun<-list()
contig_col<-list()
for(j in 1:length(chr)){
	contig_col[[j]]<-subset(databins,databins$contig==chr[j]&databins$col==c)
	if(nrow(contig_col[[j]])>1){
	vect_list[[j]]<-vector()
	vect_list_[[j]]<-vector()
	vect_name[[j]]<-vector()
	vect_name_[[j]]<-vector()
	v=1
	vect<-contig_col[[j]]$vector[v]
	vect_<-contig_col[[j]]$vector_[v]
	vect_name[[j]]<-paste0(contig_col[[j]]$contig[v],'/',contig_col[[j]]$bin[v])
	vect_name_[[j]]<-paste0(contig_col[[j]]$contig[v],'/',contig_col[[j]]$bin[v],'_')
	while(v<nrow(contig_col[[j]])){
		vect_list[[j]]<-c(vect_list[[j]],as.character(vect))
		vect_list_[[j]]<-c(vect_list_[[j]],as.character(vect_))
		v=v+1
		vect_t1<-contig_col[[j]]$vector[v]
		vect_t2<-contig_col[[j]]$vector_[v]
		vi<-as.character(vect)
		vj<-as.character(vect_t1)
		vk<-as.character(vect_t2)		
		r1<-0
		r2<-0
		for (k in 1:nchar(vi)){
			if(substr(vi,k,k)!=substr(vj,k,k)){r1<-r1+1}
			if(substr(vi,k,k)!=substr(vk,k,k)){r2<-r2+1}
		}
		if(r1<r2){vect_name[[j]]<-c(vect_name[[j]],paste0(contig_col[[j]]$contig[v],'/',contig_col[[j]]$bin[v]))
			vect_name_[[j]]<-c(vect_name_[[j]],paste0(contig_col[[j]]$contig[v],'/',contig_col[[j]]$bin[v],'_'))
			vect<-vect_t1
			vect_<-vect_t2
			}else if(r1>r2){vect_name[[j]]<-c(vect_name[[j]],paste0(contig_col[[j]]$contig[v],'/',contig_col[[j]]$bin[v],'_'))
			vect_name_[[j]]<-c(vect_name_[[j]],paste0(contig_col[[j]]$contig[v],'/',contig_col[[j]]$bin[v]))
			vect<-vect_t2
			vect_<-vect_t1}
		}
	if(v==nrow(contig_col[[j]])){
		if(grepl('_',tail(vect_name[[j]],1))==TRUE){vect_list[[j]]<-c(vect_list[[j]],as.character(contig_col[[j]][nrow(contig_col[[j]]),'vector_']))
		vect_list_[[j]]<-c(vect_list_[[j]],as.character(contig_col[[j]][nrow(contig_col[[j]]),'vector']))}else{vect_list_[[j]]<-c(vect_list_[[j]],as.character(contig_col[[j]][nrow(contig_col[[j]]),'vector_']))
		vect_list[[j]]<-c(vect_list[[j]],as.character(contig_col[[j]][nrow(contig_col[[j]]),'vector']))}
		}
	d1a<-c(vect_name[[j]][1],vect_list[[j]][1])
	d1b<-c(vect_name_[[j]][1],vect_list_[[j]][1])	
	dna<-c(tail(vect_name[[j]],1),tail(vect_list[[j]],1))
	dnb<-c(tail(vect_name_[[j]],1),tail(vect_list_[[j]],1))
	a<-rbind(d1a,dna)
	b<-rbind(d1b,dnb)
	dn[[j]]<-data.frame(a,b)
	colnames(dn[[j]])<-c('marker1','vector1','marker2','vector2')
	if(dn[[j]]$vector1[1]==dn[[j]]$vector1[2]){dun[[j]]<-dn[[j]]
	dn[[j]]<-NULL}
	}else{
	vect_list[[j]]<-as.character(contig_col[[j]]$vector)
	vect_list_[[j]]<-as.character(contig_col[[j]]$vector_)
	vect_name[[j]]<-paste0(contig_col[[j]]$contig,'/',contig_col[[j]]$bin)
	vect_name_[[j]]<-paste0(contig_col[[j]]$contig,'/',contig_col[[j]]$bin,'_')
	d1a<-c(vect_name[[j]][1],vect_list[[j]][1])
	d1b<-c(vect_name_[[j]][1],vect_list_[[j]][1])	
	a<-rbind(d1a,d1a)
	b<-rbind(d1b,d1b)
	du[[j]]<-cbind(a,b)
	colnames(du[[j]])<-c('marker1','vector1','marker2','vector2')
	}
	data_vect[[j]]<-data.frame(vect_name[[j]],vect_list[[j]],vect_name_[[j]],vect_list_[[j]])	
}
dat_vect<-do.call(rbind,data_vect)
colnames(dat_vect)<-c('marker1','vector1','marker2','vector2')
datn<-do.call(rbind,dn)
datu<-do.call(rbind,du)
datun<-do.call(rbind,dun)
if(is.matrix(datu)){for(i in 1:nrow(datu)){row.names(datu)[i]<-paste0(row.names(datu)[i],i)}
datu<-as.data.frame(datu)}
vtest<-c(as.character(datn$vector1),as.character(datn$vector2))
recom<-matrix(ncol=length(vtest),nrow=length(vtest))
for (i in 1:length(vtest)){
	for (j in 1:length(vtest)){
		vi<-as.character(vtest[i])
		vj<-as.character(vtest[j])
		r<-0
		for (k in 1:nchar(vi)){
			if(substr(vi,k,k)!=substr(vj,k,k)){r<-r+1}
			}
	if(i==j){recom[i,j]<-999}else{recom[i,j]<-r}
	}	
}
rownames(recom)<-c(as.character(datn$marker1),as.character(datn$marker2))
colnames(recom)<-c(as.character(datn$marker1),as.character(datn$marker2))
copy<-vector()
for(i in 1:ncol(recom)){
	if(length(names(subset(recom[,i],recom[,i]==0)))==0){copy<-c(copy,'NA')
	}else{copy<-c(copy,names(subset(recom[,i],recom[,i]==0)))}
	}
if(length(copy)!=length(c(as.character(datn$marker1),as.character(datn$marker2)))){print(paste0('problem colony ',c,' chromosome ',z))
	next}
copyd<-data.frame(m1=copy,m2=c(as.character(datn$marker1),as.character(datn$marker2)))

dat<-list()
i=1
while (i<nrow(datn)){
	dat<-rbind(dat,c(as.character(datn$marker1[i]),as.character(datn$marker1[i+1])),c(as.character(datn$marker2[i]),as.character(datn$marker2[i+1])),c(as.character(datn$marker1[i+1]),as.character(datn$marker1[i])),c(as.character(datn$marker2[i+1]),as.character(datn$marker2[i])))
	i=i+2
	}
colnames(dat)<-c('m1','m2')

ordera<-'NA'
x<-copyd[which(copyd$m1==tail(ordera,1))[1],]
ordera<-c(ordera,as.character(x$m2))
nb_na<-as.vector(table(copyd[copyd=='NA']))
while(length(ordera)<(nrow(datn)+(nb_na/2))){
if(x$m1%in%ordera & x$m2%in%ordera) {
	x<-subset(dat,dat[,2]==tail(ordera,1))
	ordera<-c(ordera,as.character(x[1]))
}else if (is.na(x$m1) & as.vector(table(ordera[ordera=='NA']))<((nb_na/2)-1)){
	ordera<-c(ordera,'NA')
	order_comb<-c(as.character(subset(dat_vect$marker2,dat_vect$marker1%in%ordera | dat_vect$marker2%in%ordera)),as.character(subset(dat_vect$marker1,dat_vect$marker1%in%ordera | dat_vect$marker2%in%ordera)))
	copyd<-subset(copyd,!(copyd$m2%in%order_comb))
	x<-copyd[which(copyd$m1==tail(ordera,1))[1],]
	ordera<-c(ordera,'NA',as.character(x$m2))
	x<-subset(dat,dat[,2]==tail(ordera,1))
	ordera<-c(ordera,as.character(x[1]))
}else {
	ordera<-c(ordera,as.character(x$m2))
}	
x<-copyd[which(copyd$m1==tail(ordera,1))[1],]
}

order_vect<-vector()
for (i in 1:length(ordera)){
	if(ordera[i]=='NA'|is.na(ordera[i])){o<-'NA'}else{
		if(ordera[i]%in%datn$marker1){
			o<-as.character(subset(datn$vector1,datn$marker1==ordera[i]))
		}else if(ordera[i]%in%datn$marker2){
			o<-as.character(subset(datn$vector2,datn$marker2==ordera[i]))}
		}
		order_vect<-c(order_vect,o)
	}
dat_order<-data.frame(d=rep('da',length(ordera)),ordera,order_vect)
dat_order[is.na(dat_order)]<-'NA'
h<-list()
i=2
while(i<nrow(dat_order)-1){
if(dat_order[i,2]=='NA'& dat_order[i,3]=='NA'){row1<-rep(NA,5)
row2<-rep(NA,5)
h<-rbind(h,row1)}else{
row1<-which(dat_vect$marker1==as.character(dat_order[i,2]) | dat_vect$marker2==as.character(dat_order[i,2]))
row2<-which(dat_vect$marker1==as.character(dat_order[i+1,2]) | dat_vect$marker2==as.character(dat_order[i+1,2]))
hsub<-dat_vect[row1:row2,]
if(as.character(dat_order[i,2])%in%hsub$marker1 & as.character(dat_order[i+1,2])%in%hsub$marker1){
	h<-rbind(h,hsub)
}else if (as.character(dat_order[i,2])%in%hsub$marker2 & as.character(dat_order[i+1,2])%in%hsub$marker2){
	hsubn<-data.frame(marker1=hsub$marker2,vector1=hsub$vector2,marker2=hsub$marker1,vector2=hsub$vector1)
	h<-rbind(h,hsubn)}}
i=i+2
}
h[h=='NA']<-NA
nan<-which(is.na(h[,1]))
if(length(datu)>0){
i=1
while (i<nrow(datu)){
	dir1<-datu[c(i,i+1),]
	dir2<-data.frame(marker1=datu[c(i,i+1),]$marker2,vector1=datu[c(i,i+1),]$vector2,marker2=datu[c(i,i+1),]$marker1,vector2=datu[c(i,i+1),]$vector1)
	s1<-subset(h,h$vector1==as.character(dir1$vector1[1]))
	s2<-subset(h,h$vector1==as.character(dir2$vector1[1]))
	if(nrow(s1)>0){
	d1<-data.frame(marker1=dir1[1],vector1=dir1[2],marker2=dir1[3],vector2=dir1[4])
	if(nrow(s1)>=2){h<-rbind(h[1:which(h$marker1==s1$marker1[1]),],d1[1,],h[which(h$marker1==s1$marker1[2]):nrow(h),])
		}else if (which(h$vector1==as.character(dir1$vector1[1]))==1){h<-rbind(d1[1,],h)
		}else if (which(h$vector1==as.character(dir1$vector1[1]))==nrow(h)){h<-rbind(h,d1[1,])
		}else if (is.na(h[which(h$vector1==as.character(dir1$vector1[1]))+1,'marker1'])) {h1<-h[1:which(h$vector1==as.character(dir1$vector1[1])),]
			hna<-d1[1,]
			h2<-h[(which(h$vector1==as.character(dir1$vector1[1]))+1):nrow(h),]
			h<-rbind(h1,hna,h2)}
	}else if(nrow(s2)>0){
	d2<-data.frame(marker1=dir2[1],vector1=dir2[2],marker2=dir2[3],vector2=dir2[4])
	if(nrow(s2)>=2){h<-rbind(h[1:which(h$marker1==s2$marker1[1]),],d2[1,],h[which(h$marker1==s2$marker1[2]):nrow(h),])
	}else if (which(h$vector1==as.character(dir2$vector1[1]))==1){h<-rbind(d2[1,],h)
	}else if (which(h$vector1==as.character(dir2$vector1[1]))==nrow(h)){h<-rbind(h,d2[1,])
	}else if (is.na(h[which(h$vector1==as.character(dir2$vector1[1]))+1,'marker1'])) {h1<-h[1:which(h$vector1==as.character(dir2$vector1[1])),]
			hna<-d2[1,]
			h2<-h[(which(h$vector1==as.character(dir2$vector1[1]))+1):nrow(h),]
			h<-rbind(h1,hna,h2)}
	}else if (nrow(s1)==0 & nrow(s2)==0){
		vtest<-c(as.character(datn$vector1),as.character(datn$vector2))
		vname<-c(as.character(datn$marker1),as.character(datn$marker2))
		vi<-as.character(dir1$vector1[1])
		re<-vector()
		for (j in 1:length(vtest)){
			vj<-vtest[j]
			r<-0
			for (k in 1:nchar(vi)){if(substr(vi,k,k)!=substr(vj,k,k)){r<-r+1}}
		re<-c(re,r)}
		dat_re<-data.frame(vname,re)
		k1<-subset(dat_re,re==1)
		if(nrow(k1)==0){h<-h
		}else if(k1[2,1]%in%h$marker1){h1<-h[1:which(h$marker1==as.character(k1[2,1]))-1,]
			hna<-dir1[1,]
			h2<-h[which(h$marker1==as.character(k1[2,1])):nrow(h),]
			h<-rbind(h1,hna,h2)
		}else if(k1[2,1]%in%h$marker2){h1<-h[1:which(h$marker2==as.character(k1[2,1]))-1,]
			hna<-dir2[1,]
			h2<-h[which(h$marker2==as.character(k1[2,1])):nrow(h),]
			h<-rbind(h1,hna,h2)}	
		}
	i=i+2}
}
h[h=='NA']<-NA
if(length(datun)>0){
i=1
while (i<nrow(datun)){
	dir<-dat_vect[which(dat_vect$marker2%in%c(as.character(datun$marker1[i]),as.character(datun$marker2[i]))):which(dat_vect$marker1%in%c(as.character(datun$marker1[i+1]),as.character(datun$marker2[i+1]))),]
	s1<-subset(h,h$vector1==as.character(dir$vector1[1]))
	s2<-subset(h,h$vector1==as.character(dir$vector2[1]))
	if(nrow(s1)>=2){
	h<-rbind(h[1:which(h$marker1==s1$marker1[1]),],dir,h[which(h$marker1==s1$marker1[2]):nrow(h),])
	}else if(nrow(s2)>0){
	diro<-data.frame(marker1=dir$marker2,vector1=dir$vector2,marker2=dir$marker1,vector2=dir$vector1)	
	h<-rbind(h[1:which(h$marker1==s2$marker1[1]),],diro,h[which(h$marker1==s2$marker1[2]):nrow(h),])}
	i=i+2}}	
h[h=='NA']<-NA
nan<-which(is.na(h[,1]))
if(length(nan)>0){
row_nb<-vector()
vtest<-vector()
vname<-vector()
for(i in 1:length(nan)){
	row_nb<-c(row_nb,nan[i]-1,nan[i]+1)
	vtest<-c(vtest,as.character(h[row_nb,]$vector1),as.character(h[row_nb,]$vector2))
	vname<-c(vname,as.character(h[row_nb,]$marker1),as.character(h[row_nb,]$marker2))
	}
row_nb<-c(1,row_nb,nrow(h))	
info_row<-data.frame(r=row_nb,m1=as.vector(h[row_nb,'marker1']),m2=as.vector(h[row_nb,'marker2']),end=rep(c('row_start','row_end'),length(row_nb)/2))
vtest<-c(as.character(h[1,]$vector1),as.character(h[1,]$vector2),vtest,as.character(h[nrow(h),]$vector1),as.character(h[nrow(h),]$vector2))
vname<-c(as.character(h[1,]$marker1),as.character(h[1,]$marker2),vname,as.character(h[nrow(h),]$marker1),as.character(h[nrow(h),]$marker2))
keep<-list()
for(i in 1:length(vname)){
vi<-vtest[i]
re<-vector()
	for (j in 1:length(vtest)){
		vj<-vtest[j]
		r<-0
		for (k in 1:nchar(vi)){if(substr(vi,k,k)!=substr(vj,k,k)){r<-r+1}}
	re<-c(re,r)	
	}
dat_id<-data.frame(vtest=rep(vname[i],length(vtest)),vname,re)
dat_id2<-dat_id[dat_id$re!=0,]
keep<-rbind(keep,dat_id2[which.min(dat_id2$re),])
}
k1<-subset(keep,keep$re==min(re))
hl<-list()
i=1
while(i<(length(nan)+2)){
if(i==1){n<-vname[!(vname%in%c(as.character(k1[,1]),as.character(k1[,2])))][1]}else{na<-tail(hl[[i-1]]$marker1,2)[1]
n<-subset(k1$vname,k1[,1]==as.character(na))[1]
}
infoi<-which(info_row$m1==as.character(n)|info_row$m2==as.character(n))
pos1<-subset(info_row$r,info_row$m1==as.character(n)|info_row$m2==as.character(n))
if(subset(info_row$end,info_row$m1==as.character(n)|info_row$m2==as.character(n))=='row_start'){
	pos2<-info_row[infoi+1,'r']}else{pos2<-info_row[infoi-1,'r']}
if(n%in%info_row$m1){hl[[i]]<-rbind(h[pos1:pos2,],rep(NA,5))
}else if(n%in%info_row$m2){ho<-data.frame(marker1=h[pos1:pos2,'marker2'],vector1=h[pos1:pos2,'vector2'],marker2=h[pos1:pos2,'marker1'],vector2=h[pos1:pos2,'vector1'])
	hl[[i]]<-rbind(ho,rep(NA,5))}
i=i+1}
h<-do.call(rbind,hl)}
l<-unlist(strsplit(unlist(strsplit(as.character(h$marker1),'/')),'_'))
ln<-vector()
i=0
if(TRUE%in%is.na(l)){while(i<length(l)){
	i=i+1
	if(is.na(l[i])){ln<-c(ln,NA,NA)}else{ln<-c(ln,l[i])}
	}
	}else{ln<-l}	
n<-length(ln)/2
h$contig<-ln[2*(1:n)-1]
h$bin<-ln[2*(1:n)]
contig_col<-do.call(rbind,contig_col)
contig_col[c('vector','vector_','col')]<-list(NULL)
for(i in 1:nrow(h)){
	if(is.na(h$contig[i])){
	h$nb_markers_all[i]<-NA
	h$size_all[i]<-NA
	h$position_all[i]<-NA
	h$dist_mean[i]<-NA
	h$dist_median[i]<-NA
	h$dist_var[i]<-NA
	h$density_markers[i]<-NA}
	else{sub<-subset(contig_col,contig_col$contig==h$contig[i] & contig_col$bin==h$bin[i])
	h$nb_markers_all[i]<-sub$nb_markers_all
	h$size_all[i]<-sub$size_all
	h$position_all[i]<-as.character(sub$position_all)
	h$dist_mean[i]<-sub$dist_mean
	h$dist_median[i]<-sub$dist_median
	h$dist_var[i]<-sub$dist_var
	h$density_markers[i]<-sub$density_markers}}
h[h=='NA']<-NA
h<-subset(h,!is.na(h[,1]))
if(length(unique(h$contig[!(is.na(h$contig))]))==nb_cont){print('all sorted')}else{print(paste0('missing contig: ',setdiff(chr,unique(h$contig[!(is.na(h$contig))]))))}
nchange<-vector()
for (j in 1:nrow(h)){
	vi<-as.character(h$vector1[j])
	vj<-as.character(h$vector1[j+1])
	if(is.na(vi)|is.na(vj)){r<-0}else{r<-0
	for (k in 1:nchar(vi)){
		if(substr(vi,k,k)!=substr(vj,k,k)){r<-r+1}
		}}
		nchange<-c(nchange,r)}
		h$recombination<-nchange
		h$chromosome<-z
write.table(h,paste(pattern,c,'/chromosome/chromosome',z,'.txt',sep=''),col.names=T,row.names=F,quote=F,sep=',')
hd[[z]]<-h
}
h_all<-do.call(rbind,hd)
write.table(h_all,paste(pattern,c,'/chromosome/all_chromosome.txt',sep=''),col.names=T,row.names=F,quote=F,sep=',')
} 
}else{
datadistance<-do.call(rbind,datdist2)
datadistance$col<-c(rep('1',nrow(datdist2[[1]])),rep('2',nrow(datdist2[[2]])),rep('3',nrow(datdist2[[3]])))
col<-c(rep(1,nrow(dat[[1]])),rep(2,nrow(dat[[2]])),rep(3,nrow(dat[[3]])))

databins<-rbind(dat[[1]],dat[[2]],dat[[3]])
databins$col<-col
colnames(databins)[1]<-'contig'
databins$density_markers<-databins$nb_markers_all/databins$size_all
databins$density_markers[!is.finite(databins$density_markers)]<-0
databins$dist_mean[is.na(databins$dist_mean)]<-0
databins$dist_median[is.na(databins$dist_median)]<-0
databins$dist_var[is.na(databins$dist_var)]<-0

#datacontig<-read.table(args[2],sep=',',header=T)
idx<-datacontig[grepl('geno',colnames(datacontig))]
idx<-idx[rowSums(is.na(idx)) != ncol(idx), ]
datacontig<-subset(datacontig,rownames(datacontig)%in%rownames(idx))

nb_markers1<-vector()
nb_bin1<-vector()
size_markers1<-vector()
nb_markers2<-vector()
nb_bin2<-vector()
size_markers2<-vector()
nb_markers3<-vector()
size_markers3<-vector()
nb_bin3<-vector()
for(i in 1:nrow(datacontig)){
	s1<-subset(databins,as.numeric(databins$contig)==datacontig$cont[i] & databins$col==1)
	s2<-subset(databins,as.numeric(databins$contig)==datacontig$cont[i] & databins$col==2)
	s3<-subset(databins,as.numeric(databins$contig)==datacontig$cont[i] & databins$col==3)
	nb_markers1<-c(nb_markers1,sum(s1$nb_markers_all))
	size_markers1<-c(size_markers1,sum(s1$size_all))
	nb_bin1<-c(nb_bin1,nrow(s1))
	nb_markers2<-c(nb_markers2,sum(s2$nb_markers_all))
	size_markers2<-c(size_markers2,sum(s2$size_all))
	nb_bin2<-c(nb_bin2,nrow(s2))	
	nb_markers3<-c(nb_markers3,sum(s3$nb_markers_all))
	size_markers3<-c(size_markers3,sum(s3$size_all))
	nb_bin3<-c(nb_bin3,nrow(s3))
	}
datacontig$nb_marker1<-nb_markers1
datacontig$size_markers1<-size_markers1
datacontig$nb_bin1<-nb_bin1
datacontig$nb_marker2<-nb_markers2
datacontig$size_markers2<-size_markers2
datacontig$nb_bin2<-nb_bin2
datacontig$nb_marker3<-nb_markers3
datacontig$size_markers3<-size_markers3
datacontig$nb_bin3<-nb_bin3

l1<-lm(datacontig$nb_bin1~datacontig$size_markers1)
tr1<-(l1$coef[2]/nchar(as.character(dat[[1]]$vector[1])))*10^6*10^2
l2<-lm(datacontig$nb_bin2~datacontig$size_markers2)
tr2<-(l2$coef[2]/nchar(as.character(dat[[2]]$vector[1])))*10^6*10^2
l3<-lm(datacontig$nb_bin3~datacontig$size_markers3)
tr3<-(l3$coef[2]/nchar(as.character(dat[[3]]$vector[1])))*10^6*10^2

print('colony 1')
print(paste('number of bins:',nrow(databins[databins$col==1,])),sep=' ')
print(paste('number of contigs:',length(unique(databins[databins$col==1,]$contig))),sep=' ')
print(paste('number of markers:',sum(databins[databins$col==1,]$nb_markers_all)))
print('colony 2')
print(paste('number of bins:',nrow(databins[databins$col==2,])),sep=' ')
print(paste('number of contigs:',length(unique(databins[databins$col==2,]$contig))),sep=' ')
print(paste('number of markers:',sum(databins[databins$col==2,]$nb_markers_all)))
print('colony 3')
print(paste('number of bins:',nrow(databins[databins$col==3,])),sep=' ')
print(paste('number of contigs:',length(unique(databins[databins$col==3,]$contig))),sep=' ')
print(paste('number of markers:',sum(databins[databins$col==3,]$nb_markers_all)))

pdf('bin_description.pdf',width=10,height=10)
par(mfrow=c(1,3))
plot(datacontig$size_markers1,datacontig$nb_bin1)
abline(l1,col='orange')
legend('topleft',legend=paste(round(tr1,digit=2),'cM/Mb',sep=' '))
plot(datacontig$size_markers2,datacontig$nb_bin2)
abline(l2,col='orange')
legend('topleft',tr2,legend=paste(round(tr2,digit=2),'cM/Mb',sep=' '))
plot(datacontig$size_markers3,datacontig$nb_bin3)
abline(l3,col='orange')
legend('topleft',tr3,legend=paste(round(tr3,digit=2),'cM/Mb',sep=' '))
library(ggplot2)
datadistance[is.na(datadistance)]<-0
datadistance$contig<-as.numeric(datadistance$contig)
#qplot(contig,distance,data=datadistance[datadistance$col==1,],color=bin)
#qplot(contig,distance,data=datadistance[datadistance$col==2,],color=bin)
#qplot(contig,distance,data=datadistance[datadistance$col==3,],color=bin)
par(mfrow=c(3,3))
plot(databins[databins$col==1,]$size_all,databins[databins$col==1,]$nb_markers_all,ylim=c(0,2500))
plot(databins[databins$col==2,]$size_all,databins[databins$col==2,]$nb_markers_all,ylim=c(0,2500))
plot(databins[databins$col==3,]$size_all,databins[databins$col==3,]$nb_markers_all,ylim=c(0,2500))
plot(databins[databins$size_all<10000&databins$col==1,]$size_all,databins[databins$size_all<10000&databins$col==1,]$nb_markers_all,ylim=c(0,50))
plot(databins[databins$size_all<10000&databins$col==2,]$size_all,databins[databins$size_all<10000&databins$col==2,]$nb_markers_all,ylim=c(0,50))
plot(databins[databins$size_all<10000&databins$col==3,]$size_all,databins[databins$size_all<10000&databins$col==3,]$nb_markers_all,ylim=c(0,50))
plot(databins[databins$size_all<2000&databins$col==1,]$size_all,databins[databins$size_all<2000&databins$col==1,]$nb_markers_all,ylim=c(0,50))
plot(databins[databins$size_all<2000&databins$col==2,]$size_all,databins[databins$size_all<2000&databins$col==2,]$nb_markers_all,ylim=c(0,50))
plot(databins[databins$size_all<2000&databins$col==3,]$size_all,databins[databins$size_all<2000&databins$col==3,]$nb_markers_all,ylim=c(0,50))
dev.off()

hd<-list()
for(c in 1:nb_colony){
	print(paste0('colony ',c))
for(z in 1:16){
print(paste0('chromosome ',z))
vect_name<-list()
vect_name_<-list()
vect_list<-list()
vect_list_<-list()
data_vect<-list()
dn<-list()
du<-list()
dun<-list()
contig_col<-subset(databins,databins$contig==z&databins$col==c)
	if(nrow(contig_col)>1){
	vect_list<-vector()
	vect_list_<-vector()
	vect_name<-vector()
	vect_name_<-vector()
	v=1
	vect<-contig_col$vector[v]
	vect_<-contig_col$vector_[v]
	vect_name<-paste0(contig_col$contig[v],'/',contig_col$bin[v])
	vect_name_<-paste0(contig_col$contig[v],'/',contig_col$bin[v],'_')
	while(v<nrow(contig_col)){
		vect_list<-c(vect_list,as.character(vect))
		vect_list_<-c(vect_list_,as.character(vect_))
		v=v+1
		vect_t1<-contig_col$vector[v]
		vect_t2<-contig_col$vector_[v]
		vi<-as.character(vect)
		vj<-as.character(vect_t1)
		vk<-as.character(vect_t2)		
		r1<-0
		r2<-0
		for (k in 1:nchar(vi)){
			if(substr(vi,k,k)!=substr(vj,k,k)){r1<-r1+1}
			if(substr(vi,k,k)!=substr(vk,k,k)){r2<-r2+1}
		}
		if(r1<r2){vect_name<-c(vect_name,paste0(contig_col$contig[v],'/',contig_col$bin[v]))
			vect_name_<-c(vect_name_,paste0(contig_col$contig[v],'/',contig_col$bin[v],'_'))
			vect<-vect_t1
			vect_<-vect_t2
			}else if(r1>r2){vect_name<-c(vect_name,paste0(contig_col$contig[v],'/',contig_col$bin[v],'_'))
			vect_name_<-c(vect_name_,paste0(contig_col$contig[v],'/',contig_col$bin[v]))
			vect<-vect_t2
			vect_<-vect_t1}
		}
	if(v==nrow(contig_col)){
		if(grepl('_',tail(vect_name,1))==TRUE){vect_list<-c(vect_list,as.character(contig_col[nrow(contig_col),'vector_']))
		vect_list_<-c(vect_list_,as.character(contig_col[nrow(contig_col),'vector']))}else{vect_list_<-c(vect_list_,as.character(contig_col[nrow(contig_col),'vector_']))
		vect_list<-c(vect_list,as.character(contig_col[nrow(contig_col),'vector']))}
		}
	d1a<-c(vect_name[1],vect_list[1])
	d1b<-c(vect_name_[1],vect_list_[1])	
	dna<-c(tail(vect_name,1),tail(vect_list,1))
	dnb<-c(tail(vect_name_,1),tail(vect_list_,1))
	a<-rbind(d1a,dna)
	b<-rbind(d1b,dnb)
	datn<-data.frame(a,b)
	colnames(datn)<-c('marker1','vector1','marker2','vector2')
	if(datn$vector1[1]==datn$vector1[2]){datun<-datn
	datn<-NULL}
	}else{
	vect_list<-as.character(contig_col$vector)
	vect_list_<-as.character(contig_col$vector_)
	vect_name<-paste0(contig_col$contig,'/',contig_col$bin)
	vect_name_<-paste0(contig_col$contig,'/',contig_col$bin,'_')
	d1a<-c(vect_name[1],vect_list[1])
	d1b<-c(vect_name_[1],vect_list_[1])	
	a<-rbind(d1a,d1a)
	b<-rbind(d1b,d1b)
	datu<-cbind(a,b)
	colnames(datu)<-c('marker1','vector1','marker2','vector2')
	}
dat_vect<-data.frame(vect_name,vect_list,vect_name_,vect_list_)	
colnames(dat_vect)<-c('marker1','vector1','marker2','vector2')

if(exists('datu')==TRUE){
if(is.matrix(datu)){for(i in 1:nrow(datu)){row.names(datu)[i]<-paste0(row.names(datu)[i],i)}
datu<-as.data.frame(datu)}}
vtest<-c(as.character(datn$vector1),as.character(datn$vector2))
recom<-matrix(ncol=length(vtest),nrow=length(vtest))
for (i in 1:length(vtest)){
	for (j in 1:length(vtest)){
		vi<-as.character(vtest[i])
		vj<-as.character(vtest[j])
		r<-0
		for (k in 1:nchar(vi)){
			if(substr(vi,k,k)!=substr(vj,k,k)){r<-r+1}
			}
	if(i==j){recom[i,j]<-999}else{recom[i,j]<-r}
	}	
}
rownames(recom)<-c(as.character(datn$marker1),as.character(datn$marker2))
colnames(recom)<-c(as.character(datn$marker1),as.character(datn$marker2))
copy<-vector()
for(i in 1:ncol(recom)){
	if(length(names(subset(recom[,i],recom[,i]==0)))==0){copy<-c(copy,'NA')
	}else{copy<-c(copy,names(subset(recom[,i],recom[,i]==0)))}
	}
if(length(copy)!=length(c(as.character(datn$marker1),as.character(datn$marker2)))){print(paste0('problem colony ',c,' chromosome ',z))
	next}
copyd<-data.frame(m1=copy,m2=c(as.character(datn$marker1),as.character(datn$marker2)))

dat<-list()
i=1
while (i<nrow(datn)){
	dat<-rbind(dat,c(as.character(datn$marker1[i]),as.character(datn$marker1[i+1])),c(as.character(datn$marker2[i]),as.character(datn$marker2[i+1])),c(as.character(datn$marker1[i+1]),as.character(datn$marker1[i])),c(as.character(datn$marker2[i+1]),as.character(datn$marker2[i])))
	i=i+2
	}
colnames(dat)<-c('m1','m2')

ordera<-'NA'
x<-copyd[which(copyd$m1==tail(ordera,1))[1],]
ordera<-c(ordera,as.character(x$m2))
nb_na<-as.vector(table(copyd[copyd=='NA']))
while(length(ordera)<(nrow(datn)+(nb_na/2))){
if(x$m1%in%ordera & x$m2%in%ordera) {
	x<-subset(dat,dat[,2]==tail(ordera,1))
	ordera<-c(ordera,as.character(x[1]))
}else if (is.na(x$m1) & as.vector(table(ordera[ordera=='NA']))<((nb_na/2)-1)){
	ordera<-c(ordera,'NA')
	order_comb<-c(as.character(subset(dat_vect$marker2,dat_vect$marker1%in%ordera | dat_vect$marker2%in%ordera)),as.character(subset(dat_vect$marker1,dat_vect$marker1%in%ordera | dat_vect$marker2%in%ordera)))
	copyd<-subset(copyd,!(copyd$m2%in%order_comb))
	x<-copyd[which(copyd$m1==tail(ordera,1))[1],]
	ordera<-c(ordera,'NA',as.character(x$m2))
	x<-subset(dat,dat[,2]==tail(ordera,1))
	ordera<-c(ordera,as.character(x[1]))
}else {
	ordera<-c(ordera,as.character(x$m2))
}	
x<-copyd[which(copyd$m1==tail(ordera,1))[1],]
}

order_vect<-vector()
for (i in 1:length(ordera)){
	if(ordera[i]=='NA'|is.na(ordera[i])){o<-'NA'}else{
		if(ordera[i]%in%datn$marker1){
			o<-as.character(subset(datn$vector1,datn$marker1==ordera[i]))
		}else if(ordera[i]%in%datn$marker2){
			o<-as.character(subset(datn$vector2,datn$marker2==ordera[i]))}
		}
		order_vect<-c(order_vect,o)
	}
dat_order<-data.frame(d=rep('da',length(ordera)),ordera,order_vect)
dat_order[is.na(dat_order)]<-'NA'
h<-list()
i=2
while(i<nrow(dat_order)-1){
if(dat_order[i,2]=='NA'& dat_order[i,3]=='NA'){row1<-rep(NA,5)
row2<-rep(NA,5)
h<-rbind(h,row1)}else{
row1<-which(dat_vect$marker1==as.character(dat_order[i,2]) | dat_vect$marker2==as.character(dat_order[i,2]))
row2<-which(dat_vect$marker1==as.character(dat_order[i+1,2]) | dat_vect$marker2==as.character(dat_order[i+1,2]))
hsub<-dat_vect[row1:row2,]
if(as.character(dat_order[i,2])%in%hsub$marker1 & as.character(dat_order[i+1,2])%in%hsub$marker1){
	h<-rbind(h,hsub)
}else if (as.character(dat_order[i,2])%in%hsub$marker2 & as.character(dat_order[i+1,2])%in%hsub$marker2){
	hsubn<-data.frame(marker1=hsub$marker2,vector1=hsub$vector2,marker2=hsub$marker1,vector2=hsub$vector1)
	h<-rbind(h,hsubn)}}
i=i+2
}
h[h=='NA']<-NA
nan<-which(is.na(h[,1]))
if(exists('datu')){
if(length(datu)>0){
i=1
while (i<nrow(datu)){
	dir1<-datu[c(i,i+1),]
	dir2<-data.frame(marker1=datu[c(i,i+1),]$marker2,vector1=datu[c(i,i+1),]$vector2,marker2=datu[c(i,i+1),]$marker1,vector2=datu[c(i,i+1),]$vector1)
	s1<-subset(h,h$vector1==as.character(dir1$vector1[1]))
	s2<-subset(h,h$vector1==as.character(dir2$vector1[1]))
	if(nrow(s1)>0){
	d1<-data.frame(marker1=dir1[1],vector1=dir1[2],marker2=dir1[3],vector2=dir1[4])
	if(nrow(s1)>=2){h<-rbind(h[1:which(h$marker1==s1$marker1[1]),],d1[1,],h[which(h$marker1==s1$marker1[2]):nrow(h),])
		}else if (which(h$vector1==as.character(dir1$vector1[1]))==1){h<-rbind(d1[1,],h)
		}else if (which(h$vector1==as.character(dir1$vector1[1]))==nrow(h)){h<-rbind(h,d1[1,])
		}else if (is.na(h[which(h$vector1==as.character(dir1$vector1[1]))+1,'marker1'])) {h1<-h[1:which(h$vector1==as.character(dir1$vector1[1])),]
			hna<-d1[1,]
			h2<-h[(which(h$vector1==as.character(dir1$vector1[1]))+1):nrow(h),]
			h<-rbind(h1,hna,h2)}
	}else if(nrow(s2)>0){
	d2<-data.frame(marker1=dir2[1],vector1=dir2[2],marker2=dir2[3],vector2=dir2[4])
	if(nrow(s2)>=2){h<-rbind(h[1:which(h$marker1==s2$marker1[1]),],d2[1,],h[which(h$marker1==s2$marker1[2]):nrow(h),])
	}else if (which(h$vector1==as.character(dir2$vector1[1]))==1){h<-rbind(d2[1,],h)
	}else if (which(h$vector1==as.character(dir2$vector1[1]))==nrow(h)){h<-rbind(h,d2[1,])
	}else if (is.na(h[which(h$vector1==as.character(dir2$vector1[1]))+1,'marker1'])) {h1<-h[1:which(h$vector1==as.character(dir2$vector1[1])),]
			hna<-d2[1,]
			h2<-h[(which(h$vector1==as.character(dir2$vector1[1]))+1):nrow(h),]
			h<-rbind(h1,hna,h2)}
	}else if (nrow(s1)==0 & nrow(s2)==0){
		vtest<-c(as.character(datn$vector1),as.character(datn$vector2))
		vname<-c(as.character(datn$marker1),as.character(datn$marker2))
		vi<-as.character(dir1$vector1[1])
		re<-vector()
		for (j in 1:length(vtest)){
			vj<-vtest[j]
			r<-0
			for (k in 1:nchar(vi)){if(substr(vi,k,k)!=substr(vj,k,k)){r<-r+1}}
		re<-c(re,r)}
		dat_re<-data.frame(vname,re)
		k1<-subset(dat_re,re==1)
		if(nrow(k1)==0){h<-h
		}else if(k1[2,1]%in%h$marker1){h1<-h[1:which(h$marker1==as.character(k1[2,1]))-1,]
			hna<-dir1[1,]
			h2<-h[which(h$marker1==as.character(k1[2,1])):nrow(h),]
			h<-rbind(h1,hna,h2)
		}else if(k1[2,1]%in%h$marker2){h1<-h[1:which(h$marker2==as.character(k1[2,1]))-1,]
			hna<-dir2[1,]
			h2<-h[which(h$marker2==as.character(k1[2,1])):nrow(h),]
			h<-rbind(h1,hna,h2)}	
		}
	i=i+2}
}}
h[h=='NA']<-NA
if(exists('datun')){
if(length(datun)>0){
i=1
while (i<nrow(datun)){
	dir<-dat_vect[which(dat_vect$marker2%in%c(as.character(datun$marker1[i]),as.character(datun$marker2[i]))):which(dat_vect$marker1%in%c(as.character(datun$marker1[i+1]),as.character(datun$marker2[i+1]))),]
	s1<-subset(h,h$vector1==as.character(dir$vector1[1]))
	s2<-subset(h,h$vector1==as.character(dir$vector2[1]))
	if(nrow(s1)>=2){
	h<-rbind(h[1:which(h$marker1==s1$marker1[1]),],dir,h[which(h$marker1==s1$marker1[2]):nrow(h),])
	}else if(nrow(s2)>0){
	diro<-data.frame(marker1=dir$marker2,vector1=dir$vector2,marker2=dir$marker1,vector2=dir$vector1)	
	h<-rbind(h[1:which(h$marker1==s2$marker1[1]),],diro,h[which(h$marker1==s2$marker1[2]):nrow(h),])}
	i=i+2}}}	
h[h=='NA']<-NA
nan<-which(is.na(h[,1]))
if(length(nan)>0){
row_nb<-vector()
vtest<-vector()
vname<-vector()
for(i in 1:length(nan)){
	row_nb<-c(row_nb,nan[i]-1,nan[i]+1)
	vtest<-c(vtest,as.character(h[row_nb,]$vector1),as.character(h[row_nb,]$vector2))
	vname<-c(vname,as.character(h[row_nb,]$marker1),as.character(h[row_nb,]$marker2))
	}
row_nb<-c(1,row_nb,nrow(h))	
info_row<-data.frame(r=row_nb,m1=as.vector(h[row_nb,'marker1']),m2=as.vector(h[row_nb,'marker2']),end=rep(c('row_start','row_end'),length(row_nb)/2))
vtest<-c(as.character(h[1,]$vector1),as.character(h[1,]$vector2),vtest,as.character(h[nrow(h),]$vector1),as.character(h[nrow(h),]$vector2))
vname<-c(as.character(h[1,]$marker1),as.character(h[1,]$marker2),vname,as.character(h[nrow(h),]$marker1),as.character(h[nrow(h),]$marker2))
keep<-list()
for(i in 1:length(vname)){
vi<-vtest[i]
re<-vector()
	for (j in 1:length(vtest)){
		vj<-vtest[j]
		r<-0
		for (k in 1:nchar(vi)){if(substr(vi,k,k)!=substr(vj,k,k)){r<-r+1}}
	re<-c(re,r)	
	}
dat_id<-data.frame(vtest=rep(vname[i],length(vtest)),vname,re)
dat_id2<-dat_id[dat_id$re!=0,]
keep<-rbind(keep,dat_id2[which.min(dat_id2$re),])
}
k1<-subset(keep,keep$re==min(re))
hl<-list()
i=1
while(i<(length(nan)+2)){
if(i==1){n<-vname[!(vname%in%c(as.character(k1[,1]),as.character(k1[,2])))][1]}else{na<-tail(hl[[i-1]]$marker1,2)[1]
n<-subset(k1$vname,k1[,1]==as.character(na))[1]
}
infoi<-which(info_row$m1==as.character(n)|info_row$m2==as.character(n))
pos1<-subset(info_row$r,info_row$m1==as.character(n)|info_row$m2==as.character(n))
if(subset(info_row$end,info_row$m1==as.character(n)|info_row$m2==as.character(n))=='row_start'){
	pos2<-info_row[infoi+1,'r']}else{pos2<-info_row[infoi-1,'r']}
if(n%in%info_row$m1){hl[[i]]<-rbind(h[pos1:pos2,],rep(NA,5))
}else if(n%in%info_row$m2){ho<-data.frame(marker1=h[pos1:pos2,'marker2'],vector1=h[pos1:pos2,'vector2'],marker2=h[pos1:pos2,'marker1'],vector2=h[pos1:pos2,'vector1'])
	hl[[i]]<-rbind(ho,rep(NA,5))}
i=i+1}
h<-do.call(rbind,hl)}
l<-unlist(strsplit(unlist(strsplit(as.character(h$marker1),'/')),'_'))
ln<-vector()
i=0
if(TRUE%in%is.na(l)){while(i<length(l)){
	i=i+1
	if(is.na(l[i])){ln<-c(ln,NA,NA)}else{ln<-c(ln,l[i])}
	}
	}else{ln<-l}	
n<-length(ln)/2
h$contig<-ln[2*(1:n)-1]
h$bin<-ln[2*(1:n)]
contig_col<-do.call(rbind,contig_col)
contig_col[c('vector','vector_','col')]<-list(NULL)
for(i in 1:nrow(h)){
	if(is.na(h$contig[i])){
	h$nb_markers_all[i]<-NA
	h$size_all[i]<-NA
	h$position_all[i]<-NA
	h$dist_mean[i]<-NA
	h$dist_median[i]<-NA
	h$dist_var[i]<-NA
	h$density_markers[i]<-NA}
	else{sub<-subset(contig_col,contig_col$contig==h$contig[i] & contig_col$bin==h$bin[i])
	h$nb_markers_all[i]<-sub$nb_markers_all
	h$size_all[i]<-sub$size_all
	h$position_all[i]<-as.character(sub$position_all)
	h$dist_mean[i]<-sub$dist_mean
	h$dist_median[i]<-sub$dist_median
	h$dist_var[i]<-sub$dist_var
	h$density_markers[i]<-sub$density_markers}}
h[h=='NA']<-NA
h<-subset(h,!is.na(h[,1]))
if(length(unique(h$contig[!(is.na(h$contig))]))==1){print('all sorted')}else{print(paste0('missing contig: ',setdiff(chr,unique(h$contig[!(is.na(h$contig))]))))}
nchange<-vector()
for (j in 1:nrow(h)){
	vi<-as.character(h$vector1[j])
	vj<-as.character(h$vector1[j+1])
	if(is.na(vi)|is.na(vj)){r<-0}else{r<-0
	for (k in 1:nchar(vi)){
		if(substr(vi,k,k)!=substr(vj,k,k)){r<-r+1}
		}}
		nchange<-c(nchange,r)}
		h$recombination<-nchange
		h$chromosome<-z
write.table(h,paste(pattern,c,'/chromosome/chromosome',z,'.txt',sep=''),col.names=T,row.names=F,quote=F,sep=',')
hd[[z]]<-h
}
h_all<-do.call(rbind,hd)
write.table(h_all,paste(pattern,c,'/chromosome/all_chromosome.txt',sep=''),col.names=T,row.names=F,quote=F,sep=',')
} 
}

################################################
### script to order contigs when 3 colonies in 1
################################################
#!/usr/bin/env R
args<-commandArgs(TRUE)
contig_info<-read.table(args[1],sep=',',header=T)
cont<-subset(contig_info,contig_info$geno=='y')
cont<-cont[complete.cases(cont),]
options(warn=2)

order<-list()
lwrong<-vector()
for(z in 1:max(cont$chr_a)){
lcont<-cont[cont$chr_a==z,'cont']
bkeep<-list()
recomb<-vector()
orient<-vector()
orient2<-vector()
bins<-vector()
size<-vector()
for(x in 1:length(lcont)){
b<-read.table(paste('bins2/bins_',lcont[x],'.txt',sep=''),sep=';',header=T,colClasses=c(rep('numeric',4),'character','numeric','numeric','character'))	
if(nrow(b)>1){
rec<-vector()
for(i in 1:(nrow(b)-1)){
r<-0
vi<-b$vector[i]
vj<-b$vector[i+1]
for (k in 1:nchar(vi)){
if(substr(vi,k,k)!=substr(vj,k,k)){r<-r+1}}
rec<-c(rec,r)}
recomb[x]<-sum(rec)
bins[x]<-nrow(b)
marker<-paste0(b$list_markers,collapse=',')
marker<-gsub('\\/*','',marker)
mmin<-unlist(strsplit(marker,','))[1]
mmax<-unlist(strsplit(marker,','))[length(unlist(strsplit(marker,',')))]
size[x]<-abs(as.numeric(mmin)-as.numeric(mmax))
if(as.numeric(mmin)<as.numeric(mmax)){orient[x]<-'+'}else if (as.numeric(mmin)>as.numeric(mmax)){orient[x]<-'-'}else{orient[x]<-'0'} 
bkeep[[x]]<-data.frame('contig'=rep(lcont[x],2),'start_end'=c(1,2),'vector'=b[c(1,nrow(b)),'vector'],'orient'=orient[x])
}else{
rec<-vector()
for(i in 1:(nrow(b)-1)){
recomb[x]<-0
bins[x]<-1
orient[x]<-'0'
size[x]<-1
bkeep[[x]]<-data.frame('contig'=lcont[x],'start_end'=1,'vector'=b[,'vector'],'orient'=orient[x])}
}}
data_keep<-do.call(rbind,bkeep)
data_keep$id<-paste0(data_keep$contig,'_',data_keep$start_end)
mrec<-matrix(ncol=nrow(data_keep),nrow=nrow(data_keep))
for(i in 1:nrow(data_keep)){
	for(j in 1:nrow(data_keep)){
		r<-0
		vi<-as.character(data_keep[i,'vector'])
		vj<-as.character(data_keep[j,'vector'])
		for (k in 1:nchar(vi)){if(substr(vi,k,k)!=substr(vj,k,k)){r<-r+1}}
		mrec[i,j]<-r}}
diag(mrec)<-999
m<-vector()
for(i in 1:nrow(data_keep)){m<-c(m,min(mrec[,i]))}
data_newstart<-data_keep[which.max(m),'contig']
dk<-data_keep[data_keep$contig==data_newstart,]
if(nrow(bkeep[[which(lcont==data_newstart)]])==2){datanew<-rbind(dk[rownames(dk)==which.max(m),],dk[rownames(dk)!=which.max(m),])
dataleft<-data_keep[!data_keep$contig==data_newstart,]
i=nrow(datanew)}else{datanew<-data_keep[1,]
dataleft<-data_keep[!data_keep$contig==data_newstart,]
i=nrow(datanew)}

while(i<nrow(data_keep)){
if(datanew[i,'contig']%in%subset(dataleft$contig,!(dataleft$id%in%datanew$id))){
	datanew<-rbind(datanew,dataleft[dataleft$contig==datanew[i,'contig'],])
	m<-as.numeric(rownames(dataleft[dataleft$contig==datanew[i,'contig'],]))
	dataleft<-dataleft[rownames(dataleft)!=m,]
	i=nrow(datanew)}
else{vi<-as.character(datanew[i,'vector'])
	rec<-vector()
	for(j in 1:nrow(dataleft)){
		vj<-as.character(dataleft[j,'vector'])
		r<-0		
		for (k in 1:nchar(vi)){if(substr(vi,k,k)!=substr(vj,k,k)){r<-r+1}}
		rec<-c(rec,r)}
	n<-which.min(rec)
	datanew<-rbind(datanew,dataleft[n,])
	dataleft<-dataleft[-n,]
	i=nrow(datanew)}}
data_orient<-subset(datanew,datanew$orient!=0)
data0<-subset(datanew,datanew$orient==0)

datanew2<-data_orient
if(nrow(data0)!=0){
for(o in 1:nrow(data0)){v<-as.character(data0[o,'vector'])
	b<-which(datanew2$vector==v)
	if(length(b)==2){datanew2<-rbind(datanew2[1:b[1],],data0[o,],datanew2[b[2]:nrow(datanew2),])
	}else if (length(b)==1){
		if(datanew2[(b-1),'contig']==datanew2[b,'contig'] & b==nrow(datanew2)){datanew2<-rbind(datanew2[1:b,],data0[o,])
		}else if(datanew2[(b-1),'contig']==datanew2[b,'contig'] & b<nrow(datanew2)){datanew2<-rbind(datanew2[1:b,],data0[o,],datanew2[(b+1):nrow(datanew2),])
		}else{datanew2<-rbind(datanew2[1:(b-1),],data0[o,],datanew2[b:nrow(datanew2),])}		
	}else if (length(b)==0){
			vi<-as.character(data0[o,'vector'])
			rec<-vector()
			for(j in 1:nrow(datanew2)){
			vj<-as.character(datanew2[j,'vector'])
			r<-0		
			for (k in 1:nchar(vi)){if(substr(vi,k,k)!=substr(vj,k,k)){r<-r+1}}
			rec<-c(rec,r)}
			n<-which(rec==min(rec))
			if(length(n)==2){datanew2<-rbind(datanew2[1:n[1],],data0[o,],datanew2[n[2]:nrow(datanew2),])}else if (length(n)==1){datanew2<-rbind(datanew2[1:(n-1),],data0[o,],datanew2[n:nrow(datanew2),])}
	}else{
		if(datanew2[(min(b)-1),'contig']==datanew2[min(b),'contig'] & datanew2[(max(b)+1),'contig']==datanew2[max(b),'contig']){datanew2<-rbind(datanew2[1:b[length(b)-1],],data0[o,],datanew2[b[length(b)]:nrow(datanew2),])
		}else if (min(b)==1){datanew<-rbind(data0[o,],datanew)
		}else if (max(b)==nrow(datanew)){datanew<-rbind(datanew,data0[o,])
		}else{datanew2<-rbind(datanew2[1:(min(b)-1),],data0[o,],datanew[min(b):nrow(datanew),])}
		}
}
}

l<-cont[cont$chr_a==z,]
lcont_order<-l[order(l$pos),'cont']
if(datanew2$contig[1]==lcont_order[1]|datanew2$contig[nrow(datanew2)]==lcont_order[length(lcont_order)]){datanew2<-datanew2}else{datanew2<-datanew2[nrow(datanew2):1,]}

lcont<-unique(datanew2$contig)
for(x in 1:length(lcont)){
a<-subset(datanew2,datanew2$contig==lcont[x])
if(nrow(a)>1 & a$start_end[1]==1 & a$start_end[2]==2){orient2[x]<-'ok'} else if (nrow(a)>1 & a$start_end[1]==2 & a$start_end[2]==1){orient2[x]<-'wrong'}else if (nrow(a)==1){orient2[x]<-'0'}}
orientation<-vector()
for(x in 1:length(lcont)){
	o<-unique(datanew2$orient[datanew2$contig==lcont[x]])
	if(orient2[x]=='ok'){orientation<-c(orientation,as.character(o))
	}else if (orient2[x]=='wrong' & as.character(o)=='+'){orientation<-c(orientation,'-')
														lwrong<-c(lwrong,lcont[x])
	}else if (orient2[x]=='wrong' & as.character(o)=='-'){orientation<-c(orientation,'+')
														lwrong<-c(lwrong,lcont[x])
	}else if (as.character(o)=='0'&orient2[x]=='0'){orientation<-c(orientation,'0')}}
dat_nonorder<-data.frame('chr'=z,'contig'=lcont,'orientation'=orientation,'nb_bins'=bins,'size'=size,'nb_recom'=recomb)
dato<-dat_nonorder[order(match(dat_nonorder$contig,unique(datanew2$contig))),]
order[[z]]<-data.frame(dato,'order'=seq(1:length(unique(datanew2$contig))))
}
	
final1<-do.call(rbind,order)
colnames(final1)[2]<-'cont'	
final<-merge(cont,final1,by='cont')
write.table(final,'bins2/contigonchromosomes.txt',sep=',',col.names=T,row.names=F,quote=F)	

order_all<-do.call(rbind,order)
dat<-data.frame()
for(i in 1:nrow(order_all)){
	n<-order_all[i,'contig']
	bin<-read.table(paste('bins2/bins_',n,'.txt',sep=''),sep=';',header=T,colClasses=c(rep('numeric',4),'character','numeric','numeric','character'))	
	if(n %in% lwrong){bin<-bin[rev(rownames(bin)),]}
	bin$marker1<-paste(bin$contig,bin$bin,sep='/')
	bin$vector1<-bin$vector
	bin$marker2<-paste(bin$contig,bin$bin,sep='/')
	bin$vector2<-bin$vector
	colnames(bin)[6]<-'nb_markers_all'
	colnames(bin)[7]<-'size_all'	
	colnames(bin)[3]<-'chromosome'	
	bin$density_markers<-bin$nb_markers_all/bin$size_all
	for(j in 1:nrow(bin)){
		bin$list_markers[j]<-gsub('\\/*','',bin$list_markers[j])
		l<-as.numeric(unlist(strsplit(bin$list_markers[j],',')))
		dist<-vector()
		for(k in 1:(length(l)-1)){dist<-c(dist,abs(l[k]-l[k+1]))}
		bin$dist_mean[j]<-mean(dist)
		bin$dist_median[j]<-median(dist)
		bin$dist_var[j]<-var(dist)
		vi<-as.character(bin$vector[j-1])
		vj<-as.character(bin$vector[j])
		if(j==1){bin$recombination[j]<-0}else{
		r<-0		
		for (k in 1:nchar(vi)){if(substr(vi,k,k)!=substr(vj,k,k)){r<-r+1}}
		bin$recombination[j]<-r}
		bin$pos_m[j]<-unlist(strsplit(bin$list_markers[j],','))[1]
		bin$pos_M[j]<-unlist(strsplit(bin$list_markers[j],','))[length(unlist(strsplit(bin$list_markers[j],',')))]}
dat<-rbind(dat,bin)}
dat<-dat[,c('marker1','vector1','marker2','vector2','contig','bin','nb_markers_all','size_all','dist_mean','dist_median','dist_var','density_markers','recombination','chromosome','pos_m','pos_M')]

rec<-list()
for(z in 1:max(cont$chr_a)){
datz<-dat[dat$chromosome==z,]
reco<-vector()
for(j in 1:(nrow(datz)-1)){	
vj<-datz$vector1[j]
vk<-datz$vector1[j+1]
r<-0
for (k in 1:nchar(vj)){if(substr(vj,k,k)!=substr(vk,k,k)){r<-r+1}}
reco<-c(reco,r)}
rec[[z]]<-c(0,reco)
}
dat$recombination<-unlist(rec)
write.table(dat,'bins2/chromosome/all_chr_new.txt',sep=',',col.names=T,row.names=F,quote=F)

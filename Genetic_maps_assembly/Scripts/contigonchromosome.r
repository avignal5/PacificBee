################################################
### script to order contigs on chromosome according to prior knowledge
################################################
#!/usr/bin/env R
library('data.table')
library('gplots')
library('ggplot2')
library('colorout')
args<-commandArgs(TRUE)
nb_colony<-as.numeric(args[1])
pattern<-args[3]

#contonchr<-read.table(args[2],header=T,sep=',',stringsAsFactor=F)
#contsize<-read.table(args[3])
#colnames(contsize)<-c('contig_id','size')
#contorient<-read.table(args[5],header=T,sep='\t')
#colnames(contorient)<-c('contig_id','chr_lg','orientation_a','cont_start','cont_end','cont_length','chr_start','chr_end','chr_length')
#conti<-merge(contonchr,contorient,by='contig_id',all=TRUE)
#conti<-merge(conti,contorient,by='contig_id',all=TRUE)
#conti$chr_a[is.na(conti$chr_a)]<-'Un'
conti<-read.table(args[2],header=T,sep=',',stringsAsFactors=F)
nbchr<-max(as.numeric(conti$chr_a[!is.na(conti$chr_a)&grepl('Un',conti$chr_a)==FALSE]))

if(nb_colony>1){
data<-list()
for(c in 1:nb_colony){
chr<-read.table(paste(pattern,c,'/chromosome/all_chr_new.txt',sep=''),header=T,sep=',',colClasses=c('character','character','character','character','character','character','character','character','character','character','character','character','character','character','character','character'))	
chr<-chr[!rowSums((is.na(chr))),]
chr<-subset(chr,chr$chromosome%in%seq(1:nbchr))
print(paste0('colony ',c))
chrn<-list()
for (i in 1:length(unique(chr$chromosome))){
datac<-subset(chr,chr$chromosome==i)
datac1<-subset(datac,datac$contig%in%conti$cont[conti$chr==i])
datac2<-subset(datac,!(datac$contig%in%conti$cont[conti$chr==i]))
ord_cont1<-as.numeric(unique(datac1$contig))
contic<-subset(conti,conti$cont%in%ord_cont1 & conti$chr==i)
contico<-contic[order(contic$pos),'cont']
if(contico[1]==tail(ord_cont1,1)|ord_cont1[1]==tail(contico,1)){datac1<-datac1[rev(rownames(datac1)),]
	datac2<-datac2[rev(rownames(datac2)),]}
ord_cont2<-as.numeric(unique(datac2$contig))
contic<-subset(conti,conti$cont%in%ord_cont2)
datacm<-datac1
if(nrow(datac2)>0){
for(k in 1:length(ord_cont2)){
		kd<-subset(datac2,datac2$contig==ord_cont2[k])
		if(nrow(kd)==1){n<-which(datacm$vector1==kd$vector1)
			if(length(n)==2){datacm<-rbind(datacm[1:n[1],],kd,datacm[n[2]:nrow(datacm),])
			}else if (length(n)==1){if(n==nrow(datacm)){datacm<-rbind(datacm,kd)}else if (n==1){datacm<-rbind(kd,datacm)}
			}else if (length(n)>2){if(nrow(datacm)%in%n){datacm<-rbind(datacm,kd)}else{datacm<-rbind(datacm[1:n[2],],kd,datacm[n[3]:nrow(datacm),])}}
		}else{
		na<-which(datacm$vector1==kd$vector1[1])
		nb<-which(datacm$vector1==kd$vector1[nrow(kd)])
		n<-min(na,nb)
		N<-max(na,nb)
		if(n%in%nb){n1<-nb;n2<-na}else{n1<-na;n2<-nb}
		if(length(n1)==0){if(max(n2)==nrow(datacm)){datacm<-rbind(datacm,kd[rev(rownames(kd)),])}else if (min(n2)==1){datacm<-rbind(kd[rev(rownames(kd)),],datacm)}
		}else if (length(n2)==0){if(max(n1)==nrow(datacm)){datacm<-rbind(datacm,kd)}else if (min(n1)==1){datacm<-rbind(kd,datacm)}	
		}else if (length(n1)==1&length(n2)==1){if(n2==(n1+1)){datacm<-rbind(datacm[1:n1,],kd,datacm[n2:nrow(datacm),])}else if (n1==(n2+1)){datacm<-rbind(datacm[1:n2,],kd[rev(rownames(kd)),],datacm[n1:nrow(datacm),])}
		}else if (max(n1)<min(n2)){datacm<-rbind(datacm[1:max(n1),],kd,datacm[min(n2):nrow(datacm),])
		}else if (max(n1)>min(n2)){datacm<-rbind(datacm[1:min(n1),],kd[rev(rownames(kd)),],datacm[(min(n1)+1):nrow(datacm),])}
		}
}}
rec<-vector()
j<-1
while(j<nrow(datacm)){
vj<-datacm$vector1[j]
vk<-datacm$vector1[j+1]
r<-0
for (k in 1:nchar(vj)){if(substr(vj,k,k)!=substr(vk,k,k)){r<-r+1}}
rec<-c(rec,r)
j=j+1}
rec<-c(0,rec)
datacm$recombination<-rec
chrn[[i]]<-datacm}
chrn<-do.call(rbind,chrn)
fwrite(chrn,paste(pattern,c,'/chromosome/all_chr_new2.txt',sep=''),na=NA)

pdf(paste0(pattern,c,'/chromosome/heatmap_colony',c,'.pdf'),width=10,height=10)
for (i in 1:length(unique(chrn$chromosome))){
chri<-subset(chrn,chrn$chromosome==i)
testchr<-do.call(rbind,strsplit(chri[,'vector1'],''))
space<-vector()
space[1]<-as.numeric(chri$size_all[1])
for(i in 1:(nrow(chri)-1)){
	space<-c(space,space[i]+as.numeric(chri$size_all[i+1]))
}
nind<-nchar(chri$vector1[1])
nmark<-nrow(testchr)
indc<-vector()
for(i in 1:nind){
	indc<-c(indc,rep(paste0('ind',i),nmark))
	}
markc<-rep(seq(1:nrow(testchr)),nind)
dat<-data.frame(indc,markc,c(testchr))
print(ggplot(dat,aes(x=dat[,1],y=dat[,2],fill=dat[,3]))+geom_tile()+scale_y_reverse())
}
dev.off()

d<-list()
for(i in 1:length(unique(chrn$chromosome))){
	sub<-subset(chrn,chrn$chromosome==i)
	di<-list()
	for(j in 1:length(unique(sub$contig))){
		sub2<-subset(sub,sub$contig==unique(sub$contig)[j])
		chromosome<-i
		contig<-unique(sub$contig)[j]
		if (sub2[1,'vector1']==sub2[nrow(sub2),'vector1']){orientation<-'0'}else if(as.numeric(sub2[1,'pos_m'])<as.numeric(sub2[nrow(sub2),'pos_m'])){orientation<-'+'}else if(as.numeric(sub2[1,'pos_m'])>as.numeric(sub2[nrow(sub2),'pos_m'])){orientation<-'-'}else if(nrow(sub2)==1){orientation<-'0'}
		nb_bins<-nrow(sub2)
		nb_markers<-sum(as.numeric(sub2$nb_markers_all))
		size<-sum(as.numeric(sub2$size_all))
		recombination<-sum(as.numeric(sub2$recombination))
		chr_cont<-paste0(chromosome,'/',contig)
		di[[j]]<-c(chromosome,contig,orientation,nb_bins,nb_markers,size,recombination,chr_cont)
		}
	d[[i]]<-do.call(rbind,di)
	}
data[[c]]<-do.call(rbind,d)
colnames(data[[c]])<-c('chromosome','contig','orientation','nb_bins','nb_markers','size','recombination','chr_cont')
data[[c]]<-as.data.frame(data[[c]])

for(i in 1:nrow(conti)){
	contig_n<-conti$cont[i]
	if(contig_n%in%data[[c]]$contig){
	d<-subset(data[[c]],data[[c]]$contig==contig_n)	
	conti[i,paste0('chr',c)]<-as.numeric(as.character(d$chromosome))
	conti[i,paste0('size',c)]<-as.numeric(as.character(d$size))
	conti[i,paste0('orientation',c)]<-as.character(d$orientation)
	}else {
	conti[i,paste0('chr',c)]<-'Un'
	conti[i,paste0('size',c)]<-0
	conti[i,paste0('orientation',c)]<-0}}
}
datafinal<-Reduce(function(x, y) merge(x, y,all=TRUE,by='chr_cont'), data)
datafinal[datafinal=='NULL']<-NA
#fwrite(datafinal,'contigonchromosomes.txt',na=NA)

conti$chr_lg[is.na(conti$chr_lg)]<-'Un'
conti$orientation_a[is.na(conti$orientation_a)]<-0
conti$cont_length[is.na(conti$cont_length)]<-0
conti$chr_length[is.na(conti$chr_length)]<-0
conti$chr_a[is.na(conti$chr_a)]<-'Un'
csub<-conti[grepl('orientation',colnames(conti))]
csub<-csub[,3:ncol(csub)]
for(i in 1:nrow(csub)){
	u<-unique(unname(unlist(csub[i,])))
	if(length(u)==1){conti[i,'orientation_']<-u[1]}else if(length(u)>1 & '0'%in%u){conti[i,'orientation_']<-u[u!='0']}else if (length(u)>1 & '+'%in%u & '-'%in%u){conti[i,'orientation_']<-'Un'}
	}
csub<-conti[grepl('chr',colnames(conti))]
csub$chr_lg<-NULL
csub$chr_start<-NULL
csub$chr_end<-NULL
csub$chr_length<-NULL
csub$chr_c<-NULL
csub$chr_a<-NULL
for(i in 1:nrow(csub)){
	u<-unique(unname(unlist(csub[i,])))
	if(length(u)==1 & conti$orientation_[i]!=0){conti[i,'chr']<-u[1]
	}else if (length(u)==1 & conti$orientation_[i]==0 & TRUE%in%grepl('Un',u)){conti[i,'chr']<-'Un'
	}else if (length(u)==1 & conti$orientation_[i]==0 & FALSE%in%grepl('Un',u)){conti[i,'chr']<-paste0('Un',u[1])
	}else if(length(u)>1 & TRUE%in%grepl('Un',u) & conti$orientation_[i]!=0){conti[i,'chr']<-u[grepl('Un',u)==F]
	}else if(length(u)>1 & FALSE%in%grepl('Un',u) & conti$orientation_[i]==0){conti[i,'chr']<-paste0('Un',u[grepl('Un',u)==F])
	}else{conti[i,'chr']<-'Un'}}
	
ord_c<-list()
for(c in 1:nb_colony){
ord<-list()
for(z in 1:nbchr){
sub_ordre<-subset(data[[c]],data[[c]]$contig%in%conti$cont[conti$chr==z]|data[[c]]$contig%in%conti$cont[conti$chr==paste0('Un',z)])
sub_conti<-subset(conti,conti$chr==z|conti$chr==paste0('Un',z))
scontio<-sub_conti[match(unlist(as.numeric(as.character(sub_ordre$contig))),sub_conti$cont),]
scontio_un<-which(scontio$chr==paste0('Un',z))
a<-split(scontio_un, cumsum(c(1, diff(scontio_un) != 1)))
n<-vector()
for(l in 1:length(a)){if(length(a[[l]])>1){n<-c(n,l)
m<-min(a[[n]])
M<-max(a[[n]])
s_un<-scontio[m:M,]
s_un<-s_un[order(s_un$pos),]
if(m==1){scontion<-rbind(s_un,scontio[(M+1):nrow(scontio),])}else if(M==nrow(scontio)){scontion<-rbind(scontio[1:(m-1),],s_un)}else{scontion<-rbind(scontio[1:(m-1),],s_un,scontio[(M+1):nrow(scontio),])}
scontio<-scontion}}
scontio_nun<-scontio[scontio$chr==z,]
if(scontio_nun$pos[1]<scontio_nun$pos[nrow(scontio_nun)]){ord[[z]]<-data.frame(chr_nb=z,chr=scontio$chr,cont=scontio$cont,ordre=seq(1:nrow(scontio)),orientation_correct=unlist(scontio$orientation_))
}else{orientation_correcte<-vector()
	for(k in 1:nrow(scontio)){
		if(unlist(scontio$orientation_)[k]=='+'){orientation_correcte[k]<-'-'}else if (unlist(scontio$orientation_)[k]=='-'){orientation_correcte[k]<-'+'} else if (unlist(scontio$orientation_)[k]=='0'){orientation_correcte[k]<-'0'}
		}
	ord[[z]]<-data.frame(chr_nb=z,chr=scontio$chr,cont=scontio$cont,ordre=rev(seq(1:nrow(scontio))),orientation_correct=orientation_correcte)}
}
ord_c[[c]]<-do.call(rbind,ord)}
ord_c<-Reduce(function(x, y) merge(x, y,all=TRUE,by='cont'), ord_c)
names(ord_c)[names(ord_c) == "chr_nb"] <-"chr_nb.z"
names(ord_c)[names(ord_c) == "chr"] <-"chr.z"
names(ord_c)[names(ord_c) == "orientation_correct"] <-"orientation_correct.z"
names(ord_c)[names(ord_c) == "ordre"] <-"ordre.z"
contif<-merge(conti,ord_c,by='cont',all=T)

options(warn=2)
ordre_keep<-list()
for(z in 1:nbchr){
chrnb<-colnames(contif)[grepl('chr_nb',colnames(contif))]
s<-data.frame()
for(y in 1:nrow(contif)){if(z%in%contif[y,chrnb]){s<-rbind(s,contif[y,])}}
col<-c(which(grepl('orientation',colnames(s))),which(grepl('ordre',colnames(s))))
s1<-s[,c(1,col)]
x<-which(grepl('ordre',colnames(s1)))[1]
s1<-s1[order(s1[,x]),]
s1$orientation<-NULL
s1$orientation_<-NULL
coly<-which(grepl('orientation_correct',colnames(s1)))
s1[,coly]<-NULL
ordre_mean<-vector()
for(y in 1:nrow(s1)){
v<-as.numeric(s1[y,grepl('ordre',colnames(s1))])
if(NA%in%v){om<-NA}else{om<-as.numeric(median(v[!is.na(v)]))}
ordre_mean<-c(ordre_mean,om)}
s1$ordre_mean<-ordre_mean
dup<-as.numeric(names(table(s1$ordre_mean)[table(s1$ordre_mean)>1]))
if(length(dup)>0){
s_dup<-subset(s1,round(s1$ordre_mean,digit=2)==round(dup,digit=2))
order_col<-vector()
l<-c('x','y','z')
for(col in 1:nb_colony){if(FALSE%in%(unique(s_dup[,paste0('orientation',col)])==0)){order_col<-s_dup[,c('cont',paste0('ordre.',l[col]))]}}
ordcol<-order_col[order(order_col[,2]),'cont']
ordcol_ord<-seq(from=dup,to=dup+1,by=1/(length(ordcol)-1))
ord_c<-data.frame(ordcol,ordcol_ord)
for(p in 1:nrow(ord_c)){s1[s1$cont==ord_c[p,1],'ordre_mean']<-ord_c[p,2]}}
if(NA%in%s1$ordre_mean){
s1nna<-subset(s1,!is.na(s1$ordre_mean))
s1nna<-s1nna[order(s1nna$ordre_mean),]
s1n<-s1nna
s1na<-subset(s1,is.na(s1$ordre_mean))
for(y in 1:nrow(s1na)){
	s_cont<-subset(contif,contif$cont==s1na$cont[y])
	v<-as.numeric(s1na[y,grepl('ordre',colnames(s1na))])
	n<-which(!is.na(v))
	v_x<-v[!is.na(v)]	
	m<-vector()
	M<-vector()
	if(1%in%v_x){s1n<-rbind(s1na[y,],s1n)
	}else if (max(s1n[,grepl('ordre',colnames(s1n))])%in%v_x){s1n<-rbind(s1n,s1na[y,])
	}else{
		for(u in 1:length(n)){
		m<-c(m,max(which(s1n[,grepl('ordre',colnames(s1n))][,n[u]]<v[n[u]])))	
		M<-c(M,min(which(s1n[,grepl('ordre',colnames(s1n))][,n[u]]>v[n[u]])))
		}
		m<-sort(m)
		M<-sort(M)	
	if(length(unique(m))==1&length(unique(M))==1){m<-m
		M<-M
	}else{if(length(unique(m))>1 & length(unique(M))==1){
		m<-M-1
		#a<-s1n[m,grepl('orientation',colnames(s1n))]
		#a0<-vector()
		#for(b in 1:nrow(a)){a0<-c(a0,length(grep('0',a[b,])))}
		#if(length(unique(a0))==1){m<-M-1}else{m<-m[which.min(a0)]}
		}else if (length(unique(m))==1 & length(unique(M))>1){
		M<-m+1
		#a<-s1n[M,grepl('orientation',colnames(s1n))]
		#a0<-vector()
		#for(b in 1:nrow(a)){a0<-c(a0,length(grep('0',a[b,])))}
		#if(length(unique(a0))==1){M<-m+1}else{M<-M[which.min(a0)]}
		}else if (length(unique(m))>1 & length(unique(M))>1){}
		a<-s1n[m,grepl('orientation',colnames(s1n))]
		a0<-vector()
		for(b in 1:nrow(a)){a0<-c(a0,length(grep('0',a[b,])))}
		m<-m[which.min(a0)]
		a<-s1n[M,grepl('orientation',colnames(s1n))]
		a0<-vector()
		for(b in 1:nrow(a)){a0<-c(a0,length(grep('0',a[b,])))}
		M<-M[which.min(a0)]
		}
s1n<-rbind(s1n[1:unique(m),],s1na[y,],s1n[unique(M):nrow(s1n),])}}		
ordre_keep[[z]]<-data.frame(cont=s1n$cont,ordre_consensus=seq(from=1,to=nrow(s1n),by=1))
}else{ordre_keep[[z]]<-data.frame(cont=s1$cont,ordre_consensus=seq(from=1,to=nrow(s1),by=1))}
}
ordre<-do.call(rbind,ordre_keep)
contif_2<-merge(contif,ordre,by='cont',all=T)
chr_nb_consensus<-vector()
for(l in 1:nrow(contif_2)){
	c<-unlist(contif_2[l,grepl('chr_nb',colnames(contif_2))])
	if(TRUE%in%!is.na(c)){chr_nb_consensus<-c(chr_nb_consensus,unique(c)[!is.na(unique(c))])}else{chr_nb_consensus<-c(chr_nb_consensus,NA)}}	
contif_2$chr_nb_consensus<-chr_nb_consensus	

cna<-subset(contif_2,is.na(contif_2$chr_nb_consensus))
cz<-list()
options(warn=2)
for(z in 1:nbchr){
	s<-subset(contif_2,contif_2$chr_nb_consensus==z)
	s<-s[order(s$ordre_consensus),]
	for(i in 1:nrow(s)){if(grepl('Un',s[i,'chr_c'])==TRUE & s[i,'chr_a']=='Un'){s[i,'chr_nb_consensus']<-paste0('Un',z)}}
	snun<-subset(s,s$chr_nb_consensus!=paste0('Un',z))
	snun<-snun[order(as.numeric(snun$ordre_consensus)),]
	sun<-subset(s,s$chr_nb_consensus==paste0('Un',z))
	sun<-sun[order(as.numeric(sun$ordre_consensus)),]
	if(nrow(sun)>0){
		if(max(sun$ordre_consensus)<min(snun$ordre_consensus)){
			sun$ordre_consensus<-paste0('0_1_',sun$ordre_consensus)
			snun$ordre_consensus<-seq(1:nrow(snun))
		}else if (min(sun$ordre_consensus)>max(snun$ordre_consensus)){
			sun$ordre_consensus<-paste0(nrow(snun),'_',nrow(snun)+1,'_',sun$ordre_consensus)
			snun$ordre_consensus<-seq(1:nrow(snun))
		}else{		
		for(i in 1:nrow(sun)){if(sun$ordre_consensus[i]==1){m<-'x'
		M<-snun[snun$ordre_consensus==as.numeric(sun$ordre_consensus[i])+1,'cont']
		}else if (sun$ordre_consensus[i]>=max(snun$ordre_consensus)){
		m<-snun[snun$ordre_consensus==as.numeric(sun$ordre_consensus[i])-1,'cont']
		M<-'x'
		}else{m<-snun[snun$ordre_consensus==as.numeric(sun$ordre_consensus[i])-1,'cont']
		M<-snun[snun$ordre_consensus==as.numeric(sun$ordre_consensus[i])+1,'cont']}
		sun$ordre_consensus[i]<-paste0(m,'_',M)
		}
		snun$ordre_consensus<-seq(1:nrow(snun))
		for(i in 1:nrow(sun)){
			c1<-unlist(strsplit(sun$ordre_consensus[i],'_'))[1]
			c2<-unlist(strsplit(sun$ordre_consensus[i],'_'))[2]
			if(c1=='x'){sun$ordre_consensus[i]<-paste0('0_',snun$ordre_consensus[snun$cont==c2])
			}else if (c2=='x'){sun$ordre_consensus[i]<-paste0(snun$ordre_consensus[snun$cont==c1],'_',(nrow(snun)+1))
			}else{sun$ordre_consensus[i]<-paste0(snun$ordre_consensus[snun$cont==c1],'_',snun$ordre_consensus[snun$cont==c2])}}}	
	if(nrow(snun)==1){snun$orientation_<-'+'}
	cz[[z]]<-rbind(snun,sun)
	}else{cz[[z]]<-snun}
	}
contif_2n<-do.call(rbind,cz)
contif_2n<-rbind(cna,contif_2n)
contif_2n$chr<-NULL
contif_2n$chr1<-NULL
contif_2n$chr2<-NULL
contif_2n$chr3<-NULL
contif_2n$orientation1<-NULL
contif_2n$orientation2<-NULL
contif_2n$orientation3<-NULL
names(contif_2n)[names(contif_2n) == "pos"] <-"pos_c"
names(contif_2n)[names(contif_2n) == "orientation"] <-"orientation_c"
names(contif_2n)[names(contif_2n) == "orientation_correct.x"] <-"orientation1"
names(contif_2n)[names(contif_2n) == "orientation_correct.y"] <-"orientation2"
names(contif_2n)[names(contif_2n) == "orientation_correct.z"] <-"orientation3"
names(contif_2n)[names(contif_2n) == "orientation_"] <-"orientation_consensus"
names(contif_2n)[names(contif_2n) == "chr_nb.x"] <-"chr_nb1"
names(contif_2n)[names(contif_2n) == "chr_nb.y"] <-"chr_nb2"
names(contif_2n)[names(contif_2n) == "chr_nb.z"] <-"chr_nb3"
names(contif_2n)[names(contif_2n) == "chr.x"] <-"chr1"
names(contif_2n)[names(contif_2n) == "chr.y"] <-"chr2"
names(contif_2n)[names(contif_2n) == "chr.z"] <-"chr3"
names(contif_2n)[names(contif_2n) == "ordre.x"] <-"ordre1"
names(contif_2n)[names(contif_2n) == "ordre.y"] <-"ordre2"
names(contif_2n)[names(contif_2n) == "ordre.z"] <-"ordre3" 
contif_3<-contif_2n[,c("cont","contig_id","chr_c","pos_c","orientation_c",
"chr_lg","orientation_a","cont_start","cont_end","cont_length","chr_start","chr_end","chr_length","remarques","chr_a",
"geno1","geno2","geno3","geno","size",
"size1","chr1","chr_nb1","orientation1","ordre1",
"size2","chr2","chr_nb2","orientation2","ordre2",
"size3","chr3","chr_nb3","orientation3","ordre3",
"chr_nb_consensus","ordre_consensus","orientation_consensus")]
orient_final<-vector()
for(i in 1:nrow(contif_3)){if(contif_3$orientation_consensus[i]!=0){orient_final<-c(orient_final,contif_3$orientation_consensus[i])
	}else if(contif_3$orientation_consensus[i]==0 & contif_3$orientation_a[i]==0){orient_final<-c(orient_final,contif_3$orientation_c[i])
	}else if (contif_3$orientation_consensus[i]==0 & contif_3$orientation_a[i]!=0){orient_final<-c(orient_final,contif_3$orientation_a[i])}}
contif_3$orientation_final<-orient_final
nb_cons<-vector()
for(i in 1:nrow(contif_3)){
	if(!is.na(contif_3$chr_nb_consensus[i])){nb_cons<-c(nb_cons,contif_3$chr_nb_consensus[i])
	}else{
		if(contif_3$chr_c[i]=='Un' & contif_3$chr_a[i]=='Un'){nb_cons<-c(nb_cons,'Un')
		}else if (contif_3$chr_c[i]=='Un' & contif_3$chr_a[i]!='Un'){nb_cons<-c(nb_cons,paste0('Un',contif_3$chr_a[i]))
		}else if(contif_3$chr_c[i]!='Un' & contif_3$chr_a[i]=='Un'){nb_cons<-c(nb_cons,'Un')
		}else {nb_cons<-c(nb_cons,'Un')}	
	}
}
contif_3$chr_nb_consensus<-nb_cons
contif_3$ordre_consensus[is.na(contif_3$ordre_consensus)]<-'Un'
write.table(contif_3,'combined_info.txt',sep=',',col.names=T,row.names=F,quote=F)
}else{
chr<-read.table(paste(pattern,'/chromosome/all_chr_new.txt',sep=''),header=T,sep=',',colClasses=c('character','character','character','character','character','character','character','character','character','character','character','character','character','character','character','character'))	
chr[is.na(chr)]<-0
#chr<-chr[!rowSums((is.na(chr))),]
chr<-subset(chr,chr$chromosome%in%seq(1:nbchr))
chrn<-list()
for (i in 1:length(unique(chr$chromosome))){
datac<-subset(chr,chr$chromosome==i)
datac1<-subset(datac,datac$contig%in%conti$cont[conti$chr==i])
datac2<-subset(datac,!(datac$contig%in%conti$cont[conti$chr==i]))
ord_cont1<-as.numeric(unique(datac1$contig))
contic<-subset(conti,conti$cont%in%ord_cont1 & conti$chr==i)
contico<-contic[order(contic$pos),'cont']
if(contico[1]==tail(ord_cont1,1)){datac1<-datac1[rev(rownames(datac1)),]}
ord_cont2<-as.numeric(unique(datac2$contig))
contic<-subset(conti,conti$cont%in%ord_cont2)
datacm<-datac1
if(nrow(datac2)>0){
for(k in 1:length(ord_cont2)){
		kd<-subset(datac2,datac2$contig==ord_cont2[k])
		if(nrow(kd)==1){n<-which(datacm$vector1==kd$vector1)
			if(length(n)==2){datacm<-rbind(datacm[1:n[1],],kd,datacm[n[2]:nrow(datacm),])
			}else if (length(n)==1){if(n==nrow(datacm)){datacm<-rbind(datacm,kd)}else if (n==1){datacm<-rbind(kd,datacm)}
			}else if (length(n)>2){if(nrow(datacm)%in%n){datacm<-rbind(datacm,kd)}else{datacm<-rbind(datacm[1:n[2],],kd,datacm[n[3]:nrow(datacm),])}}
		}else{
		n1<-which(datacm$vector1==kd$vector1[1])
		n2<-which(datacm$vector1==kd$vector1[nrow(kd)])
		if(length(n1)==0){if(max(n2)==nrow(datacm)){datacm<-rbind(datacm,kd[rev(rownames(kd)),])}else if (min(n2)==1){datacm<-rbind(kd[rev(rownames(kd)),],datacm)}
		}else if (length(n2)==0){if(max(n1)==nrow(datacm)){datacm<-rbind(datacm,kd)}else if (min(n1)==1){datacm<-rbind(kd,datacm)}	
		}else if (length(n1)==1&length(n2)==1){if(n2==(n1+1)){datacm<-rbind(datacm[1:n1,],kd,datacm[n2:nrow(datacm),])}else if (n1==(n2+1)){datacm<-rbind(datacm[1:n2,],kd[rev(rownames(kd)),],datacm[n1:nrow(datacm),])}
		}else if (max(n1)<min(n2)){datacm<-rbind(datacm[1:max(n1),],kd,datacm[min(n2):nrow(datacm),])}else if (max(n1)>min(n2)){datacm<-rbind(datacm[1:min(n1),],kd,datacm[(min(n1)+1):nrow(datacm),])}
		}
}}
rec<-vector()
j<-1
while(j<nrow(datacm)){
vj<-datacm$vector1[j]
vk<-datacm$vector1[j+1]
r<-0
for (k in 1:nchar(vj)){if(substr(vj,k,k)!=substr(vk,k,k)){r<-r+1}}
rec<-c(rec,r)
j=j+1}
rec<-c(0,rec)
datacm$recombination<-rec
chrn[[i]]<-datacm}
chrn<-do.call(rbind,chrn)
fwrite(chrn,paste(pattern,'/chromosome/all_chr_new2.txt',sep=''),na=NA)

pdf(paste0(pattern,'/chromosome/heatmap_colony.pdf'),width=10,height=10)
for (i in 1:length(unique(chrn$chromosome))){
chri<-subset(chrn,chrn$chromosome==i)
testchr<-do.call(rbind,strsplit(chri[,'vector1'],''))
space<-vector()
space[1]<-as.numeric(chri$size_all[1])
for(i in 1:(nrow(chri)-1)){
	space<-c(space,space[i]+as.numeric(chri$size_all[i+1]))
}
nind<-nchar(chri$vector1[1])
nmark<-nrow(testchr)
indc<-vector()
for(i in 1:nind){
	indc<-c(indc,rep(paste0('ind',i),nmark))
	}
markc<-rep(seq(1:nrow(testchr)),nind)
dat<-data.frame(indc,markc,c(testchr))
print(ggplot(dat,aes(x=dat[,1],y=dat[,2],fill=dat[,3]))+geom_tile()+scale_y_reverse())
}
dev.off()

d<-list()
for(i in 1:length(unique(chrn$chromosome))){
	sub<-subset(chrn,chrn$chromosome==i)
	di<-list()
	for(j in 1:length(unique(sub$contig))){
		sub2<-subset(sub,sub$contig==unique(sub$contig)[j])
		chromosome<-i
		contig<-unique(sub$contig)[j]
		if (sub2[1,'vector1']==sub2[nrow(sub2),'vector1']){orientation<-'0'}else if(sub2[1,'bin']<sub2[nrow(sub2),'bin']){orientation<-'+'}else if(sub2[1,'bin']>sub2[nrow(sub2),'bin']){orientation<-'-'}else if(nrow(sub2)==1){orientation<-'0'}
		nb_bins<-nrow(sub2)
		nb_markers<-sum(as.numeric(sub2$nb_markers_all))
		size<-sum(as.numeric(sub2$size_all))
		recombination<-sum(as.numeric(sub2$recombination))
		chr_cont<-paste0(chromosome,'/',contig)
		di[[j]]<-c(chromosome,contig,orientation,nb_bins,nb_markers,size,recombination,chr_cont)
		}
	d[[i]]<-do.call(rbind,di)
	}
data<-do.call(rbind,d)
colnames(data)<-c('chromosome','contig','orientation','nb_bins','nb_markers','size','recombination','chr_cont')
data<-as.data.frame(data)

for(i in 1:nrow(conti)){
	contig_n<-conti$cont[i]
	if(contig_n%in%data$contig){
	d<-subset(data,data$contig==contig_n)	
	conti[i,'chr_']<-as.numeric(as.character(d$chromosome))
	conti[i,'size_']<-as.numeric(as.character(d$size))
	conti[i,'orientation_']<-as.character(d$orientation)
	}else {
	conti[i,'chr_']<-'Un'
	conti[i,'size_']<-0
	conti[i,'orientation_']<-0}}
#fwrite(data,'contigonchromosomes.txt',na=NA)

csub<-conti[grepl('chr',colnames(conti))]
csub$chr_lg<-NULL
csub$chr_start<-NULL
csub$chr_end<-NULL
csub$chr_length<-NULL
csub$chr_c<-NULL
csub$chr_a<-NULL
for(i in 1:nrow(csub)){
	u<-unique(unname(unlist(csub[i,])))
	if(length(u)==1 & conti$orientation_[i]!=0){conti[i,'chr']<-u[1]
	}else if (length(u)==1 & conti$orientation_[i]==0 & TRUE%in%grepl('Un',u)){conti[i,'chr']<-'Un'
	}else if (length(u)==1 & conti$orientation_[i]==0 & FALSE%in%grepl('Un',u)){conti[i,'chr']<-paste0('Un',u[1])
	}else if(length(u)>1 & TRUE%in%grepl('Un',u) & conti$orientation_[i]!=0){conti[i,'chr']<-u[grepl('Un',u)==F]
	}else if(length(u)>1 & FALSE%in%grepl('Un',u) & conti$orientation_[i]==0){conti[i,'chr']<-paste0('Un',u[grepl('Un',u)==F])
	}else{conti[i,'chr']<-'Un'}}

ord<-list()
for(z in 1:nbchr){
sub_ordre<-subset(data,data$contig%in%conti$cont[conti$chr==z])
sub_conti<-subset(conti,conti$chr==z)
scontio<-sub_conti[match(unlist(as.numeric(as.character(sub_ordre$contig))),sub_conti$cont),]
if(scontio$pos[1]<scontio$pos[nrow(scontio)]){ord[[z]]<-data.frame(cont=scontio$cont,ordre=seq(1:nrow(scontio)),orientation_correct=unlist(scontio$orientation_))}
else{
	orientation_correcte<-vector()
	for(k in 1:nrow(scontio)){
		if(unlist(scontio$orientation_)[k]=='+'){orientation_correcte[k]<-'-'}else if (unlist(scontio$orientation_)[k]=='-'){orientation_correcte[k]<-'+'} else if (unlist(scontio$orientation_)[k]=='0'){orientation_correcte[k]<-'0'}
		}
	ord[[z]]<-data.frame(cont=scontio$cont,ordre=rev(seq(1:nrow(scontio))),orientation_correct=orientation_correcte)}
}
ord_c<-do.call(rbind,ord)
contif<-merge(conti,ord_c,by='cont',all=T)
contif$ordre<-as.numeric(contif$ordre)
contif$orientation_correct<-as.character(contif$orientation_correct)
contif$orientation_correct[is.na(contif$orientation_correct)]<-0

for(z in 1:nbchr){
s1<-subset(contif,contif$chr==paste0('Un',z))
if(nrow(s1)>0){
for(i in 1:nrow(s1)){
s11<-s1[i,'cont']
o<-vector()
ri<-as.numeric(rownames(data[as.numeric(as.character(data$contig))==s11,]))
if(length(ri)!=0){
ds1<-data[(ri-1):(ri+1),]
cs1<-subset(contif,contif$cont%in%ds1$contig)
if(!FALSE%in%grepl('Un',cs1$chr)){
contif$ordre[contif$cont==s11]<-NA
}else if (FALSE%in%grepl('Un',cs1$chr) & grepl(z,cs1$chr[1]) & grepl(z,cs1$chr[2]) & grepl(z,cs1$chr[3])){	
d1<-unlist(ds1[1,'contig'])
d3<-unlist(ds1[3,'contig'])
o<-c(o,paste0(cs1$ordre[cs1$cont==d1],'/',cs1$ordre[cs1$cont==d3]))}}else{if(length(o)!=0){o<-o}else{contif$ordre[contif$cont==s11]<-NA}}
if(length(unique(o))==1){contif$ordre[contif$cont==s11]<-unique(o)}else{contif$ordre[contif$cont==s11]<-NA}
}}}
contif$orientation_<-NULL
names(contif)[ncol(contif)]<-'orientation'
n<-which(grepl('/NA/',contif$ordre))
for(i in 1:length(n)){contif[n[i],'ordre']<-paste0(unlist(strsplit(contif[n[i],'ordre'],'/'))[1],'/',unlist(strsplit(contif[n[i],'ordre'],'/'))[3])}
write.table(contif,'combined_info.txt',sep=',',col.names=T,row.names=F,quote=F)

datao<-list()
for(k in 1:nrow(data)){chrk<-subset(chr,chr$contig==data$contig[k])
	if(data$orientation[k]=='-'&as.numeric(chrk$pos_m[1])<as.numeric(chrk$pos_m[nrow(chrk)])){datao[[k]]<-chrk[rev(rownames(chrk)),]
	}else if(data$orientation[k]=='-'&as.numeric(chrk$pos_m[1])>as.numeric(chrk$pos_m[nrow(chrk)])){datao[[k]]<-chrk
	}else if(data$orientation[k]=='+'&as.numeric(chrk$pos_m[1])<as.numeric(chrk$pos_m[nrow(chrk)])){datao[[k]]<-chrk
	}else if(data$orientation[k]=='+'&as.numeric(chrk$pos_m[1])>as.numeric(chrk$pos_m[nrow(chrk)])){datao[[k]]<-chrk[rev(rownames(chrk)),]
	}else if (data$orientation[k]=='0'){datao[[k]]<-chrk}}
data_order<-do.call(rbind,datao)
write.table(data_order,paste0(pattern,'/chromosome/all_chr_order.txt'),sep=',',col.names=T,row.names=F,quote=F)
}
	

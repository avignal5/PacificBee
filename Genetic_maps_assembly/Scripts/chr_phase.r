################################################
### script to phase chromosome
################################################
#!/usr/bin/env R
args<-commandArgs(TRUE)
nb_colony<-as.numeric(args[1])
pattern<-args[2]

hd<-list()
for(c in 1:nb_colony){
print(paste0('colony ',c))
filenames<-list.files(paste(pattern,c,sep=''),pattern='*.txt',full.names=T)
n<-unlist(strsplit(filenames,'/|_|.txt'))
nb<-n[seq(4,length(n),4)]
for(z in 1:length(nb)){
chrnb<-paste0('bins_',nb[z],'.txt')
f<-which(grepl(chrnb,filenames))	
chr<-read.table(filenames[f],header=T,sep=';',colClasses=c('numeric','numeric','character','character','numeric','character','numeric','character'))	
chrn<-chr
rec<-vector()
i<-1
while(i<nrow(chr)){
vtest<-chrn$vector0[i]
vj<-chr$vector0[i+1]
vk<-chr$vector1[i+1]
r1<-0
r2<-0
for (k in 1:nchar(vtest)){
if(substr(vtest,k,k)!=substr(vj,k,k)){r1<-r1+1}
if(substr(vtest,k,k)!=substr(vk,k,k)){r2<-r2+1}
}
if(r1>r2){
chrn[i+1,'bin']<-paste0(chr[i+1,'bin'],'_')
chrn[i+1,'vector0']<-chr[i+1,'vector1']
chrn[i+1,'vector1']<-chr[i+1,'vector0']}else{chrn[i+1,]<-chr[i+1,]}
rec<-c(rec,min(r1,r2))
i=i+1}
rec<-c(rec,0)
chrn$recombination<-rec
chrn$chromosome<-nb[z]
write.table(chrn,paste(pattern,c,'/chromosome/chromosome',nb[z],'.txt',sep=''),col.names=T,row.names=F,quote=F,sep=';')
hd[[z]]<-chrn
}
h_all<-do.call(rbind,hd)
write.table(h_all,paste(pattern,c,'/chromosome/all_chromosome.txt',sep=''),col.names=T,row.names=F,quote=F,sep=';')
}




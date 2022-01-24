################################################
### script to complete contig list by bins previously missing (either not genotyped for 3 colonies or without prior information on chromosome affiliation)
################################################
#!/usr/bin/env R
args<-commandArgs(TRUE)
nb_colony<-as.numeric(args[1])
contonchr<-read.table(args[2],header=T,sep=',',stringsAsFactor=F)
contonchr$contig<-as.numeric(gsub('tig','',contonchr$contig_id))
pattern<-args[3]

for(c in 1:nb_colony){
print(paste0('colony ',c))
all_chr<-read.table(paste0(pattern,c,'/chromosome/all_chromosome.txt'),header=T,sep=',',colClasses=c('character','character','character','character','numeric','numeric','numeric','numeric','character','numeric','numeric','numeric','numeric','numeric','numeric'),stringsAsFactor=F)
posc<-as.numeric(do.call(rbind,strsplit(all_chr$position_all,'_'))[,1])
posC<-as.numeric(do.call(rbind,strsplit(all_chr$position_all,'_'))[,2])
all_chr$position_all<-NULL
all_chr$pos_m<-posc
all_chr$pos_M<-posC

contmissing_c<-subset(contonchr,(contonchr[,paste0('geno',c)]=='y' & contonchr$geno=='n')|(is.na(contonchr$chr_a) & contonchr$geno=='y'))
all_bins<-read.table(paste0('bin_description_colony',c,'.txt'),header=T,sep=',',colClasses=c('numeric','numeric','numeric','numeric','numeric','character','character','character','numeric','numeric','numeric'),stringsAsFactor=F)
bins_missing<-subset(all_bins,all_bins$group%in%contmissing_c$contig)
for(i in 1:nrow(bins_missing)){bins_missing$chr[i]=subset(contmissing_c$chr,contmissing_c$contig==bins_missing$group[i])}
posm<-as.numeric(do.call(rbind,strsplit(bins_missing$position_all,'_'))[,1])
posM<-as.numeric(do.call(rbind,strsplit(bins_missing$position_all,'_'))[,2])
bins_missing$chr[is.na(bins_missing$chr)]<-'Un'
bins_missing[is.na(bins_missing)]<-1
b_miss<-data.frame(marker1=paste0(bins_missing$group,'/',bins_missing$bin,'/m'),vector1=bins_missing$vector,marker2=paste0(bins_missing$group,'/',bins_missing$bin,'_/m'),vector2=bins_missing$vector_,contig=bins_missing$group,
bin=bins_missing$bin,nb_markers_all=bins_missing$nb_markers_all,size_all=bins_missing$size_all,dist_mean=bins_missing$dist_mean,dist_median=bins_missing$dist_median,dist_var=bins_missing$dist_var,
density_markers=bins_missing$nb_markers_all/bins_missing$size_all,recombination=0,chromosome=bins_missing$chr,pos_m=posm,pos_M=posM)
b_miss$chromosome<-as.character(b_miss$chromosome)
all_chr_new<-all_chr
for(i in 1:length(unique(b_miss$vector1))){
	bsub<-subset(b_miss,b_miss$vector1==unique(b_miss$vector1)[i])
	bsubo<-bsub[order(bsub$pos_m),]
	v_b_sub<-bsubo$vector1[1]
	pos1<-which(all_chr_new$vector1==v_b_sub)
	pos2<-which(all_chr_new$vector2==v_b_sub)
	if(length(pos1)==0){
		bsubon<-bsubo
		bsubon$marker1<-bsubo$marker2
		bsubon$marker2<-bsubo$marker1
		bsubon$vector1<-bsubo$vector2
		bsubon$vector2<-bsubo$vector1
		bsubo<-bsubon	
		pos<-pos2
		}else{bsubo<-bsubo
			pos<-pos1}
	if(length(pos)==1){
			chrp<-all_chr_new[pos,'chromosome']
			for(j in 1:nrow(bsubo)){
			if(bsubo$chromosome[j]=='Un'){bsubo$chromosome[j]=all_chr_new[pos,'chromosome']}else if(bsubo$chromosome[j]==all_chr_new[pos,'chromosome']) {bsubo$chromosome[j]=all_chr_new[pos,'chromosome']} else{bsubo$chromosome[j]=paste0(bsubo[j,'chromosome'],'/',all_chr_new[pos,'chromosome'])}}
			if(all_chr_new[pos+1,'chromosome']!=all_chr_new[pos,'chromosome']){all_chr_new<-rbind(all_chr_new[1:pos,],bsubo,all_chr_new[(pos+1):nrow(all_chr_new),])
			}else if(all_chr_new[pos+1,'chromosome']==all_chr_new[pos,'chromosome']&!is.na(all_chr_new$marker1[pos-1])&all_chr_new$contig[pos-1]==all_chr_new$contig[pos]){all_chr_new<-rbind(all_chr_new[1:pos,],bsubo,all_chr_new[(pos+1):nrow(all_chr_new),])
			}else{all_chr_new<-rbind(all_chr_new[1:(pos-1),],bsubo,all_chr_new[pos:nrow(all_chr_new),])}
	}else if (length(pos)==2){
		if(all_chr_new[pos[1],'chromosome']==all_chr_new[pos[2],'chromosome']){chrp<-all_chr_new[pos[1],'chromosome']}else{chrp<-paste0(all_chr_new[pos[1],'chromosome'],'/',all_chr_new[pos[2],'chromosome'])}
		for(j in 1:nrow(bsubo)){
		if(bsubo$chromosome[j]=='Un'){bsubo$chromosome[j]=chrp}else if(bsubo$chromosome[j]==chrp) {bsubo$chromosome[j]=chrp} else{bsubo$chromosome[j]=paste0(bsubo[j,'chromosome'],'/',chrp)}}
		if(all_chr_new[pos[1],'contig']!=all_chr_new[pos[2],'contig']){all_chr_new<-rbind(all_chr_new[1:pos[1],],bsubo,all_chr_new[(pos[2]):nrow(all_chr_new),])}else {all_chr_new<-rbind(all_chr_new[1:(pos[1]-1),],bsubo,all_chr_new[pos[1]:nrow(all_chr_new),])}
	}else if (length(pos)==0){print(paste0('problem contig ',bsubo[1,'contig']))
	}else if (length(pos)>2){
		if(length(unique(all_chr_new[pos,'chromosome']))==1){chrp<-all_chr_new[pos[1],'chromosome']}else{chrp<-paste0(c(unique(all_chr_new[pos,'chromosome'])),collapse='/')}
		for(j in 1:nrow(bsubo)){
		if(bsubo$chromosome[j]=='Un'){bsubo$chromosome[j]=chrp}else if(bsubo$chromosome[j]==chrp) {bsubo$chromosome[j]=chrp} else{bsubo$chromosome[j]=paste0(bsubo[j,'chromosome'],'/',chrp)}}
		if(length(unique(all_chr_new[pos,'chromosome']))==1 & all_chr_new$contig[min(pos)-1]%in%all_chr_new$contig[pos] & all_chr_new$contig[max(pos)-1]%in%all_chr_new$contig[pos]){all_chr_new<-rbind(all_chr_new[1:(max(pos)-1),],bsubo,all_chr_new[max(pos):nrow(all_chr_new),])
		} else if (length(unique(all_chr_new[pos,'chromosome']))==1 & !(all_chr_new$contig[min(pos)-1]%in%all_chr_new$contig[pos]) & all_chr_new$contig[max(pos)-1]%in%all_chr_new$contig[pos]){all_chr_new<-rbind(all_chr_new[1:(max(pos)-1),],bsubo,all_chr_new[max(pos):nrow(all_chr_new),])
		}else if (length(unique(all_chr_new[pos,'chromosome']))==1 & all_chr_new$contig[min(pos)-1]%in%all_chr_new$contig[pos] & !(all_chr_new$contig[max(pos)-1]%in%all_chr_new$contig[pos])){all_chr_new<-rbind(all_chr_new[1:max(pos),],bsubo,all_chr_new[(max(pos)+1):nrow(all_chr_new),])
		}else{print(paste0('problem contig ',bsubo$contig))}}
}
all_chr_new[is.na(all_chr_new$marker1),'recombination']<-NA
all_chr_new[is.na(all_chr_new$marker1),'chromosome']<-NA
write.table(all_chr_new,paste0(pattern,c,'/chromosome/all_chr_new.txt'),col.names=T,row.names=F,quote=F,sep=',')

}



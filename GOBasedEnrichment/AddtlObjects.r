#A file describing two different objects used in the script 'Define_ASDmeQTLtargets.r'

#############################################
#############################################
#####Object 1: Recombination hotspot file####
#############################################
#############################################

#Downloaded from: ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/pilot_data/technical/reference/ (genetic_map_b36.tar.gz)
map36<-read.table("/amber2/scratch/claddaco/SEEDmethylation/Gaphunting/dbSNP138/Map36/hotspots_b36.txt",header=TRUE)

#Format a table for input to LiftOver software (https://genome.ucsc.edu/cgi-bin/hgLiftOver) and send in. 
tablines<-c()
for (i in 1:32991){
	templine<-paste(map36$Chromosome[i],":",map36$Start[i],"-",map36$End[i],sep="")
	tablines<-rbind(tablines,templine)
}
write.table(tablines,file="/amber2/scratch/claddaco/SEEDmethylation/Gaphunting/dbSNP138/Map36/tablines.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

#Read in LiftOver output.  
map37<-read.table(file="/amber2/scratch/claddaco/SEEDmethylation/Gaphunting/dbSNP138/Map36/hglft_genome_17ef_756770.bed",header=FALSE)

#Format output, define everything between recombination hot spots as LD blocks/haplotypes. 
#Will generate a isRHS column, which will be 0/1 indicator variable.
	#0: is RHS
	#1: is LD block
map37RHS<-data.frame(0,nrow(map37),3)
for (i in 1:nrow(map37)){
#for (i in 1:1){
	map37RHS[i,3]<-unlist(strsplit(as.character(map37[i,1]),"-"))[2]
	temper<-unlist(strsplit(as.character(map37[i,1]),"-"))[1]
	map37RHS[i,1:2]<-unlist(strsplit(temper,":"))
}
colnames(map37RHS)<-c("CHR","Start","Stop")	

chrnum<-as.numeric(ifelse(nchar(map37RHS$CHR)==5,substr(map37RHS$CHR,4,5),substr(map37RHS$CHR,4,4)))
map37RHS.order<-cbind(map37RHS,chrnum)
map37RHS.order.fin<-map37RHS.order[order(as.numeric(map37RHS.order$chrnum),as.numeric(map37RHS.order$Start)),]

RRspots37<-data.frame(0,(2*nrow(map37RHS.order.fin)-1),4)
keptsequence<-seq(1,(2*nrow(map37RHS.order.fin)-1),by=2)
RRspots37[keptsequence,1:3]<-c(map37RHS.order.fin[,1:3])
RRspots37[keptsequence,4]<-1
colnames(RRspots37)<-c("CHR","Start","Stop","isRHS")	
changesequence<-seq(2,(2*nrow(map37RHS)-1),by=2)
for (k in changesequence){
	RRspots37[k,1]<-RRspots37[(k-1),1]
	RRspots37[k,2]<-RRspots37[(k-1),3]
	RRspots37[k,3]<-RRspots37[(k+1),2]
	RRspots37[k,4]<-0
}
corRemlist<-which(as.numeric(RRspots37$Start)>as.numeric(RRspots37$Stop))
identical(corRemlist,remlist)

RRspots37<-RRspots37[-corRemlist,]
save(RRspots37,file="RRspots37.rda")

#############################################
#############################################
#####Object 2: RefSeq Gene File##############
#############################################
#############################################

#We downloaded the complete list of RefSeq genes from the UCSC genome browser (http://genome.ucsc.edu/cgi-bin/hgTables - refFlat table)
refFlat<-read.table("refFlat_March2016.txt",header=FALSE,stringsAsFactors=FALSE)
#Now clean the refFlat object. 
refFlat.new<-refFlat[,c(1,3,5:8)]
validchr<-paste0("chr",1:22)
refFlat.new<-refFlat.new[which(refFlat.new[,2]%in%validchr),]
refFlat.clean<-data.frame(matrix(0,length(unique(refFlat.new[,1])),6))
for (i in 1:length(unique(refFlat.new[,1]))){
#for (i in 1:2){
	truncated<-refFlat.new[which(refFlat.new[,1]==unique(refFlat.new[,1])[i]),]
	if (nrow(truncated)==1){refFlat.clean[i,]<-truncated[1,]}
	if (nrow(truncated)>1){
		refFlat.clean[i,1:2]<-(truncated[1,1:2])
		refFlat.clean[i,3]<-min(truncated[,3])
		refFlat.clean[i,4]<-max(truncated[,4])
		refFlat.clean[i,5]<-min(truncated[,5])
		refFlat.clean[i,6]<-max(truncated[,6])
	}
}
colnames(refFlat.clean)<-c("Gene","CHR","txStart","txEnd","cdsStart","cdsEnd")
#refFlat.clean<-refFlat.clean[-which(refFlat.clean[,2]=="chrX"|refFlat.clean[,2]=="chrY"),]
#refFlat.clean<-refFlat.clean[-grep("chrUn",refFlat.clean$CHR,fixed=T),]
refFlat.clean<-refFlat.clean[-grep("MIR",refFlat.clean$Gene,fixed=T),] #Remove microRNAs
#FURTHER CLEAN THIS GENE LIST TO GET RID OF STUFF
#Long intergeneic non-coding RNA: "LINC"
#Long non-coding RNA: "LOC"
#small associated RNA: "SNAR"
refFlat.clean<-refFlat.clean[-grep("LINC",refFlat.clean$Gene,fixed=T),]
refFlat.clean<-refFlat.clean[-grep("LOC",refFlat.clean$Gene,fixed=T),]
refFlat.clean<-refFlat.clean[-grep("SNAR",refFlat.clean$Gene,fixed=T),]

chrFlat<-substring(refFlat.clean[,2],4,5)
chrFlat<-ifelse(chrFlat=="6_",6,chrFlat)
chrFlat<-ifelse(chrFlat=="1_",1,chrFlat)
chrFlat<-ifelse(chrFlat=="4_",4,chrFlat)
chrFlat<-ifelse(chrFlat=="7_",7,chrFlat)
chrFlat<-as.numeric(chrFlat)
refFlat.clean.new<-data.frame(chrFlat,refFlat.clean[,c(1,3:6)])
colnames(refFlat.clean.new)<-c("CHR","Gene","txStart","txEnd","cdsStart","cdsEnd")
refFlat.clean.new.order<-refFlat.clean.new[order(refFlat.clean.new$CHR,refFlat.clean.new$txStart),]
save(refFlat.clean.new.order,file="refFlat.clean.new.order.rda")

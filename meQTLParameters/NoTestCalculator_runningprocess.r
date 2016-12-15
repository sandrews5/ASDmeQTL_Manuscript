#This script will walk through how to set up input for the function 
#'NoTestCalculator', which was used to calculate the total number of SNP-CpG
#comparisons to be made given certain cutoffs. 

#STEP 1: LOAD DNAm INFORMATION
library(minfi)
load("object.qnorm_postQC.Batch1.rda") #Of class GenomicRatioSet.
object<-object.postQC
pd<-pData(object)

#STEP 1a: Limit to white samples as determined by SNP PCs (see Methods)
seed1<-read.csv("/amber2/scratch/claddaco/SEEDmethylation/SampleAnnotation/MasterAnnotationSEEDmethylation.csv",header=TRUE)
whites<-seed1[which(seed1$Race.y=="white"),]
object.whites<-object[,which(pd$Family%in%whites$Family)] #just do this because this object was cleaned out for duplicates alread
pd.whites<-pd[which(pd$Family%in%whites$Family),]

#STEP 2: Generate the following vectors for SNP information: MAF, CHR, position
library(sqldf)
totalmaflist<-c()
totalsnppos<-c()
totalsnpchr<-c()
for (chr in 1:22){
	print(chr)
	f<-file(paste("SEED_MethGeno_",chr,".tped",sep="")) 
	SNPgeno <- sqldf("select * from f", dbname = tempfile(), file.format = list(header = F, row.names = F, sep=" "))
	closeAllConnections() 
	peoplemap<-read.table(paste("SEED_MethGeno_",chr,".tfam",sep=""))
	
	seqskip1<-seq(5,ncol(SNPgeno),by=2)
	seqskip2<-seq(6,ncol(SNPgeno),by=2)
	translated<-t((SNPgeno[,seqskip1])+(SNPgeno[,seqskip2]))
	reformedgeno<-ifelse(translated==2,0,9)
	reformedgeno<-ifelse(translated==3,1,reformedgeno)
	reformedgeno<-ifelse(translated==4,2,reformedgeno)
	rm(translated)
	
	#Match to order of samples in methylation data. Not really necessary for these purposes but will use
	#same code in places that it is. 
	reformedgeno<-reformedgeno[match(pd.whites$Family,peoplemap[,2]),]
	
	#Ensure that dosage is count of minor alleles and MAF is always <= 0.5
	#Collect and retain necessary information
	maflist<-rep(0,ncol(reformedgeno))
	reformedgeno.recode<-reformedgeno
	for (what in 1:ncol(reformedgeno)){
		zeros<-length(which(reformedgeno[,what]==0))
		ones<-length(which(reformedgeno[,what]==1))
		twos<-length(which(reformedgeno[,what]==2))
		total<-nrow(reformedgeno)-length(which(reformedgeno[,what]==9))
		p<-((2*zeros)+ones)/(2*total)
		if (p<(1-p)){reformedgeno.recode[,what]<-ifelse(reformedgeno.recode[,what]==0,2,reformedgeno.recode[,what])
				reformedgeno.recode[,what]<-ifelse(reformedgeno[,what]==2,0,reformedgeno.recode[,what])
			maflist[what]<-p
		} else {maflist[what]<-(1-p)}
	}
	totalmaflist<-c(totalmaflist,maflist)
	totalsnppos<-c(totalsnppos,SNPgeno[,4])
	totalsnpchr<-c(totalsnpchr,rep(paste("chr",chr,sep=""),ncol(reformedgeno.recode)))
}
#STEP 2a: Save objects. 
save(totalmaflist,file="totalmaflist.rda")
save(totalsnppos,file="totalsnppos.rda")
save(totalsnpchr,file="totalsnpchr.rda")

#STEP 3: Load SD, pos, and chr vectors for DNAm information 
chrnames<-as.character(seqnames(object.whites))
object.whites<-object.whites[-which(chrnames=="chrX"|chrnames=="chrY"),]
B<-getBeta(object.whites)
rowsds<-apply(B,1,sd)
rowmeans<-apply(B,1,mean)
chrnames<-as.character(seqnames(object.whites))
pos<-as.numeric(start(object.whites))

#STEP 3a: Calculate sd for all probes
sdvector=seq(0.15,0.30,by=0.05)
sdcutoffs<-quantile(rowsds,unique(sdvector))

#STEP 3b: Plot historgrams of sd and means for all probes to better inform values input in Quanto power calculator. 
pdf("SEED_MeansandSD.pdf")
for (i in 1:length(sdvector)){
	par(mfrow=c(2,2))
	tempsds<-rowsds[which(rowsds>sdcutoffs[i])]
	tempmean<-rowmeans[which(rowsds>sdcutoffs[i])]
	lowmeth<-which(tempmean<0.4)
	highmeth<-which(tempmean>0.7)
	
	lowsds<-tempsds[lowmeth]
	lowmeans<-tempmean[lowmeth]
	highsds<-tempsds[highmeth]
	highmeans<-tempmean[highmeth]
	
	hist(lowsds,breaks=25,main=paste("Mean ",round(mean(lowsds),digits=2),sep=""))
	abline(v=mean(lowsds),col="red")
	hist(lowmeans,breaks=25,main=paste("Mean ",round(mean(lowmeans),digits=2),sep=""))
	abline(v=mean(lowmeans),col="red")
	hist(highsds,breaks=25,main=paste("Mean ",round(mean(highsds),digits=2),sep=""))
	abline(v=mean(highsds),col="red")
	hist(highmeans,breaks=25,main=paste("Mean ",round(mean(highmeans),digits=2),sep=""))	
	abline(v=mean(highmeans),col="red")
}
dev.off()

#STEP 4: Set up parallelization.
library(doParallel)
library(doRNG)
registerDoParallel(cores=10)	
verbose=TRUE
if (!getDoParRegistered())
	registerDoSEQ()
workers <- getDoParWorkers()
backend <- getDoParName()
version <- getDoParVersion()

#STEP 5: Run function with generated vectors. 
sdvector=0.15
MAFvector=seq(0.02,0.04,0.0025)
windowsizevector=c(500000,1000000,2000000,3000000,4000000)
testresults<-NoTestCalculator(chrnames,pos,rowsds,totalsnpchr,totalsnppos,totalmaflist,sdvector=sdvector,MAFvector=MAFvector,windowsizevector=windowsizevector)
write.csv(testresults,file="/amber2/scratch/claddaco/SEEDmethylation/meQTLPaper/PowerCalculations/SEED_numtestresults1redo.csv",col.names=TRUE)

#STEP 6: Run parameters in Quanto to calculate power. See 'QuantoTemplate.qpp'
#Perform the meQTL analysis. 

args<-commandArgs(trailingOnly=TRUE) #Send in 1 script per chromosome via $SGE_TASK_ID
chr<-as.numeric(args[1])

#Define the cutoffs as calculated via the scripts in meQTLParameters folder
MAFcutoff<-0.0275
sdcutoff<-0.15
windowsize<-1000000
pvaluethreshold<-1E-4 #Printing threshold for matrixEQTL, see below. 
resultsdirectory<-"" #Where to put matrixEQTL output objects
genotypedirectory<-"" #Where to grab the genotype data from

library(minfi)
#Load methylation object and phenotype data. 
load("object.qnorm_postQC.Batch1.rda")
object<-object.postQC
pd<-pData(object)

#Limit methylation data to white individuals (see Methods).  
seed1<-read.csv("MasterAnnotationSEEDmethylation.csv",header=TRUE)
whites<-seed1[which(seed1$Race.y=="white"),]
object.whites<-object[,which(pd$Family%in%whites$Family)] #just do this because this object was cleaned out for duplicates alread
pd.whites<-pd[which(pd$Family%in%whites$Family),]

#Calculate beta values and only include those that pass methylation sd threshold 
#defined earlier and calculated via scripts in meQTLParameters folder.
#Also remove X and Y chromosomes. 
chrnames<-as.character(seqnames(object.whites))
object.whites<-object.whites[-which(chrnames=="chrX"|chrnames=="chrY"),]
B<-getBeta(object.whites)
rowsds<-apply(B,1,sd)
chrnames<-as.character(seqnames(object.whites))
pos<-as.numeric(start(object.whites))
sdcutoffs<-quantile(rowsds,sdcutoff)
B<-B[-which(rowsds<sdcutoffs),]
chrnames<-chrnames[-which(rowsds<sdcutoffs)]
pos<-pos[-which(rowsds<sdcutoffs)]

#Read in race PCs and format them as numeric variables. 
firstpheno<-read.csv("SEED_Full_Phenotype_Smoke_Alc_B2R_SSRI_Annotation_20150416.csv",header=TRUE)
racepcs<-firstpheno[,5:14]
rownames(racepcs)<-firstpheno$FamilyID
racepcs<-racepcs[match(pd.whites$Family,firstpheno$FamilyID),]
racepcs<-lapply(1:ncol(racepcs),function(x){(as.character(racepcs[,x]))})
racepcs<-t(do.call("rbind",racepcs))
getrid<-racepcs=="no_genetic_ancestry_data" #Anything missing? No, good to go. 
racepcs<-lapply(1:ncol(racepcs),function(x){(as.numeric(racepcs[,x]))})
racepcs<-t(do.call("rbind",racepcs))

#Define the model matrix to be used in meQTL regressions. 
mod=model.matrix(~as.factor(pd.whites$predictedSex)+racepcs[,1:2]) #Covariates: First 2 PCs and sex.

#Compute and extract residuals. We will use residuals as outcome in matrixEQTL invocation. 
library(limma)
fit<-lmFit(B,mod)
resids<-residuals(fit,B)
regressionobj<-list("resids"=resids,"chrnames"=chrnames,"pos"=pos)
save(regressionobj,file=paste(resultsdirectory,"regressionobj.rda",sep="")) #Save once and load again for use on all chromosomes. 

load(paste(resultsdirectory,"regressionobj.rda",sep=""))
resids<-regressionobj$resids
chrnames<-regressionobj$chrnames
pos<-regressionobj$pos

#Read in genotype data in .tpd format
library(sqldf)
f<-file(paste(genotypedirectory,"SEED_MethGeno_",chr,".tped",sep=""))
SNPgeno <- sqldf("select * from f", dbname = tempfile(), file.format = list(header = F, row.names = F, sep=" "))
closeAllConnections() 
peoplemap<-read.table(paste(genotypedirectory,"SEED_MethGeno_",chr,".tfam",sep="")) #CHANGE

#Format genotype in terms of dosage of minor allele. Reorder samples so that they match
#the same order in the DNAm data.
seqskip1<-seq(5,ncol(SNPgeno),by=2)
seqskip2<-seq(6,ncol(SNPgeno),by=2)
translated<-t((SNPgeno[,seqskip1])+(SNPgeno[,seqskip2]))
reformedgeno<-ifelse(translated==2,0,9)
reformedgeno<-ifelse(translated==3,1,reformedgeno)
reformedgeno<-ifelse(translated==4,2,reformedgeno)
rm(translated)

reformedgeno<-reformedgeno[match(pd.whites$Family,peoplemap[,2]),]

#Ensure that dosage is count of minor alleles
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

#Limit SNPs to those greater than the input MAF cutoff
reformedgeno.recode<-reformedgeno.recode[,which(maflist>MAFcutoff)] 

#Ordering around genotypes by increasing position in case it's necessary, and 
#collect position information for SNPs to be used in matrixEQTL. 
allsnps.ord<-t(reformedgeno.recode)
allpos<-SNPgeno[which(maflist>MAFcutoff),4]
allsnps.ord<-allsnps.ord[order(allpos),]
allpos.ord<-allpos[order(allpos)]
allsnps.ord<-ifelse(allsnps.ord==9,NA,allsnps.ord)

#Remove some objects to free up some memory. 
rm(SNPgeno,reformedgeno,reformedgeno.recode)
rownames(allsnps.ord)<-allpos.ord

#Load matrixEQTL and format SNPs and methylation data in appropriate format to be used by the package. 
library(MatrixEQTL)
usesnps<-SlicedData$new(allsnps.ord) #format SNPs in sliced data class the function needs
mysnppos<-data.frame(allpos.ord,rep(chr,nrow(allsnps.ord)),allpos.ord)
colnames(mysnppos)<-c("snpid","chr","pos")
relevantmeth<-resids[which(chrnames==paste("chr",chr,sep="")),] #probes on this chr
relevantmethpos<-data.frame(rownames(relevantmeth),rep(chr,nrow(relevantmeth)),pos[which(chrnames==paste("chr",chr,sep=""))],pos[which(chrnames==paste("chr",chr,sep=""))])
colnames(relevantmethpos)<-c("geneid","chr","left","right")

#The object breakpoint well send in only "breakpoint" DNAm probes at a time to be used in matrixEQTL engine. 
#It will output results according to indices creating by the number of sections needed to complete the total number of probes
#on this chromosome in 'breakpoint' increments. 
breakpoint<-2000
breaks<-ceiling(length(which(chrnames==paste("chr",chr,sep="")))/breakpoint) #breaking up probes into breakpoint group segments

#Remove all uncessary objects to free up more memory. 
rm(list=ls()[-(which(ls()%in%c("pvaluethreshold","resultsdirectory","windowsize","usesnps","relevantmeth","breaks","chr","relevantmethpos","mysnppos","breakpoint")==TRUE))])

for (j in 1:breaks){
	if (j>1){rm(ret)} #again get rid of previous 'ret' objects for memory purposes
	#Perform indexing of DNAm data differently for last section of object. 
		if (j<breaks){
		usemeth<-SlicedData$new(relevantmeth[(((breakpoint*j)-(breakpoint-1)):(j*breakpoint)),]) #different indexing of probes
		usepos<-relevantmethpos[(((breakpoint*j)-(breakpoint-1)):(j*breakpoint)),]
		} else {
			usemeth<-SlicedData$new(relevantmeth[(((breakpoint*j)-(breakpoint-1)):nrow(relevantmeth)),])
			usepos<-relevantmethpos[(((breakpoint*j)-(breakpoint-1)):nrow(relevantmeth)),]
		}
	#Invoke matrixEQTL function.
	system.time(ret<-Matrix_eQTL_main(usesnps,usemeth,pvOutputThreshold=0,pvOutputThreshold.cis=pvaluethreshold,
		 cisDist=windowsize,snpspos=mysnppos,genepos=usepos,output_file_name=NULL,output_file_name.cis=NULL))
	#Retain and save relevant output. 
	ret<-ret$cis$eqtls[,c("snps","gene","pvalue","beta")]
	save(ret,file=paste(resultsdirectory,"restotal",chr,"_",j,".rda",sep=""))
}
rm(list=ls())
q()
n

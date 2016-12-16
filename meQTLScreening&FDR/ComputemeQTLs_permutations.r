#Generate the null object outputs to be used in FDR estimation. 
#See methods and Shi 2015 for description of approach. 

args<-commandArgs(trailingOnly=TRUE) #Send in 1 script per chromosome via $SGE_TASK_ID
chr<-as.numeric(args[1])

#Define the cutoffs as calculated via the scripts in meQTLParameters folder
MAFcutoff<-0.0275
sdcutoff<-0.15
windowsize<-1000000
pvaluethreshold<-1E-4 #Printing threshold for matrixEQTL, see below. 
resultsdirectory<-"" #Where to put matrixEQTL output objects
permresultsdirectory<-"" #Where to put permutation results. 
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

#Load residuals object generated in ComputemeQTLs.r
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

#The object "myscrambles" should be generated once, and loaded each time for 
#every chromsome. That way for the nth permutation, the scrambled order of samples for chromosome I 
#is the same as the scrambled order for chromosome J, even though they are being run in different scripts
	# mypermutes<-matrix(rep(seq(1,ncol(allsnps.ord)),100),nrow=100,ncol=(ncol(allsnps.ord)),byrow=TRUE)
	# myscrambles<-apply(mypermutes,1,sample)
	# save(myscrambles,file=paste(permresultsdirectory,"myscrambles.rda",sep=""))

load(paste(permresultsdirectory,"myscrambles.rda",sep=""))

#Set up parallelization. 
library(doParallel)
library(doRNG)
registerDoParallel(cores=10)	
verbose=TRUE
if (!getDoParRegistered())
	registerDoSEQ()
workers <- getDoParWorkers()
backend <- getDoParName()
version <- getDoParVersion()
chunksize <- ceiling(ncol(myscrambles)/workers) 

#Run the meQTL query for each permutation. 
tmp <- foreach(probes = iter(myscrambles, by = "column", chunksize = chunksize)) %dorng%
	{apply(probes,2,function(k){
		#Load matrixEQTL and format SNPs and methylation data in appropriate format to be used by the package. 	
		library(MatrixEQTL)
		usesnps<-SlicedData$new(allsnps.ord[,k]) #format SNPs in sliced data class the function eneds
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

		#Initialize an object that will collect the significnat cpg sites and their p-values; this is all that we need for FDR calculation. 
		runningnames<-matrix(0,0,2)
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
				ret<-Matrix_eQTL_main(usesnps,usemeth,pvOutputThreshold=0,pvOutputThreshold.cis=pvaluethreshold,
						cisDist=windowsize,snpspos=mysnppos,genepos=usepos,output_file_name=NULL,output_file_name.cis=NULL)
				ret<-ret$cis$eqtls[,c("gene","pvalue")]
				runningnames<-rbind(runningnames,ret)
			}
		return((runningnames))
		})
	}
#Collect the significant cpg sites in the null sets but be sure to keep information on which 
#turned up in which permutation. Save. 
cgnames<-unlist(tmp,recursive=FALSE)	
save(cgnames,file=paste(permresultsdirectory,"cgnames_",chr,".rda",sep=""))
q()
n
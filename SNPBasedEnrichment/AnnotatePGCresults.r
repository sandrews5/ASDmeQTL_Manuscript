#Here we will do the following:
	#Create objects of the PGC lists that survive pruning (i.e. their direct overlap or proxy survived pruning)
	#Attach several columns to them
		#is asd p-value less than 1E-3,1E-4,1E-5 (3 columns)
		#is meQTL at FDR of 10%, 5%, 1%
		#Count of CpGs within window size

pgcresults<- #Location of 'truncres' objects generated in 'PGCoverlap.r'. 
myfiles<-list.files(path="") #Where meQTL results live. 		
pvaluecutoffs<-c(3.1E-5,1E-5,3E-7) #p-value cutoffs for FDR 10%, 5%, 1%

#Parameters used in meQTL query. 
MAFcutoff<-0.0275 ###
sdcutoff<-0.15
windowsize<-1000000 ###

#Load objects generated in power calculations. 
load("totalmaflist.rda")
load("totalsnppos.rda")
load("totalsnpchr.rda")

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

asdenrichmentdirectory<-"" #Where to put objects relevant to performing SNP based enrichment 
pruningdirectory<-"" #Location of Priority pruner input tables, as well as the pruned genotype files generated next in 'RunPriorityPruner.sh'
overlapdirectory<-"" #DIRECTORY OF SNP FILES THAT JUST HAVE THOSE SNPS THAT OVERLAP PGC; set in 'PriorityPruner_MakeGenotypes.sh'


for (chr in 1:22){
	#Load the 'truncres' objects and identify those snps that survive pruning.
	#'prunein' object will contain only those SNPs
	load(paste(pgcresults,"truncres_",chr,".rda",sep=""))
	truncres$matchme<-as.character(truncres$SNP)
	truncres$matchme[which(truncres$indicator==1)]<-truncres$Proxy[which(truncres$indicator==1)]
	library(sqldf)
	f<-file(paste(overlapdirectory,"SEED_Overlap_",chr,".tped",sep=""))
	SNPgeno <- sqldf("select * from f", dbname = tempfile(), file.format = list(header = F, row.names = F, sep=" "))
	closeAllConnections() 
	SNPgeno<-SNPgeno[,c(2,4)]
	truncres$matchpos<-SNPgeno[match(truncres$matchme,SNPgeno[,1]),2]
	pruneres.asd<-read.table(paste(pruningdirectory,"SEED_Pruned_",chr,".results",sep=""),header=TRUE) 
	kept<-pruneres.asd[which(pruneres.asd$selected==1),]
	prunein<-truncres[which(truncres$matchpos%in%kept$pos),]
	
	#Annotate with ASD status (0/1 variable)
	prunein$isASD3<-ifelse(prunein$P<=1E-3,1,0)
	prunein$isASD4<-ifelse(prunein$P<=1E-4,1,0)
	
	#Annotate with meQTL status based on 3 FDR-based cutoffs. 
	thischr<-grep(paste("restotal",chr,"_",sep=""),myfiles,fixed=TRUE) #or grep for 'ret'
	myres<-lapply(1:length(thischr),function(x){load(paste("",myfiles[thischr[x]],sep="")) #Add root file path for meQTL results. 
		ret1<-ret[which(ret$pvalue<pvaluecutoffs[1]),] 
		ret2<-ret[which(ret$pvalue<pvaluecutoffs[2]),] 
		ret3<-ret[which(ret$pvalue<pvaluecutoffs[3]),] 
		return(list(cut1=prunein$matchpos%in%ret1$snps,cut2=prunein$matchpos%in%ret2$snps,cut3=prunein$matchpos%in%ret3$snps))
	})
	cut1<-lapply(myres,function(x){return(x$cut1)})
	prunein$ismeQTL10<-rowSums(do.call("cbind",cut1))
	cut2<-lapply(myres,function(x){return(x$cut2)})
	prunein$ismeQTL5<-rowSums(do.call("cbind",cut2))
	cut3<-lapply(myres,function(x){return(x$cut3)})
	prunein$ismeQTL1<-rowSums(do.call("cbind",cut3))

	#Now get information on study-based MAF, and count of CpGs in the window. 	
	snpmaf<-totalmaflist[which(totalsnpchr==paste("chr",chr,sep=""))]
	snppos<-totalsnppos[which(totalsnpchr==paste("chr",chr,sep=""))]
	methpos<-pos[which(chrnames==paste("chr",chr,sep=""))]
	
	snpmaf<-snpmaf[match(prunein$matchpos,snppos)]
	range(snpmaf) #CHECK THIS
	snppos<-snppos[match(prunein$matchpos,snppos)] #IDENTICAL TO prunein$matchpos
	identical(snppos,prunein$matchpos)
	
	#Count the number of 450k CpGs within the meQTL detection limit of each SNP. 
	#We match on this later to account for differential opporuntity for being a meQTL that proximity to CpG sites introduces. 
	breakpoint<-500
	breaks<-ceiling(nrow(prunein)/breakpoint) #breaking up probes into breakpoint group segments
	runningsum<-c()
	for (j in 1:breaks){
		message("chr ",chr," on break ",j," out of ",breaks,sep="")
			if (j<breaks){
			usesnp<-snppos[(((breakpoint*j)-(breakpoint-1)):(j*breakpoint))] 
			} else {
				usesnp<-snppos[(((breakpoint*j)-(breakpoint-1)):length(snppos))]
			}

		res<-abs(outer(usesnp,methpos,"-"))
		res<-res<=windowsize
		runningsum<-c(runningsum,rowSums(res))
	}
	prunein$CountCGs<-runningsum	
	prunein$StudyMAF<-snpmaf
	save(prunein,file=paste(asdenrichmentdirectory,"prunein_asd_",chr,".rda",sep=""))
}			

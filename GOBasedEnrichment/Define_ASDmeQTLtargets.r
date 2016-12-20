#In this script, we define the meQTL targets of ASD-related SNPs and 
#their proxies. There are 2 main parts to this:
		#1) Identifying ASD SNPs in PGC
		#2) Identifiying which ASD SNPs and their proxies are meausred/imputed in respective study sets
		#3) Identifying the meQTL targets of those SNPs. 

#####################################################################
#####################################################################
###########################PART 1####################################
#####################################################################
#####################################################################

#Read in PGC data. Define ASD-related SNPs via p-value threshold of 1E-4. 
pvaluethresh<-1E-4
library(sqldf)
f<-file("/legacy/amber2/scratch/claddaco/seedgewis/PGC_ASD/onqah.pgcasdeuro")
pgcnew <- sqldf("select * from f", dbname = tempfile(), file.format = list(header = T, row.names = F, sep="\t"))

collect<-lapply(seq(1:22),function(chr){
	thispgc<-pgcnew[which(pgcnew$CHR==chr),]
	truncres<-thispgc[,c("SNP","BP","EUR_FRQ","OR","P")]
	colnames(truncres)<-c("SNP","BP","MAF","OR","P")
	truncres$CHR<-chr
	
	truncres<-truncres[which(truncres$P<pvaluethresh),]
	return(truncres)
	}
)
collect<-do.call("rbind",collect)
save(collect,file="collect.rda")

load("collect.rda")
#Load recombination hotspot file. See 'AddtlObjects.r'
load("RRspots37.rda")
#Load RefSeq gene list. See 'AddtlObjects.r'
load("refFlat.clean.new.order.rda")

#We will create this 'collectplus' object which will delimit a variety of important information for each
#ASD-related (PGC p-value < 1E-4) SNP. This information includes: 
	#Genes that harbor these SNPs (if any)
	#Indicator variables to see if that SNP is an meQTL
#We only need the SNP ID info for purposes of the GO analysis, but we will use this collectplus object
#later on in the GWAS loci expansion analysis. 

#Load the fetal brain meQTL results and limit to similar window size as before. 
windowsize<-1000000
fetalbrain<-read.csv("All_Imputed_BonfSignificant_mQTLs.csv",header=TRUE)
fetalbrain<-fetalbrain[which(fetalbrain$SNP_Chr==fetalbrain$DNAm_CHR),]
fetalbrain$Space<-abs(fetalbrain$SNP_BP-fetalbrain$DNAm_BP)
fetalbrain<-fetalbrain[which(fetalbrain$Space<=windowsize),]

#Make the collectplus object. Collect genes that these SNPs map to using GenomicRanges::findOverlaps, 
#and create meQTL indicator variables. 
library(GenomicRanges)
asdgenelist<-list() #Initialize gene list
meQTLfilelist<-c("/SEED_meQTLs/","/EARLI_meQTLs/","fetalbrain")
templabels<-c("BloodmeQTL","CBmeQTL","BrainmeQTL")
pvaluecutoffs<-c(1E-5,2.7E-6,1E-8)
collectplus<-matrix(0,0,13)
for (chr in 1:22){
	print(chr)
	temp<-collect[which(collect$CHR==chr),]
	thischr<-RRspots37[which(RRspots37$CHR==paste("chr",chr,sep="")),]
	gr.RHS<-GRanges(seqnames = Rle(rep(chr,nrow(thischr))), ranges = IRanges(start=as.numeric(thischr$Start), end=as.numeric(thischr$Stop)))
	gr.temp<-GRanges(seqnames = Rle(rep(chr,nrow(temp))), ranges = IRanges(start=as.numeric(temp$BP), end=as.numeric(temp$BP)))		
	ovl<-findOverlaps(gr.RHS,gr.temp)
	temp$Assignment<-as.numeric(as.factor(queryHits(ovl)))
	temp$isRHS<-as.numeric(thischr$isRHS[queryHits(ovl)])
	temp$BlockStart<-as.numeric(thischr$Start[queryHits(ovl)])
	temp$BlockStop<-as.numeric(thischr$Stop[queryHits(ovl)])
	genelist<-refFlat.clean.new.order[which(refFlat.clean.new.order$CHR==chr),]
	gr.genes<-GRanges(seqnames = Rle(rep(chr,nrow(genelist))), ranges = IRanges(start=as.numeric(genelist$txStart), end=as.numeric(genelist$txEnd)))		
	ovl<-findOverlaps(gr.genes,gr.temp)
	
	#Attach gene information. Find the genes implicated by each peak. 
	assigngenes<-lapply(unique(temp$Assignment),function(x){
		peaksnps<-which(temp$Assignment==x)
		relevanthits<-which(subjectHits(ovl)%in%peaksnps)
		if(length(relevanthits)>0){
		return(unique(genelist$Gene[queryHits(ovl)[relevanthits]]))
		} else{return(NA)}				
	})
	asdgenelist[[chr]]<-assigngenes
	
	#See if each SNP has a direct overlap meQTL in blood, cord blood and brains. 
	for (y in 1:3){
		if(y<3){
			thischr<-grep(paste("restotal",chr,"_",sep=""),list.files(meQTLfilelist[y]),fixed=TRUE)
			myres<-lapply(1:length(thischr),function(x){load(paste(meQTLfilelist[y],list.files(meQTLfilelist[y])[thischr[x]],sep="")) #CHANGE HERE
				ret<-ret[which(ret$pvalue<pvaluecutoffs[y]),] #FDR 5%
				return(temp$BP%in%ret$snps)})
			myres<-do.call("cbind",myres)
			temp$New<-rowSums(myres)
			colnames(temp)[which(colnames(temp)=="New")]<-templabels[y]
		} else{
			thischr<-fetalbrain[which(fetalbrain$SNP_Chr==chr),]
			temp$BrainmeQTL<-temp$BP%in%thischr$SNP_BP
		}
	}
	collectplus<-rbind(collectplus,temp)
}
save(collectplus,file="collectplus.rda")
save(asdgenelist,file="asdgenelist.rda")			
		
		
#####################################################################
#####################################################################
###########################PART 2####################################
#####################################################################
#####################################################################

load("collectplus.rda")

#Define parameters. 
MAFthresh<-0.0275 #Determined previously; MAF threshold for SNPs. 
genotypedirectory<-"" #Where the genotype files are. 
pgcresults<-"" #A FOLDER TO RETAIN PGC OBJECTS THAT HAVE A DIRECT OR PROXY OVERLAP
pruningfolder<-"" #A FOLDER TO CONTAIN THE SNP LISTS THAT SHOULD BE RETAINED IN THE STUDY SNPs 

#Load objects generated in power calculations. 
load("totalmaflist.rda")
load("totalsnppos.rda")
load("totalsnpchr.rda")

#Location of SNAP output files. 
myfiles<-list.files(path="")

bloodsnps<-list() #Initialize list of SNPs to be looked at for their meQTL status in peripheral blood. 
for(chr in 1:22){
	truncres<-collectplus[which(collectplus$CHR==chr),]
	if(nrow(truncres)>0){
		#Load genotype data. 
		library(sqldf)
		f<-file(paste(genotypedirectory,"SEED_MethGeno_",chr,".tped",sep="")) 
		SNPgeno <- sqldf("select * from f", dbname = tempfile(), file.format = list(header = F, row.names = F, sep=" "))
		closeAllConnections() 
		SNPgeno<-SNPgeno[,c(2,4)]
		thissnppos<-totalsnppos[which(totalsnpchr==paste("chr",chr,sep=""))]
		thismaflist<-totalmaflist[which(totalsnpchr==paste("chr",chr,sep=""))]
		SNPgeno<-SNPgeno[which(thismaflist>MAFthresh),] #Limit SNPs to those past MAF threshold. 
		
		#Initialize separate variables for direct and proxy SNPs. mSNP FOR SNP DATA FROM METHYATLION DATA
		thischr<-grep(paste("SNAP_",chr,"_",sep=""),myfiles,fixed=TRUE)
		mSNPproxy<-matrix(0,nrow(SNPgeno),(length(thischr)))
		PGCproxy<-matrix(0,nrow(truncres),(length(thischr)))
		
		#COMPUTE DIRECT OVERLAP EASILY
		mSNPdirect<-SNPgeno[,1]%in%truncres$SNP
		PGCdirect<-truncres$SNP%in%SNPgeno[,1]
		
		#Unlike in 'PGCOverlap.r', we aren't just interested in the best proxy that we have measured/imputed, but rather all of them. 
		collectproxies<-list()
		for (x in 1:length(thischr)){ #LOOP THROUGH SNAP RESULTS FOR THIS CHR
			message("Processing SNAP result ",x," for chr ",chr)
			temp<-read.table(paste("filepath",myfiles[thischr[x]],sep=""),header=TRUE,stringsAsFactors=FALSE) #LOAD SNAP RESULT
			temp<-temp[-which(temp$Distance==0),] #GET RID OF ROWS WHERE SNP IS PROXY OF ITSELF
			temp<-temp[which(temp$Proxy%in%SNPgeno[,1]),] #LIMIT TO SNPs THAT ARE ACTUALLY IN OUT STUDY POPULATION
			if(sum(temp$SNP%in%truncres$SNP)>0){
				temp<-temp[which(temp$SNP%in%truncres$SNP),]
				collectproxies[[x]]<-unique(temp$Proxy)
			}
			PGCproxy[,x]<-truncres$SNP%in%temp$SNP
		}
		collectproxies<-unlist(collectproxies) #These are the SNPs that are proxies of the PGC SNPs that we have in our data
		PGCproxy<-rowSums(PGCproxy)
		collectproxies<-unique(c(truncres$SNP[which(PGCdirect==1)],collectproxies)) #Final list to send out
		collectproxies<-data.frame(rep(chr,length(collectproxies)),collectproxies,SNPgeno[match(collectproxies,SNPgeno[,1]),2])
		colnames(collectproxies)<-c("CHR","SNP","POS")
		bloodsnps[[chr]]<-collectproxies
	}
}
bloodsnps<-do.call("rbind",bloodsnps)
save(bloodsnps,file="bloodsnps.rda")

#We then do the same in cord blood and fetal brain tissues to define the following objects:
	#cordsnps.rda		A vector of SNPs that are direct or proxy overlaps of PGC SNPs that are available in the cord blood data. 
	#brainsnps.rda		A vector of SNPs that are direct or proxy overlaps of PGC SNPs that are available in the fetal brain data.  	


#####################################################################
#####################################################################
###########################PART 3####################################
#####################################################################
#####################################################################

#Load 'totalbrain' object - all SNPs that passed QC in fetal brain data
#(obtained from study authors) which is needed in this instance just to map to meQTL results. 
library(GenomicRanges)
load("totalbrain.rda")
#Load SNPs that are direct or proxy overlaps of PGC SNPs that are available in each tissue. 
load("bloodsnps.rda")
load("cordsnps.rda") 
load("brainsnps.rda")

#Go through SNP lists and collect which CpG sites are SNP controlled in each tissue. 
meQTLfilelist<-c("/SEED_meQTLs/","/EARLI_meQTLs/","fetalbrain") #Location of meQTL results
snplists<-list(bloodsnps,cordsnps,brainsnps)
templabels<-c("BloodmeQTL","CBmeQTL","BrainmeQTL")
pvaluecutoffs<-c(1E-5,2.7E-6,1E-8) #meQTL p-value cutoffs
mappingcgs.blood<-list()
mappingcgs.cord<-list()
mappingcgs.brain<-list()
totalmappings<-list(mappingcgs.blood,mappingcgs.cord,mappingcgs.brain)
for (chr in 1:22){
	print(chr)
	for (y in 1:2){
	mappingobj<-totalmappings[[y]]
	temp<-snplists[[y]]
	temp<-temp[which(temp$CHR==chr),]
	temp$POS<-as.numeric(as.character(temp$POS))
	temp$SNP<-as.character(temp$SNP)
	if(nrow(temp)>0){
		if(y<3){
			thischr<-grep(paste("restotal",chr,"_",sep=""),list.files(meQTLfilelist[y]),fixed=TRUE)
			myres<-lapply(1:length(thischr),function(x){load(paste(meQTLfilelist[y],list.files(meQTLfilelist[y])[thischr[x]],sep=""))
				ret<-ret[which(ret$pvalue<pvaluecutoffs[y]),]
				iteration<-c()
				for (i in 1:nrow(temp)){
					matches<-which(ret$snps==temp$POS[i])
					#print(length(matches))
					if (length(matches)>0){
					iteration<-c(iteration,unique(as.character(ret$gene[which(ret$snps==temp$POS[i])])))} 
				}
				return(unique(iteration))
			}
			)
			mappingobj[[chr]]<-unlist(myres)
		} 
		if(y==3){
			thisbrain<-fetalbrain[which(fetalbrain$SNP_Chr==chr),]
			thisbrain<-thisbrain[which(thisbrain$p.value<pvaluecutoffs[y]),]
			thisbrain$SNP<-totalbrain[match(thisbrain$SNP_BP,totalbrain[,3]),1]
			iteration<-c()
			for (i in 1:nrow(temp)){
				matches<-which(thisbrain$SNP==temp$SNP[i])
				if (length(matches)>0){
				iteration<-c(iteration,unique(as.character(thisbrain$ProbeID[which(thisbrain$SNP==temp$SNP[i])])))} 
			}		
			mappingobj[[chr]]<-unique(iteration)
		}
	
	}
	totalmappings[[y]]<-mappingobj
	}
}
totalmappings[[1]]<-unlist(totalmappings[[1]])
totalmappings[[2]]<-unlist(totalmappings[[2]])
totalmappings[[3]]<-unique(fetalbrain$ProbeID[which(fetalbrain$SNP_BP%in%brainsnps)])

taggedcgs<-totalmappings
save(taggedcgs,file="taggedcgs.rda")
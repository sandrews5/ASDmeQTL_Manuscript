#Here we will annotate ASD-related PGC SNPs (p < 1 E-4) with information
#regarding their meQTL status, genes implicated by SNPs, and genes implicated by
#meQTL targets. 

#######################################################################
#######################Supplementary Table 5###########################
#######################################################################
#Here we want to report the degree of evidence for meQTLs in ASD-related (p < 1e-4)
#SNPs. We create mutually exclusive categories binning meQTL evidence. For example,
#how many SNPs are meQTLs in all three tissues? How many are meQTLs in just peripheral blood?
#We stress here that these results simply constitute the extent of evidence so far, and not
# a concrete statement as to which SNPs are definitively meQTLs. This intepretation is based
# in large part on the limitations and/or characterisitics of cross-tissue meQTL studies; please see
#the Discussion for more thoughts on this. 
load("collectplus.rda") #Produced in the script /GOBasedEnrichment/Define_ASDmeQTLtargets.r

#First, collect information on a SNP level. 
templabels<-c("BloodmeQTL","CBmeQTL","BrainmeQTL")
meQTLclass<-collectplus[,templabels]
meQTLclass$BloodmeQTL<-ifelse(meQTLclass$BloodmeQTL>0,1,0)
meQTLclass$CBmeQTL<-ifelse(meQTLclass$CBmeQTL>0,1,0)
table(meQTLclass)

#Next, collect information on an independent site/loci level. Recall that we binned
#SNPs into independent sites based on 1000 Genomes recombination hot spot data. 
locibased<-(matrix(0,0,3))
for (chr in 1:22){
	chrloci<-collectplus[which(collectplus$CHR==chr),]
	compressed<-lapply(unique(chrloci$Assignment),function(mull){
		region<-chrloci[which(chrloci$Assignment==mull),templabels]
		region<-colSums(region)
		return(region)
	})
	compressed<-do.call("rbind",compressed)
	locibased<-rbind(locibased,compressed)
}
locibased<-ifelse(locibased>0,1,0)
locibased<-data.frame(locibased)
colnames(locibased)<-templabels
table(locibased)

#######################################################################
#######################Supplementary Table 6###########################
#######################################################################
#Next we want to generate meQTL evidence for every ASD-associated (PGC p-value < 1E-4)
#locus. This table contains the following locus:
	#Locus: Genomic coordinates for locus
	#MinASD_pvalue: The most significant SNP to ASD p-value as identified by PGC in this locus
	#SNP_Genes: The genes implicated by the ASD-associated SNPs
	#isBloodmeQTL: Indicator variable noting if an ASD-associated SNP in this locus is also a meQTL in peripheral blood at FDR 5% significance
	#isCBmeQTL: Indicator variable noting if an ASD-associated SNP in this locus is also a meQTL in cord blood at FDR 5% significance
	#isBrainmeQTL: Indicator variable noting if an ASD-associated SNP in this locus is also a meQTL in fetal brain at meQTL p-value of 1E-8
	#minBloodpval: The minimum meQTL p-value for ASD-associated SNP to CpG relationships for this locus in peripheral blood
	#CpG_Genes_Blood: The genes implicated by the CpG sites controlled by the ASD-associated SNPs at FDR = 5% in peripheral blood
	#minCordpval: The minimum meQTL p-value for ASD-associated SNP to CpG relationships for this locus in cord blood
	#CpG_Genes_Cord: The genes implicated by the CpG sites controlled by the ASD-associated SNPs at FDR = 5% in cord blood
	#minBrainpval: The minimum meQTL p-value for ASD-associated SNP to CpG relationships for this locus in fetal brain
	#CpG_Genes_Brain: The genes implicated by the CpG sites controlled by the ASD-associated SNPs at meQTL p-value of 1E-8 in fetal brain

library(minfi)
load("refFlat.clean.new.order.rda") #Produced in the script /GOBasedEnrichment/AddtlObjects.r
load("collectplus.rda") #Produced in the script /GOBasedEnrichment/Define_ASDmeQTLtargets.r
load("asdgenelist.rda") #Produced in the script /GOBasedEnrichment/Define_ASDmeQTLtargets.r
#Load 450k GenomicMethylSet so that we can positions of measured CpG sites. 
load("object.qnorm_preQC.Batch1.rda")
pos<-as.numeric(start(object))
#Load fetal brain meQLT data. 
windowsize<-1000000
fetalbrain<-read.csv("All_Imputed_BonfSignificant_mQTLs.csv",header=TRUE)
fetalbrain<-fetalbrain[which(fetalbrain$SNP_Chr==fetalbrain$DNAm_CHR),]
fetalbrain$Space<-abs(fetalbrain$SNP_BP-fetalbrain$DNAm_BP)
fetalbrain<-fetalbrain[which(fetalbrain$Space<=windowsize),]
#Vector of filepaths delimiting location of meQTL data
meQTLfilelist<-c("/SEED_meQTLs/","/meQTL_Paper/EARLI_meQTLs/","fetalbrain")
templabels<-c("BloodmeQTL","CBmeQTL","BrainmeQTL")
pvaluecutoffs<-c(1E-5,2.7E-6,1E-8)
#Initialize lists to capture meQTL target information, meQTL p-values, 
#and genes that meQTL targets harbor to across tissue type. 
mappingcgs.blood<-list()
mappingcgs.cord<-list()
mappingcgs.brain<-list()
totalmappings<-list(mappingcgs.blood,mappingcgs.cord,mappingcgs.brain)
pvaluecgs.blood<-list()
pvaluecgs.cord<-list()
pvaluecgs.brain<-list()
totalpvalues<-list(pvaluecgs.blood,pvaluecgs.cord,pvaluecgs.brain)
bloodgenes<-list()
cordgenes<-list()
braingenes<-list()
totalgenes<-list(bloodgenes,cordgenes,braingenes)
library(GenomicRanges)
for (chr in 1:22){
	#Grab results from index chromosomes and format as genomic ranges object. 
	temp<-collectplus[which(collectplus$CHR==chr),]
	genelist<-refFlat.clean.new.order[which(refFlat.clean.new.order$CHR==chr),]
	gr.genes<-GRanges(seqnames = Rle(rep(chr,nrow(genelist))), ranges = IRanges(start=as.numeric(genelist$txStart), end=as.numeric(genelist$txEnd)))		
	for (y in 1:3){
	#Working with collecting objects for this tissue. 
	mappingobj<-totalmappings[[y]]
	pvalueobj<-totalpvalues[[y]]
	geneobj<-totalgenes[[y]]
		if(y<3){
			#For the peripheral blood and cord blood meQTL results, go through and look for the ASD SNPs, and 
			#grab the position of the CpG sites that they are associated with (if any). 
			thischr<-grep(paste("restotal",chr,"_",sep=""),list.files(meQTLfilelist[y]),fixed=TRUE)
			myres<-lapply(1:length(thischr),function(x){load(paste(meQTLfilelist[y],list.files(meQTLfilelist[y])[thischr[x]],sep=""))
				ret<-ret[which(ret$pvalue<pvaluecutoffs[y]),]
				cgpos<-pos[match(ret$gene,rownames(object))]
				iteration<-list()
				for (i in 1:nrow(temp)){
					matches<-which(ret$snps==temp$BP[i])
					if (length(matches)>0){
					iteration[[i]]<-cgpos[which(ret$snps==temp$BP[i])]} else{iteration[[i]]<-NA}
				}
				return(iteration)}
			)
			compile<-list()
			for (j in 1:nrow(temp)){
				compile[[j]]<-Filter(Negate(is.na),unlist(lapply(myres,function(x){return(x[[j]])})))		
			}
			mymax<-max(unlist(lapply(compile,length)))
			finished<-do.call("rbind",lapply(compile,function(x){length(x)<-mymax; return(x)}))
			mappingobj[[chr]]<-finished	
			
			#Then do the same thing but this time grab the most significant p-value for the meQTL associations
			#(if any). 
			thischr<-grep(paste("restotal",chr,"_",sep=""),list.files(meQTLfilelist[y]),fixed=TRUE)
			myres<-lapply(1:length(thischr),function(x){load(paste(meQTLfilelist[y],list.files(meQTLfilelist[y])[thischr[x]],sep=""))
				ret<-ret[which(ret$pvalue<pvaluecutoffs[y]),]
				cgpos<-pos[match(ret$gene,rownames(object))]
				iteration<-list()
				for (i in 1:nrow(temp)){
					matches<-which(ret$snps==temp$BP[i])
					if (length(matches)>0){
					iteration[[i]]<-ret$pvalue[which(ret$snps==temp$BP[i])]} else{iteration[[i]]<-NA}
				}
				return(iteration)}
			)
			compile<-list()
			for (j in 1:nrow(temp)){
				compile[[j]]<-Filter(Negate(is.na),unlist(lapply(myres,function(x){return(x[[j]])})))		
			}
			mymax<-max(unlist(lapply(compile,length)))
			pfinished<-do.call("rbind",lapply(compile,function(x){length(x)<-mymax; return(x)}))
			pvalueobj[[chr]]<-pfinished	
		} else{
			#Here we are doing the exame same thing (grab meQTL target position information, most significant
			#p-value for the fetal brain data. We just have to write separate code to extract this information
			#because those results are formatted differently. 
			thisbrain<-fetalbrain[which(fetalbrain$SNP_Chr==chr),]
			iteration<-list()
			for (i in 1:nrow(temp)){
				matches<-which(thisbrain$SNP_BP==temp$BP[i])
				if (length(matches)>0){
				iteration[[i]]<-thisbrain$DNAm_BP[which(thisbrain$SNP_BP==temp$BP[i])]} else{iteration[[i]]<-NA}
			}
			mymax<-max(unlist(lapply(iteration,length)))
			finished<-do.call("rbind",lapply(iteration,function(x){length(x)<-mymax; return(x)}))
			mappingobj[[chr]]<-finished	
			
			thisbrain<-fetalbrain[which(fetalbrain$SNP_Chr==chr),]
			iteration<-list()
			for (i in 1:nrow(temp)){
				matches<-which(thisbrain$SNP_BP==temp$BP[i])
				if (length(matches)>0){
				iteration[[i]]<-thisbrain$p.value[which(thisbrain$SNP_BP==temp$BP[i])]} else{iteration[[i]]<-NA}
			}
			mymax<-max(unlist(lapply(iteration,length)))
			pfinished<-do.call("rbind",lapply(iteration,function(x){length(x)<-mymax; return(x)}))
			pvalueobj[[chr]]<-pfinished	
		}
	#Then we want to take the meQTL target position information and conduct a findOverlaps() query with the 
	#RefSeq gene object we previoulsy created and loaded. Finally we retain that information. 
	totalmappings[[y]]<-mappingobj
	totalpvalues[[y]]<-pvalueobj
	workwithme<-c(t(finished)) #So it attaches the 2nd row after the 1st row, and the 3rd row after that, etc.
	tempindicator<-rep(1:nrow(temp),each=ncol(finished))
	tempindicator<-tempindicator[-which(is.na(workwithme))]
	workwithme<-workwithme[-which(is.na(workwithme))]
		#workwithme: cg positions that are controlled by meQTLs
		#tempindicator: what row in temp do these belong to
	gr.cpgs<-GRanges(seqnames = Rle(rep(chr,length(workwithme))), ranges = IRanges(start=as.numeric(workwithme), end=as.numeric(workwithme)))
	ovl<-findOverlaps(gr.genes,gr.cpgs)
	assigngenes<-lapply(unique(temp$Assignment),function(x){
		peaksnps<-which(temp$Assignment==x)
		downtocgs<-which(tempindicator%in%peaksnps)
		relevanthits<-which(subjectHits(ovl)%in%downtocgs)
		if(length(relevanthits)>0){
		return(unique(genelist$Gene[queryHits(ovl)[relevanthits]]))
		} else{return(NA)}				
	})
	geneobj[[chr]]<-assigngenes
	totalgenes[[y]]<-geneobj
	}
}
save(totalmappings,file="/GLE/totalmappings.rda")
save(totalpvalues,file="/GLE/totalpvalues.rda")
save(totalgenes,file="/GLE/totalgenes.rda")

#Now format the table. 
load("collectplus.rda")
load("asdgenelist.rda")
load("/GLE/totalmappings.rda")
load("/GLE/totalpvalues.rda")
load("/GLE/totalgenes.rda")

#This table will have one row per locus/independent site. So we take advantage of the 
#binning of SNPs we performed earlier. 
locustable<-data.frame(matrix(0,0,15))
for (chr in 1:22){
	print(chr)
	temp<-collectplus[which(collectplus$CHR==chr),]
	collapseloci<-lapply(unique(temp$Assignment),function(x){ #here is where we invoke the loci that we binned the SNPs into
		template<-data.frame(matrix(0,1,15))
		colnames(template)<-c("Locus","MinASD_pvalue","SNP_Genes",
			"isBloodmeQTL","isCBmeQTL","isBrainmeQTL",
			"CpG_Genes_Blood","CpG_Genes_Cord","CpG_Genes_Brain",
			"minBloodpval","minCordpval","minBrainpval",
			"MinCGpos","MaxCGpos","FullLocus")
		thislocus<-temp[which(temp$Assignment==x),]
		template[1]<-paste("chr",chr,":",min(thislocus$BP),"-",max(thislocus$BP),sep="")
		template[2]<-min(thislocus$P)
		template[3]<-paste(paste0(unlist(asdgenelist[[chr]][[x]]),","),collapse=" ")
		template[4]<-ifelse(sum(thislocus$BloodmeQTL)>0,1,0)
		template[5]<-ifelse(sum(thislocus$CBmeQTL)>0,1,0)
		template[6]<-ifelse(sum(thislocus$BrainmeQTL)>0,1,0)
		genecollect<-lapply(1:3,function(y){
			return(paste(paste0(unlist(totalgenes[[y]][[chr]][[x]]),","),collapse=" "))
		})
		template[7:9]<-unlist(genecollect)
		template[10:12]<-NA
		if(sum(template[4:6]>0)){ #If there are meQTLs to bin, loop through the 3 different tissues and grab pvalue info.  
			poscollect<-unlist(lapply(1:3,function(y){
				out<-NA
				markers<-which(temp$Assignment==x)
				if(sum(is.na(totalpvalues[[y]][[chr]][markers,]))<length(totalpvalues[[y]][[chr]][markers,])){
				out<-min(totalpvalues[[y]][[chr]][markers,],na.rm=TRUE)}
				return(out)
			}))
			template[10:12]<-poscollect
		}
		template[13:15]<-NA
		if(sum(template[4:6]>0)){ #Again, if there are meQTLs to bin, loop through the 3 different tissues and grab meQTL target position info.  
			poscollect<-lapply(1:3,function(y){
				markers<-which(temp$Assignment==x)
				out<-c(NA,NA)
				if(sum(is.na(totalmappings[[y]][[chr]][markers,]))<length(totalmappings[[y]][[chr]][markers,])){
				out<-range(totalmappings[[y]][[chr]][markers,],na.rm=TRUE)}
				return(out)
			})
			poscollect<-do.call("rbind",poscollect)
			template[13]<-min(poscollect[,1],na.rm=TRUE)
			template[14]<-max(poscollect[,2],na.rm=TRUE)
			template[15]<-paste("chr",chr,":",min(c(thislocus$BP,as.numeric(template[13]))),"-",max(c(thislocus$BP,as.numeric(template[14]))),sep="")
		}
		return(template)
	}
	)
	collapseloci<-do.call("rbind",collapseloci)
	locustable<-rbind(locustable,collapseloci)
}

#Finally we want to to sort the table so that regions with the most meQTL information are at the top. 
#First we have regions with meQTLs across all 3 tissues, and then regions with 2, and so on. 
#For those with 2 we prioritize the instances where one of the 2 is a fetal brain meQTL. 
locustable$Count<-rowSums(locustable[,c("isBloodmeQTL","isCBmeQTL","isBrainmeQTL")])
locustable<-locustable[order(locustable$Count,decreasing=TRUE),]
locustable$Count<-ifelse((locustable$Count==2 & locustable$isBrainmeQTL ==1), locustable$Count+1,locustable$Count)
locustable<-locustable[order(locustable$Count,decreasing=TRUE),]
write.csv(locustable,file="/GLE/locustable.csv")

save(locustable,file="/GLE/locustable.rda")

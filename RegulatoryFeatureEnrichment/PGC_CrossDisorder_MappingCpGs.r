#In this script we will identify the CpG sites that are associated with
#the SNPs in each study that were determined to be direct or proxy overlaps
#of the suggestively significant SNPs from the PGC cross disorder results. 

load("/CrossDisorder/bloodsnps.rda") 
load("/CrossDisorder/cordsnps.rda")
load("/CrossDisorder/brainsnps.rda")
load("/CrossDisorder/lungsnps.rda")

#Load the fetal brain meQTL results and limit to similar window size as before. 
windowsize<-1000000
fetalbrain<-read.csv("All_Imputed_BonfSignificant_mQTLs.csv",header=TRUE)
fetalbrain<-fetalbrain[which(fetalbrain$SNP_Chr==fetalbrain$DNAm_CHR),]
fetalbrain$Space<-abs(fetalbrain$SNP_BP-fetalbrain$DNAm_BP)
fetalbrain<-fetalbrain[which(fetalbrain$Space<=windowsize),]

#Load 'totalbrain' object - all SNPs that passed QC in fetal brain data
#(obtained from study authors) which is needed in this instance just to map to meQTL results. 
load("totalbrain.rda")

#Load lung meQTL results. 
lungresults<-read.csv("Lung_meQTLs.csv",header=TRUE)

#Go through SNP lists and collect which CpG sites are SNP controlled in each tissue. 
meQTLfilelist<-c("/SEED_meQTLs/","/EARLI_meQTLs/","fetalbrain") #Location of meQTL results
snplists<-list(bloodsnps,cordsnps,brainsnps)
templabels<-c("BloodmeQTL","CBmeQTL","BrainmeQTL")
pvaluecutoffs<-c(1E-5,2.7E-6,1E-8,4E-5) #FDR = 5% p-value cutoffs or printing threshold for fetal brain results. 
mappingcgs.blood<-list()
mappingcgs.cord<-list()
mappingcgs.brain<-list()
mappingcgs.lung<-list()
totalmappings<-list(mappingcgs.blood,mappingcgs.cord,mappingcgs.brain,mappingcgs.lung)
for (chr in 1:22){
	print(chr)
	for (y in 1:4){
	mappingobj<-totalmappings[[y]]
	temp<-snplists[[y]]
	temp<-temp[which(temp$CHR==chr),]
	temp$POS<-as.numeric(as.character(temp$POS))
	temp$SNP<-as.character(temp$SNP)
	if(nrow(temp)>0){
		if(y<3){ #Gather meQTL target information for peripheral blood and cord blood. 
			thischr<-grep(paste("restotal",chr,"_",sep=""),list.files(meQTLfilelist[y]),fixed=TRUE)
			myres<-lapply(1:length(thischr),function(x){load(paste(meQTLfilelist[y],list.files(meQTLfilelist[y])[thischr[x]],sep="")) #CHANGE HERE
				ret<-ret[which(ret$pvalue<pvaluecutoffs[y]),]
				iteration<-c()
				for (i in 1:nrow(temp)){
					matches<-which(ret$snps==temp$POS[i])
					if (length(matches)>0){
					iteration<-c(iteration,unique(as.character(ret$gene[which(ret$snps==temp$POS[i])])))} 
				}
				return(unique(iteration))
			}
			)
			mappingobj[[chr]]<-unlist(myres)
		} 
		if(y==3){ #Gather meQTL target information for fetal brain. 
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
		if(y==4){ #Gather meQTL target information for lung. 
			thislung<-lung[which(lung$chr_SNP==chr),]
			thislung<-thislung[which(thislung$p_EAGLE<pvaluecutoffs[y]),]
			iteration<-c()
			for (i in 1:nrow(temp)){
				matches<-which(thislung$SNP==temp$SNP[i])
				if (length(matches)>0){
				iteration<-c(iteration,unique(as.character(thislung$CpG_probe[which(thislung$SNP==temp$SNP[i])])))} 
			}		
			mappingobj[[chr]]<-unique(iteration)
		}
		totalmappings[[y]]<-mappingobj
	}
	}
}
totalmappings[[1]]<-unlist(totalmappings[[1]]) #peripheral blood meQTL targets of PGC cross disorder SNPs
totalmappings[[2]]<-unlist(totalmappings[[2]]) #cord blood meQTL targets of PGC cross disorder SNPs
totalmappings[[3]]<-unlist(totalmappings[[3]]) #fetal brain meQTL targets of PGC cross disorder SNPs
totalmappings[[4]]<-unlist(totalmappings[[4]]) #lung meQTL targets of PGC cross disorder SNPs

taggedcgs<-totalmappings
save(taggedcgs,file="/CrossDisorder/taggedcgs_CrossDisorder.rda") #Save to invoke later on. 

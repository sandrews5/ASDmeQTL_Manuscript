#In this script, we will peform an enrichment analysis looking at degree of 
#overlap with various regulatory features. In this script, the analysis is:
	#meQTL targets vs. non-meQTL targets

#Load the vectors of tested CpG sites generated in '/GOBasedEnrichment/Define_TestedCpGs.r'
load("testedblood.rda")
load("testedcord.rda")
load("testedbrain.rda")
load("testedlung.rda")

#Load the fetal brain meQTL results and limit to similar window size as before. 
windowsize<-1000000
fetalbrain<-read.csv("measuredmeQTLS.csv",header=TRUE)
fetalbrain<-fetalbrain[which(fetalbrain$SNP_Chr==fetalbrain$DNAm_CHR),]
fetalbrain$Space<-abs(fetalbrain$SNP_BP-fetalbrain$DNAm_BP)
fetalbrain<-fetalbrain[which(fetalbrain$Space<=windowsize),]

#Load lung meQTL results. 
lungresults<-read.csv("Lung_meQTLs.csv",header=TRUE)

#first points out the locations of different meQTL results, including the 3 downsampled results.  
first<-c("/SEED_meQTLs/","/EARLI_meQTLs/","fetalbrain","lung","/SEED_DownSample1/","/SEED_DownSample2/","/SEED_DownSample3/")
#A list of tested CpGs; in the downsampled cases it must be the union (i.e. passed QC in both) of 
#peripheral blood and the other tissue. 
pickcontrols<-list(testedblood,testedcord,testedbrain,testedlung,
		testedblood[which(testedblood%in%testedcord)],
		testedblood[which(testedblood%in%testedbrain)],
		testedblood[which(testedblood%in%testedlung)])

#meQTL p-values corresponding to the dataset in question. They are (numbered by index in 'pvaluecutoffs' vector):
	#1: Peripheral blood, FDR 5% (calculated)
	#2: Cord blood, FDR 5% (calculated)
	#3: Fetal brain, as thresholded in results
	#4: Lung, FDR 5% (as reported in paper)
	#5: Peripheral blood downsampled to cord blood, FDR 5% (calculated)
	#6: Peripheral blood downsampled to fetal brain, as thresholded in fetal brain results
	#7: Peripheral blood downsampled to lung, FDR 5% (calculated)
 pvaluecutoffs<-c(1E-5,2.7E-6,1E-8,4E-5,2.3E-5,1E-8,2.4E-5)

#Initialize and store CpG names for those CpGs associated with SNPs at appropriate p-values cutoffs. 
#Note for y == 3 and y == 4 no explicit p-value thresholding is performed because these results are already thresholded
#at those values. 
mappingcgs.blood<-list()
mappingcgs.cord<-list()
mappingcgs.brain<-list()
mappingcgs.lung<-list()
mappingcgs.blood1<-list()
mappingcgs.blood2<-list()
mappingcgs.blood3<-list()
totalmappings<-list(mappingcgs.blood,mappingcgs.cord,mappingcgs.brain,
		mappingcgs.lung,mappingcgs.blood1,mappingcgs.blood2,mappingcgs.blood3)
for (chr in 1:22){
	for (y in 1:length(first)){
	mappingobj<-totalmappings[[y]]
		if(y%in%c(1:2,5:7)){
			thischr<-grep(paste("restotal",chr,"_",sep=""),list.files(first[y]),fixed=TRUE)
			myres<-unlist(lapply(1:length(thischr),function(x){load(paste(first[y],list.files(first[y])[thischr[x]],sep="")) #CHANGE HERE
				ret<-ret[which(ret$pvalue<pvaluecutoffs[y]),] #BONFERONNI
				return(unique(as.character(ret$gene)))}
			))
			mappingobj[[chr]]<-unique(myres)
		} 
		if(y==3){
			thisbrain<-fetalbrain[which(fetalbrain$SNP_Chr==chr),]
			thisbrain<-thisbrain[which(thisbrain$p.value<pvaluecutoffs[y]),]
			mappingobj[[chr]]<-unique(as.character(thisbrain$ProbeID))
		}
		if(y==4){
			thislung<-lungresults[which(lungresults$chr_SNP==chr),]
			thislung<-thislung[which(thislung$p_EAGLE<pvaluecutoffs[y]),]
			mappingobj[[chr]]<-unique(as.character(thislung$CpG_probe))
		}		
	totalmappings[[y]]<-mappingobj
	}
}

#Load list of probes with annotated SNPs from Chen et al. 2013. These were filtered out in the fetail brain study and so we do that in the 
#peripheral blood and cord blood data as well, because we know that these probes can affect enrichment results (see Methods and McClay 2015).
#We load them here again as a sanity check even though they were used already in producing 'testedxxx'.  
	polycg<-read.csv("Chen_PolyCpGandSBE.csv",header=TRUE,stringsAsFactors=FALSE)
	probelength<-read.csv("Chen_ProbeLength.csv",header=TRUE,stringsAsFactors=FALSE)
	probelength<-probelength[which(probelength$BASE_FROM_SBE<=10),]
remme<-unique(c(polycg$PROBE,probelength$PROBE))

#Finally here we create indicator vectors to say if the tested CpGs (as defined in 
#the list 'pickcontrols') are meQTL targets. For the overlap determination, we only designate
#a CpG site as not being a meQTL target if it was not a meQTL targe in EITHER tissue. 
#The end result will be the list 'cacoindicators.strict' which is identical in size to 'pickcontrols'
#(7 elements, corresponding elements are equal in length). 

cacoindicators.strict<-list()
grabov<-c(2,3,4)
count<-1
for (level in 1:length(totalmappings)){
	if (level < 5){
	current<-totalmappings[[level]]
	current<-unlist(current)
	current<-current[which(!current%in%remme)]
	print(table(current%in%pickcontrols[[level]]))
	thetested<-pickcontrols[[level]]%in%current
	cacoindicators.strict[[level]]<-thetested
	} else{
		#Get the overlap between the down samples blood lists and whatever we are comparing to
		withov<-grabov[count]
		print(withov)
		count<-count+1
		current<-totalmappings[[level]]
		current<-unlist(current)
		current2<-totalmappings[[withov]]
		current2<-unlist(current2)
		current<-current[which(current%in%current2)]
			current<-current[which(!current%in%remme)]
		print(table(current%in%pickcontrols[[level]]))

		thetested<-pickcontrols[[level]]%in%current
		against1<-pickcontrols[[level]]%in%unlist(totalmappings[[level]])
		against2<-pickcontrols[[level]]%in%unlist(totalmappings[[withov]])
		thenewtested<-ifelse((thetested==0 & (against1 ==1 | against2==1)),2,thetested) 
		cacoindicators.strict[[level]]<-thenewtested
	}
}
save(cacoindicators.strict,file="/CrossTissue/cacoindicators.strict.rda")

####################################################################
#Actual Enrichment testing
library(minfi)
load("RGset.rda") #Red green channel set, just so that we can load the 450k manifest from minfi
manifest<-getAnnotation(RGset)
manifest<-manifest[,c("Name","chr","pos")]
colnames(manifest)<-c("Name","CHR","Pos")

#Load the vectors of tested CpG sites generated in '/GOBasedEnrichment/Define_TestedCpGs.r'
load("testedblood.rda")
load("testedcord.rda")
load("testedbrain.rda")
load("testedlung.rda")

#first points out the locations of different meQTL results, including the 3 downsampled results.  
first<-c("/SEED_meQTLs/","/EARLI_meQTLs/","fetalbrain","lung","/SEED_DownSample1/","/SEED_DownSample2/","/SEED_DownSample3/")
#A list of tested CpGs; in the downsampled cases it must be the union (i.e. passed QC in both) of 
#peripheral blood and the other tissue. 
pickcontrols<-list(testedblood,testedcord,testedbrain,testedlung,
		testedblood[which(testedblood%in%testedcord)],
		testedblood[which(testedblood%in%testedbrain)],
		testedblood[which(testedblood%in%testedlung)])

load("/CrossTissue/cacoindicators.strict.rda")

#Testing against all DHS as in the manifest. We do this for all 7 types of meQTL targets we want to think about, 
full.stat<-rep(0,7)
full.p<-rep(0,7)
for(x in 1:7){
	full<-getAnnotation(RGset)
	complete<-pickcontrols[[x]]
	full<-full[match(complete,full$Name),c("Name","chr","pos","DHS")]
	full$status<-0
	full$status[which(full$Name%in%complete[which(cacoindicators.strict[[x]]==1)])]<-1
	
	#Intialize a 2x2 table called mytab which has the follow structure:
	# 1 0
	#1
	#0
	#Where 1/0 is meant to be an indicator variable. 
	#We will run a simple fisher's exact test and retain effect size and p-value. 
	mytab<-matrix(0,2,2)
	mytab[1,1]<-length(which(full$DHS==TRUE & full$status==1))
	mytab[1,2]<-length(which(full$DHS==TRUE & full$status==0))
	mytab[2,1]<-length(which(full$DHS!=TRUE & full$status==1))
	mytab[2,2]<-length(which(full$DHS!=TRUE & full$status==0))
	res<-fisher.test(mytab)
	full.stat[x]<-res$estimate
	full.p[x]<-res$p.value
}
#Names indicate what meQTl target lists we are testing. 
names(full.stat)<-names(full.p)<-c("PB","CB","FB","L","PB-CB","PB-FB","PB-L")

library(rtracklayer)
#Set to the directory containing the bed files for various chromatin marks and TFBSs, and DNaseI HS sites.
setwd("/FunctionalEnrichment")  
files <- list.files()
totallist<-c("CD14_monocytes_Pk_DNaseIHS","CD4_Pk_DNaseIHS",
	"Cerebellum_Pk_DNaseIHS","Cerebrum_Pk_DNaseIHS",
	"FrCortex_Pk_DNaseIHS","IMR90_Pk_DNaseIHS",
	"H3K36me3_Blood","H3K36me3_FetalBrain","H3K36me3_Lung",
	"H3K4me1_Blood","H3K4me1_FetalBrain","H3K4me1_Lung",
	"H3K4me3_Blood","H3K4me3_FetalBrain","H3K4me3_Lung",
	"H3K27me3_Blood","H3K27me3_FetalBrain",
	"H3K9me3_Blood","H3K9me3_FetalBrain","H3K9me3_Lung")

library(GenomicRanges)
FunctionalEnrichment<-function(regions,interest,controls,manifest){
	#gr.regions<-GRanges(seqnames = Rle(regions$CHR), ranges = IRanges(start=as.numeric(regions$Start), end=as.numeric(regions$Stop)))		
	gr.regions<-regions
	gr.interest<-GRanges(seqnames = Rle(manifest$CHR[match(interest,manifest$Name)]), ranges = IRanges(start=as.numeric(manifest$Pos[match(interest,manifest$Name)]), end=as.numeric(manifest$Pos[match(interest,manifest$Name)])))
	gr.controls<-GRanges(seqnames = Rle(manifest$CHR[match(controls,manifest$Name)]), ranges = IRanges(start=as.numeric(manifest$Pos[match(controls,manifest$Name)]), end=as.numeric(manifest$Pos[match(controls,manifest$Name)])))
	ovl.interest<-findOverlaps(gr.regions,gr.interest)
	ovl.controls<-findOverlaps(gr.regions,gr.controls)
	
	#Intialize a 2x2 table called mytab which has the follow structure:
	# 1 0
	#1
	#0
	#Where 1/0 is meant to be an indicator variable. 
	#We will run a simple fisher's exact test and retain effect size and p-value. 
	mytab<-matrix(0,2,2)
	mytab[1,1]<-length(unique(subjectHits(ovl.interest)))
	mytab[1,2]<-length(unique(subjectHits(ovl.controls)))
	mytab[2,1]<-length(interest)-length(unique(subjectHits(ovl.interest)))
	mytab[2,2]<-length(controls)-length(unique(subjectHits(ovl.controls)))
	res<-fisher.test(mytab)
	return(list(est=res$estimate,p=res$p.value))	
}

#First we'll collect results for the DNaseI HS and the chromatin marks. We run the previously defined
#function "FunctionalEnrichment" for each site/mark with each definition of meQTL targets. 
#Cases here are defined as meQTL targets
#Controls here are defined as non-meQTL targets. 
enrichmat<-matrix(0,length(totallist),length(pickcontrols))
pmat<-matrix(0,length(totallist),length(pickcontrols))
for(i in 1:length(totallist)){
#for(i in 1:6){
	print(i)
	a<-import.bed(con=totallist[i])
	thisresult<-unlist(lapply(1:length(pickcontrols),function(x){
		complete<-pickcontrols[[x]]
		cases<-complete[which(cacoindicators.strict[[x]]==1)]
		controls<-complete[which(cacoindicators.strict[[x]]==0)]
		return(FunctionalEnrichment(a,cases,controls,manifest))
	}))
	enrichmat[i,]<-thisresult[grep("est",names(thisresult))]
	pmat[i,]<-thisresult[grep("p",names(thisresult))]
}
rownames(enrichmat)<-rownames(pmat)<-totallist
colnames(enrichmat)<-colnames(pmat)<-c("PB","CB","FB","L","PB-CB","PB-FB","PB-L")
head(enrichmat)
head(pmat<(0.05/181))
save(enrichmat,file="enrichmat.rda")
save(pmat,file="pmat.rda")

#Next we'll collect results for the TFBSs. We run the previously defined
#function "FunctionalEnrichment" for each site/mark with each definition of meQTL targets. 
#Cases here are defined as meQTL targets
#Controls here are defined as non-meQTL targets. 
tfbs<-import.bed(con="TFBS")
totalfactors<-unique(tfbs$name)
enrichmat.tfbs<-matrix(0,length(totalfactors),length(pickcontrols))
pmat.tfbs<-matrix(0,length(totalfactors),length(pickcontrols))
for(i in 1:length(totalfactors)){
	print(i)
	a<-tfbs[which(tfbs$name==totalfactors[i]),]
	thisresult<-unlist(lapply(1:length(pickcontrols),function(x){
		complete<-pickcontrols[[x]]
		cases<-complete[which(cacoindicators.strict[[x]]==1)] #vector of cg namees
		controls<-complete[which(cacoindicators.strict[[x]]==0)]  #vector of cg namees
		return(FunctionalEnrichment(a,cases,controls,manifest))
	}))
	enrichmat.tfbs[i,]<-thisresult[grep("est",names(thisresult))]
	pmat.tfbs[i,]<-thisresult[grep("p",names(thisresult))]
}
rownames(enrichmat.tfbs)<-rownames(pmat.tfbs)<-totalfactors
colnames(enrichmat.tfbs)<-colnames(pmat.tfbs)<-c("PB","CB","FB","L","PB-CB","PB-FB","PB-L")
save(enrichmat.tfbs,file="enrichmat.tfbs.rda")
save(pmat.tfbs,file="pmat.tfbs.rda")

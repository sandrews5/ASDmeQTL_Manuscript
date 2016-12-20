#Here we will define a vector that indicates whether the CpGs in the 'testedxxx' vectors 
#(see Define_TestedCpGs.r) are SNP-controlled. 

#Load the vectors of tested CpG sites generated in 'Define_TestedCpGs.r'
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

#The GO analyses do not use meQTL target data outside of those from peripheral blood, cord blood,
#or fetal brain. Below we generate results for meQTL targets for an additional 4 sets of interest: lung, 
#peripheral blood - cord blood overlap, peripheral blood - fetal brain overlap, peripheral blood - lung overlap.
#These are used later in the regulatory element enrichment analyses. 

#When we were determining the overlap between peripheral blood meQTL targets and three different tissues,
#we did this from a downsampled subset of the peripheral blood data run at the same parameters used in the other
#tissues. For example, when comparing peripheral blood and cord blood, we ran the full peripheral blood meQTL query
#on a random sample of 121 individuals, and used the same MAF, window size, and CpG sd cutoffs used for the cord blood
#meQTL query. We used the same scripts included in this repository for these purposes. 

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
	#3: Fetal brain, bonferroni (as reported in paper)
	#4: Lung, FDR 5% (as reported in paper)
	#5: Peripheral blood downsampled to cord blood, FDR 5% (calculated)
	#6: Peripheral blood downsampled to fetal brain, bonferroni (as used in fetal brain paper)
	#7: Peripheral blood downsampled to lung, FDR 5% (calculated)
pvaluecutoffs<-c(1E-5,2.7E-6,3.69E-13,4E-5,2.3E-5,3.69E-13,2.4E-5)

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
			myres<-unlist(lapply(1:length(thischr),function(x){load(paste(first[y],list.files(first[y])[thischr[x]],sep="")) 
				ret<-ret[which(ret$pvalue<pvaluecutoffs[y]),]
				return(unique(as.character(ret$gene)))}
			))
			mappingobj[[chr]]<-unique(myres)
		} 
		if(y==3){
			thisbrain<-fetalbrain[which(fetalbrain$SNP_Chr==chr),]
			mappingobj[[chr]]<-unique(as.character(thisbrain$ProbeID))
		}
		if(y==4){
			thislung<-lungresults[which(lungresults$chr_SNP==chr),]
			mappingobj[[chr]]<-unique(as.character(thislung$CpG_probe))
		}		
	totalmappings[[y]]<-mappingobj
	}
}

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
		thetested<-pickcontrols[[level]]%in%current
		against1<-pickcontrols[[level]]%in%unlist(totalmappings[[level]])
		against2<-pickcontrols[[level]]%in%unlist(totalmappings[[withov]])
		thenewtested<-ifelse((thetested==0 & (against1 ==1 | against2==1)),2,thetested) 
		cacoindicators.strict[[level]]<-thenewtested
	}
}
save(cacoindicators.strict,file="cacoindicators.strict.rda")

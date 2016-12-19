#Here we will use the objects we created in the script AnnotatePGC.r to run the ASD enrichment testing.
#We test for the enrichment of meQTLs in ASD snps

#We will need to match on MAF and tagged CGs. So let's first
#find the cg range in the total SNP list (MAF range will definitely be 
#MAFcutoff to 0.5. 

asdenrichmentdirectory<-"" #Where to put objects relevant to performing SNP based enrichment 

#We're going to great bins that jointly acknowledge the maf bin and the 
#cg bin that each snp is in. The constraints will be the same for every chromosome. 

mafrange<-matrix(0,22,2)
cgrange<-matrix(0,22,2)
for (chr in 1:22){
	load(paste(asdenrichmentdirectory,"prunein_asd_",chr,".rda",sep=""))
	mafrange[chr,]<-range(prunein$StudyMAF)
	cgrange[chr,]<-range(prunein$CountCGs)
}

MAFcutoff<-0.0275 ###Study specific MAF cutoff. 
cgmax<-max(cgrange)


mymafcuts<-seq(MAFcutoff,0.5,0.05)
mymafcuts<-c(mymafcuts,0.5)
mycgcuts<-seq(1,cgmax,50)
mycgcuts<-c(mycgcuts,cgmax)

#Create the bins of MAF and number of 450k CpGs in proximity of SNP (count generated in previous script). 
#Also create a joint bin of the two, which is what we will match on (i.e. SNPs in null sets must match set of interest on both factors). 
for (chr in 1:22){
	load(paste(asdenrichmentdirectory,"prunein_asd_",chr,".rda",sep=""))
	prunein$MAFbin<-as.numeric(cut(prunein$StudyMAF,mymafcuts))
	prunein$CGbin<-as.numeric(cut(prunein$CountCGs,mycgcuts,include.lowest=TRUE))
	prunein$JointBin<-as.numeric(paste(prunein$MAFbin,prunein$CGbin,sep=""))
	PGCannot<-prunein[,c(10:20)]
	save(PGCannot,file=paste(asdenrichmentdirectory,"PGCannot_asd_",chr,".rda",sep=""))
}

###########################################################################
#Finally we're ready to start enrichment testing. 

asdenrichmentdirectory<-"" #Where to put objects relevant to performing SNP based enrichment 
asdcols<-c("isASD3","isASD4")
meqtlcols<-c("ismeQTL10","ismeQTL5","ismeQTL1")

#Set up doing all combinations of asdcols and meqtlcols. 
asdcols<-rep(asdcols,each=length(unique(meqtlcols)))
meqtlcols<-rep(meqtlcols,length(unique(asdcols)))

permcount<-1000 #number of permutations to perform
totalnullsets<-matrix(0,permcount,length(asdcols)) #initialize storage for all null set proportions for all examined combinations of ASD SNP p-value & meQTL p-value
obspropvector<-rep(0,length(asdcols)) #storage for proportion of meQTLs in observed (non-null) SNP sets. 
for (i in 1:length(asdcols)){
	
	#First grab all of the proportion of meQTLs in the observed data (for this particular ASD SNP defintion and meQTL definition). 
	#Also grab the bins of MAF and CpGs in proximity to match on later in the null set generation, and an indicator of what SNPs are meQTLs.
	obsmatrix<-matrix(0,22,2)
	runningbins<-c()
	totalmeqtl<-c()
	totalbins<-c()
	for (chr in 1:22){
		load(paste(asdenrichmentdirectory,"PGCannot_asd_",chr,".rda",sep=""))
		PGCannot<-PGCannot[,c(asdcols[i],meqtlcols[i],"JointBin")]
		obsmatrix[chr,1]<-sum(PGCannot[,1])
		obsmatrix[chr,2]<-length(which(PGCannot[,1]==1 & PGCannot[,2]>0))
		runningbins<-c(runningbins,PGCannot$JointBin[which(PGCannot[,1]==1)])
		totalmeqtl<-c(totalmeqtl,PGCannot[,2])
		totalbins<-c(totalbins,PGCannot$JointBin)
	}
	obsmatrix<-colSums(obsmatrix)
	obspropvector[i]<-obsmatrix[2]/obsmatrix[1]
	totalmeqtl<-ifelse(totalmeqtl>0,1,totalmeqtl)

	nullprops<-rep(0,permcount) #The proprotions of meQTL in the null sets. 
	for (k in 1:permcount){
	message("iteration ",k," of combination ",i)
		throughlist<-unique(runningbins)
		collectbins<-c()
		collectindicator<-c()
		#Cycle through (unique) MAF/CpG bin list and grab as many SNPs for the null sets as there are in the observed set for that bin. 
		#Choose randomly from all of the avaialble SNPs in that bin. Collect the indicator as to the SNPs meQTL status
		for (y in throughlist){ 
			takecount<-length(which(runningbins==y))
			temp.totalbins<-totalbins[which(totalbins==y)]
			temp.totalmeqtl<-totalmeqtl[which(totalbins==y)]
			if(takecount>length(temp.totalbins)){
				mychoices<-sample(1:length(temp.totalbins),takecount,replace=TRUE)
			} else {mychoices<-sample(1:length(temp.totalbins),takecount,replace=FALSE)}
			collectbins<-c(collectbins,temp.totalbins[mychoices])
			collectindicator<-c(collectindicator,temp.totalmeqtl[mychoices])
		}
	#Grab the proportion of meQTLs in this null set and store it. 
	nullprops[k]<-mean(collectindicator)
	}
	totalnullsets[,i]<-nullprops #Grab all the proportions from all the null sets for this ASD SNP p-value/meQTL p-value combination. 
}

obspropvector/colMeans(totalnullsets) #Fold enrichment statistic
unlist(lapply(1:6,function(desk){sum(totalnullsets[,desk]>=obspropvector[desk])}))/1000 #Pvalue (may need to flip sign depending upon direction of fold enrichment (want 'more exteme'))

save(totalnullsets,file=paste(asdenrichmentdirectory,"totalnullsets_asd.rda",sep=""))
save(obspropvector,file=paste(asdenrichmentdirectory,"obspropvector_asd.rda",sep=""))
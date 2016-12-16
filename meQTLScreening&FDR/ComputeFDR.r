resultsdirectory<-"" #Where the meQTL results are (observed data)
permresultsdirectory<-"" #Where the permuted meQTL results are 

myfiles<-list.files(path=resultsdirectory)
#Thresholds at which to compute FDR
thresholds<-c(1E-4,1E-5,1E-6,1E-7,1E-8,1E-9,1E-10,1E-11,1E-12,1E-13,1E-14,1E-15) 
fdrresults<-matrix(0,100,length(thresholds)) #100 for 100 permutations. 
N1<-rep(0,length(thresholds))
for (thresh in 1:length(thresholds)){	
	pval<-thresholds[thresh]
	truncresults2<-matrix(0,100,22)
	truncresults1<-rep(0,22)
	mysequence<-c(1:22) #sequence of chromosomes. Can change this here to only compute FDR on first 6 chromosomes, for example. 
	for (CHR in mysequence){
		#From permutation results, grab CpG site names that were significant past designated p value threshold. 
		load(paste(permresultsdirectory,"cgnames_",CHR,".rda",sep=""))
		truncresults2[,CHR]<-unlist(lapply(cgnames,function(x){relevant<-x[which(x$pvalue<pval),1]
				return(length(unique(relevant)))}))
		
		#Now collect the significant CpG sites from the observed data. 
		thischr<-grep(paste("restotal",CHR,"_",sep=""),myfiles,fixed=TRUE)
		myres<-lapply(1:length(thischr),function(x){load(paste(resultsdirectory,myfiles[thischr[x]],sep=""))
			ret<-ret[which(ret$pvalue<pval),]
			return(length(unique(ret$gene)))})
		truncresults1[CHR]<-sum(unlist(myres))		
	}
	fdrresults[,thresh]<-rowSums(truncresults2) #Get one result across all chromosomes for each pval threshold
	N1[thresh]<-sum(truncresults1) # The total number of significant CpG sites in the observed data. 
}
N2<-colMeans(fdrresults) #Mean across permutations. 
FDR<-N2/N1
cbind(thresholds,FDR)

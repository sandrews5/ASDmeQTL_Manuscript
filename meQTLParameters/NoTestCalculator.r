#A function to calculate the number of tests given input of
	#CpG probe information
	#SNP information
#at various cutoffs relating to:
	#CpG VMP sd filtering cutoffs
	#SNP MAF filtering cutoffs
	#Window sizes 
	
#INPUT PARAMETERS
	#methchr: a character vector indicating chromosome information for the 450k probes ("chr1", "chr1", etc.)
	#methpos: numeric vector indicating 450k probe position information
	#methsd: numeric vector indicating standard deviation of 450k probe across samples [NOTE: methchr, methpos, methsd, should all be same length]
	#snpchr: same as methchr but for genotype data
	#snppos: same as methpos but for genotype data
	#snpMAF: numeric vector indicating minor allele frequency of each SNP [NOTE: snpchr, snppos, snpMAF, should all be same length]
	#sdvector: numeric vector of percentiles of total standard deviation to consider
	#MAFvector: numeric vector of MAF cutoffs to consider
	#windowsizevector: numeric vector of window sizes to consider (defined as values on either side of a CpG site (i.e. input of 50kb means 50kb upstream and downstream)). Units of base pairs. 
#The function will calculate the number of tests based on every unique combination of elements from sdvector, MAFvector, windowsizevector

#Output will be a data frame showing all of the unique combinations of the elements of sdvector, MAFvector, and windowsizevector, 
		#and the subsequent number of tests and and Bonferroni corrected alpha (0.05/number of tests)

NoTestCalculator<-function(methchr,methpos,methsd,snpchr,snppos,snpMAF,sdvector=seq(0.05,0.30,by=0.05),MAFvector=c(0.01,0.05,0.10,0.20,0.30),windowsizevector=c(500000,1000000,2000000,3000000,4000000)){
	sdvector<-sort(sdvector)
	MAFvector<-sort(MAFvector)
	windowsizevector<-sort(windowsizevector,decreasing=TRUE)
	
	sdvector<-rep(sdvector,each=length(MAFvector))
	MAFvector<-rep(MAFvector,length(unique(sdvector)))
	sdvector<-rep(sdvector,each=length(windowsizevector))
	MAFvector<-rep(MAFvector,each=length(windowsizevector))
	parameters<-cbind(sdvector,MAFvector,windowsizevector)
	#print(dim(parameters))
	
	sdcutoffs<-quantile(methsd,unique(sdvector))
	
	sdchoice<-parameters[1,1]
	MAFchoice<-parameters[1,2]
	windowchoice<-parameters[1,3]
	
	tempmethchr<-methchr[which(methsd>sdcutoffs[which(unique(sdvector)==sdchoice)])]
	tempmethpos<-methpos[which(methsd>sdcutoffs[which(unique(sdvector)==sdchoice)])]
	tempsnpchr<-snpchr[which(snpMAF>MAFchoice)]
	tempsnppos<-snppos[which(snpMAF>MAFchoice)]
	
	tempmethsd<-methsd[which(methsd>sdcutoffs[which(unique(sdvector)==sdchoice)])]
	tempsnpmaf<-snpMAF[which(snpMAF>MAFchoice)]
		 
	methposlist<-lapply(1:length(unique(tempmethchr)),function(x){return(tempmethpos[which(tempmethchr==unique(tempmethchr)[x])])})
	snpposlist<-lapply(1:length(unique(tempsnpchr)),function(x){return(tempsnppos[which(tempsnpchr==unique(tempsnpchr)[x])])})
	methsdlist<-lapply(1:length(unique(tempmethchr)),function(x){return(tempmethsd[which(tempmethchr==unique(tempmethchr)[x])])})
	snpmaflist<-lapply(1:length(unique(tempsnpchr)),function(x){return(tempsnpmaf[which(tempsnpchr==unique(tempsnpchr)[x])])})

	parsedcalcs<-lapply(1:length(methposlist),function(k){
		message("Working on chromosome ",k)
		mymethchr<-methposlist[[k]]
		mymethsd<-methsdlist[[k]]
		mysnpchr<-snpposlist[[k]]
		mysnpmaf<-snpmaflist[[k]]
		totalres<-rep(0,nrow(parameters))
		breakpoint<-500
		breaks<-ceiling(length(mysnpchr)/breakpoint) #breaking up probes into breakpoint group segments
		#print(breaks)
		#print(workers)
		for (j in 1:breaks){
			message("chr ",k," on break ",j," out of ",breaks,sep="")
				if (j<breaks){
				usesnp<-mysnpchr[(((breakpoint*j)-(breakpoint-1)):(j*breakpoint))] 
				usemaf<-mysnpmaf[(((breakpoint*j)-(breakpoint-1)):(j*breakpoint))] 
				} else {
					usesnp<-mysnpchr[(((breakpoint*j)-(breakpoint-1)):length(mysnpchr))]
					usemaf<-mysnpmaf[(((breakpoint*j)-(breakpoint-1)):length(mysnpmaf))]
				}

				res<-abs(outer(usesnp,mymethchr,"-"))
				chunksize <- ceiling(nrow(parameters)/workers)
				#print(chunksize)
				fillparameters <- foreach(probes = iter(parameters, by = "row", chunksize = chunksize)) %dorng%
					{apply(probes,1,function(sophie){
						sdchoicehere<-sophie[1]
						MAFchoicehere<-sophie[2]
						windowchoicehere<-sophie[3]
						keepsnps<-which(usemaf>MAFchoicehere)
						keepprobes<-which(mymethsd>sdcutoffs[which(unique(sdvector)==sdchoicehere)])
						complete<-sum(res[keepsnps,keepprobes]<=windowchoicehere)
						return(complete)
					})}
			fillparameters<-unlist(fillparameters)
			totalres<-totalres+fillparameters
			}
		return(totalres)
	})
	testcounts<-do.call("cbind",parsedcalcs)
	testcounts<-rowSums(testcounts)
	
	parameters<-data.frame(parameters)
	colnames(parameters)<-c("MethPercSdCutoff","MAFThresh","WindowSize")
	parameters$TestCount<-testcounts
	parameters$BonfAlpha<-0.05/testcounts
	
	return(parameters)
	#return(list("testcounts" = parameters, "sdcutoffs" = sdcutoffs))
}

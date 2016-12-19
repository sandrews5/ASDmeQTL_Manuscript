#In this script, we look for proxies of all of the SNPs in the PGC results. 
#We use the SNAP software: https://www.broadinstitute.org/mpg/snap/ldsearch.php

#Load the PGC results downloaded from here: http://www.med.unc.edu/pgc/results-and-downloads

library(sqldf)
f<-file("onqah.pgcasdeuro")
pgcnew <- sqldf("select * from f", dbname = tempfile(), file.format = list(header = T, row.names = F, sep="\t"))

for (chr in 1:22){
	#Subset results to this chromosome.
	thispgc<-pgcnew[which(pgcnew$CHR==chr),]
	truncres<-thispgc[,c("SNP","BP","EUR_FRQ","OR","P")]
	colnames(truncres)<-c("SNP","BP","MAF","OR","P")
	
	#We will use the unix curl command to query the SNAP database. 
	#We only want to query 10,000 SNPs at a time so as to not over burden the database. 
	#So, we set up this breaks variable to only send 10,000 SNPs at a time. 
	#Then we write a file called 'output.txt' that we send to the SNAP server via curl. 
	#In this file we specify all of the arguments you would manually input if you 
	#were to use the SNAP website. 
	breakpoint<-10000
	breaks<-ceiling(nrow(truncres)/breakpoint) #breaking up snps into breakpoint group segments
	for (j in 1:breaks){
		if (j<breaks){
			temptrunc<-truncres[(((breakpoint*j)-(breakpoint-1)):(j*breakpoint)),"SNP"]
			} else {
				temptrunc<-truncres[(((breakpoint*j)-(breakpoint-1)):nrow(truncres)),"SNP"]
			}
		firstpart<-"searchPairwise=&snpList="
		lastpart<-"&hapMapRelease=onekgpilot&hapMapPanel=CEU&RSquaredLimit=0.8&distanceLimit=500000&downloadType=File&includeQuerySnp=on&arrayFilter=query&columnList[]=DP&columnList[]=GP&columnList[]=AM&submit=search"
		partition<-"%0D%0A"
		snplist<-paste0(as.character(temptrunc)[1:(length(temptrunc)-1)],partition)
		snplist<-paste(snplist,collapse="")
		snplist<-paste(snplist,temptrunc[length(temptrunc)],sep="")
	write.table(paste(firstpart,snplist,lastpart,sep=""),file="output.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
	system("curl -d @output.txt 'http://www.broadinstitute.org/mpg/snap/ldsearch.php' > results.txt")
	filename<-paste("/dcs01/arking/arkinglab/shared/PGC_meQTL/PGC_SNPNames/SNAP_",chr,"_",j,".txt",sep="")
	system(paste("grep -v 'WARNING' results.txt > ",filename))
	}
}
q()
n

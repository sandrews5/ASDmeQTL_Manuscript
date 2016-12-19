#Here we will make the 'snp input table' required by priority pruner (see http://prioritypruner.sourceforge.net/documentation.html)
	#First limit to the SNPs that are in the genotype files after overlap computation (i.e. either direct or proxy: 'snplist_x')
	#Then grab the PGC p-value for the SNP is overlaps or that it is a proxy for

overlapdirectory<-"" #DIRECTORY OF SNP FILES THAT JUST HAVE THOSE SNPS THAT OVERLAP PGC; set in 'PriorityPruner_MakeGenotypes.sh'
pgcresults<-"" #Location of 'truncres' objects generated in 'PGCoverlap.r'. 
pruningdirectory<-"" #Location of input tables generated here, as well as the pruned genotype files generated next in 'RunPriorityPruner.sh'

library(sqldf)
for (chr in 1:22){
	print(chr)
	#Load genotype data and initialize certain columns of input table. 
	f<-file(paste(overlapdirectory,"SEED_Overlap_",chr,".tped",sep=""))
	SNPgeno <- sqldf("select * from f", dbname = tempfile(), file.format = list(header = F, row.names = F, sep=" "))
	closeAllConnections() 
	peoplemap<-read.table(paste(overlapdirectory,"SEED_Overlap_",chr,".tfam",sep="")) 
	inputable<-data.frame(matrix(0,nrow(SNPgeno),8))
	inputable[,1]<-SNPgeno[,2] #SNP name
	inputable[,2]<-SNPgeno[,1] #chromosome
	inputable[,3]<-SNPgeno[,4] #position
	inputable[,4]<-1 #Allele 1
	inputable[,5]<-2 #Allele 2

	#Column matchme will indicate the study SNP to use for each PGC SNP (direct or proxy overlap)
	load(paste(pgcresults,"truncres_",chr,".rda",sep=""))
	truncres$matchme<-as.character(truncres$SNP)
	truncres$matchme[which(truncres$indicator==1)]<-truncres$Proxy[which(truncres$indicator==1)]

	#The 6th column is the p-value that is taken into account during the pruning step; here the PGC p-value. 
	inputable[,6]<-truncres$P[match(inputable[,1],truncres$matchme)] 
	print(sum(is.na(inputable[,6])))
	
	#The 7th column is the option to force include. We will leave this as 0 and opt out of force including anything. 
	inputable[,7]<-0
	#The 8th column is the design score, which we will set as 1 for all SNPs. 
	inputable[,8]<-1
	colnames(inputable)<-c("name","chr","pos","a1","a2","p","forceSelect","designScore")
	write.table(inputable,paste(pruningdirectory,"inputable_",chr,".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE)
}

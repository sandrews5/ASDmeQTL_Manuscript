#In this script, we will do the following:
	#Read in the PGC files, limit to those that overlap with SEED SNPs
		#directly or via proxy (using files generated by MakeSNPlists.r
	#Create the following columns in addition to the columns from the PGC lists
		#Overlap via proxy or directly? (0 - no, 1-direct, 2-proxy). Get rid of non overlap SNPs.
		#if 2 in variable above, what is the name of the proxy?
	#Print study SNP lists that will be used in pruning. 

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

#Load the PGC results downloaded from here: http://www.med.unc.edu/pgc/results-and-downloads
library(sqldf)
f<-file("onqah.pgcasdeuro")
pgcnew <- sqldf("select * from f", dbname = tempfile(), file.format = list(header = T, row.names = F, sep="\t"))

for (chr in 1:22){
	#Subset results to this chromosome.
	thispgc<-pgcnew[which(pgcnew$CHR==chr),]
	truncres<-thispgc[,c("SNP","BP","EUR_FRQ","OR","P")]
	colnames(truncres)<-c("SNP","BP","MAF","OR","P")
	
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
	
	#Write function to get one result per SNP. 
	#Pick the SNP with the highest rsquared. if tied, pick shortest distance. If still tied, randomly pick one. 
	cleanlist<-function(temp){
		snpids<-unique(temp$SNP)
		res<-lapply(snpids,function(y){
			thesesnps<-temp[which(temp$SNP==y),]
			thesesnps<-thesesnps[which(thesesnps$RSquared==max(thesesnps$RSquared)),]
			thesesnps<-thesesnps[which(thesesnps$Distance==min(thesesnps$Distance)),]
			if(nrow(thesesnps)>1){thesesnps<-thesesnps[1,]}
			return(thesesnps)
		})
		res<-do.call("rbind",res)
		return(res)
	}
	
	#Run through SNPs and find the best proxy for that SNP that is measured in the study population. 
	collectproxies<-rep(NA,nrow(truncres))
	for (x in 1:length(thischr)){ #LOOP THROUGH SNAP RESULTS FOR THIS CHR
		temp<-read.table(paste("/dcs01/arking/arkinglab/shared/PGC_meQTL/PGC_SNPNames/",myfiles[thischr[x]],sep=""),header=TRUE,stringsAsFactors=FALSE) #LOAD SNAP RESULT
		temp<-temp[-which(temp$Distance==0),] #GET RID OF ROWS WHERE SNP IS PROXY OF ITSELF
		temp<-temp[which(temp$Proxy%in%SNPgeno[,1]),] #LIMIT TO SNPs THAT ARE ACTUALLY IN OUT STUDY POPULATION
		PGCproxy[,x]<-truncres$SNP%in%temp$SNP #NOTE IN PGC SNPS WHAT WE HAVE PROXIES FOR
		matchme<-cleanlist(temp) #ONE RESULT PER SNP
		mSNPproxy[,x]<-SNPgeno[,1]%in%matchme$Proxy #KEEP TRACK OF WHAT STUDY SNPS THAT ARE PROXIES
		sendsnps<-matchme$Proxy[match(truncres$SNP,matchme$SNP)]
		collectproxies[which(!is.na(sendsnps))]<-sendsnps[which(!is.na(sendsnps))] #COLLECT NAMES OF PROXY SNPS
	}
	mSNPproxy<-rowSums(mSNPproxy)
	PGCproxy<-rowSums(PGCproxy)
	PGCindicator<-rep(0,nrow(truncres))
	PGCindicator[which(PGCdirect==1)]<-2 #DIRECT MATCH
	PGCindicator[which(PGCproxy==1&PGCdirect==0)]<-1 #PROXY MATCH
	collectproxies[which(PGCdirect==1)]<-NA #IF A PGC SNP HAS A DIRECT OVERLAP AND A PROXY, RETAIN THE DIRECT OVERLAP
	truncres$indicator<-PGCindicator
	truncres$Proxy<-collectproxies
	truncres<-truncres[-which(truncres$indicator==0),] #GET RID OF PGC SNPS THAT DON'T HAVE ANY OVERLAP
	#Save master PGC objects (limited to SNPs in study data and pointing out which SNP to use in study data - i.e. direct overlap or proxy)
	save(truncres,file=paste(pgcresults,"truncres_",chr,".rda",sep="")) 

	#Finally, write out the study SNPs lists of SNPs we are going to use (because they are direct or proxy overlap of 
	#a PGC SNP. These SNPs will be sent to pruning step. 
	truncres$CollectSNP<-as.character(truncres$SNP)
	truncres$CollectSNP[which(truncres$indicator==1)]<-truncres$Proxy[which(truncres$indicator==1)]
	write.table(unique(truncres$CollectSNP),file=paste(pruningfolder,"snplist_",chr,".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
}
q()
n

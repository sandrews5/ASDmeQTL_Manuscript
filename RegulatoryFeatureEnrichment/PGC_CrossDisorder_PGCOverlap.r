#Here we will demonstrate how we generated objects detailing which SNPs in 
#each respective study overlap with SNPs (of suggestive significance) from the 
#PGC cross disorder analysis or their proxies. 

load("/CrossDisorder/crosspgc.rda")
crosspgc<-crosspgc[which(crosspgc$P<1E-4),] #Define 1E4 p-value at suggestive significance. 

#Define parameters. 
MAFthresh<-0.0275 #Determined previously; MAF threshold for SNPs. 
genotypedirectory<-"" #Where the genotype files are. 

#Load objects generated in power calculations. 
load("totalmaflist.rda")
load("totalsnppos.rda")
load("totalsnpchr.rda")

#Location of SNAP output files. 
myfiles<-list.files(path="")

bloodsnps<-list() #Initialize list of SNPs to be looked at for their meQTL status in peripheral blood. 
for(chr in 1:22){
	truncres<-collectplus[which(collectplus$CHR==chr),]
	if(nrow(truncres)>0){
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
		
		#Unlike in 'PGCOverlap.r', we aren't just interested in the best proxy that we have measured/imputed, but rather all of them. 
		collectproxies<-list()
		for (x in 1:length(thischr)){ #LOOP THROUGH SNAP RESULTS FOR THIS CHR
			message("Processing SNAP result ",x," for chr ",chr)
			temp<-read.table(paste("filepath",myfiles[thischr[x]],sep=""),header=TRUE,stringsAsFactors=FALSE) #LOAD SNAP RESULT
			temp<-temp[-which(temp$Distance==0),] #GET RID OF ROWS WHERE SNP IS PROXY OF ITSELF
			temp<-temp[which(temp$Proxy%in%SNPgeno[,1]),] #LIMIT TO SNPs THAT ARE ACTUALLY IN OUT STUDY POPULATION
			if(sum(temp$SNP%in%truncres$SNP)>0){
				temp<-temp[which(temp$SNP%in%truncres$SNP),]
				collectproxies[[x]]<-unique(temp$Proxy)
			}
			PGCproxy[,x]<-truncres$SNP%in%temp$SNP
		}
		collectproxies<-unlist(collectproxies) #These are the SNPs that are proxies of the PGC SNPs that we have in our data
		PGCproxy<-rowSums(PGCproxy)
		collectproxies<-unique(c(truncres$SNP[which(PGCdirect==1)],collectproxies)) #Final list to send out
		collectproxies<-data.frame(rep(chr,length(collectproxies)),collectproxies,SNPgeno[match(collectproxies,SNPgeno[,1]),2])
		colnames(collectproxies)<-c("CHR","SNP","POS")
		bloodsnps[[chr]]<-collectproxies
	}
}
bloodsnps<-do.call("rbind",bloodsnps)
save(bloodsnps,file="/CrossDisorder/bloodsnps.rda")

#We then do the same in cord blood and fetal brain tissues to define the following objects:
	#cordsnps.rda		A vector of SNPs that are direct or proxy overlaps of PGC cross disorder SNPs that are available in the cord blood data. 
	#brainsnps.rda		A vector of SNPs that are direct or proxy overlaps of PGC cross disorder SNPs that are available in the fetal brain data.  	
	#lungsnps.rda		A vector of SNPs that are direct or proxy overlaps of PGC cross disorder SNPs that are available in the lung data.  	


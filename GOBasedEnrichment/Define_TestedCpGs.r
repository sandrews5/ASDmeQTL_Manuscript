#In this script we will simply define a vector of CpGs that were tested in each meQTL analysis
#for each tissue. 

##############################################
#Peripheral Blood 

#Load objects generated in power calculations. 
load("totalmaflist.rda")
load("totalsnppos.rda")
load("totalsnpchr.rda")

library(minfi)
#Load methylation object and phenotype data. 
load("object.qnorm_postQC.Batch1.rda")
object<-object.postQC
pd<-pData(object)

#Limit methylation data to white individuals (see Methods).  
seed1<-read.csv("MasterAnnotationSEEDmethylation.csv",header=TRUE)
whites<-seed1[which(seed1$Race.y=="white"),]
object.whites<-object[,which(pd$Family%in%whites$Family)] #just do this because this object was cleaned out for duplicates alread
pd.whites<-pd[which(pd$Family%in%whites$Family),]

#Calculate beta values and only include those that pass methylation sd threshold 
#defined earlier and calculated via scripts in meQTLParameters folder.
#Also remove X and Y chromosomes. 
chrnames<-as.character(seqnames(object.whites))
object.whites<-object.whites[-which(chrnames=="chrX"|chrnames=="chrY"),]
B<-getBeta(object.whites)
rowsds<-apply(B,1,sd)
chrnames<-as.character(seqnames(object.whites))
pos<-as.numeric(start(object.whites))
sdcutoffs<-quantile(rowsds,sdcutoff)
B<-B[-which(rowsds<sdcutoffs),]
chrnames<-chrnames[-which(rowsds<sdcutoffs)] 
pos<-pos[-which(rowsds<sdcutoffs)]

#Load list of probes with annotated SNPs from Chen et al. 2013 (http://www.sickkids.ca/Research/Weksberg-Lab/Publications/index.html). 
#These were filtered out in the fetail brain study and so we do that in the 
#peripheral blood and cord blood data as well, because we know that these probes can affect enrichment results (see Methods and McClay 2015). 
	polycg<-read.csv("Chen_PolyCpGandSBE.csv",header=TRUE,stringsAsFactors=FALSE)
	probelength<-read.csv("Chen_ProbeLength.csv",header=TRUE,stringsAsFactors=FALSE)
	probelength<-probelength[which(probelength$BASE_FROM_SBE<=10),]
remme<-unique(c(polycg$PROBE,probelength$PROBE))
chrnames<-chrnames[-which(rownames(B)%in%remme)]
pos<-pos[-which(rownames(B)%in%remme)]
testedblood<-rownames(B)[-which(rownames(B)%in%remme)]
save(testedblood,file="testedblood.rda")

#########################################
#We then do the same in cord blood and fetal brain tissues to define the following objects:
	#testedcord.rda		A vector of 450k CpG sites that were included in the cord blood meQTL query. 
	#testedbrain.rda		A vector of 450k CpG sites that were included in the fetal brain meQTL query (obtained from study authors). 	
	#testedlung.rda		A vector of 450k CpG sites that were included in the lung meQTL query (obtained from study authors). Note: we did not conduct the GO analysis 
			#in the lung data but we do use this object in the regulatory feature enrichment analyesis. 
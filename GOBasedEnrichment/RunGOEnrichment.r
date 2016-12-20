#In this script we will use all of the objects we have created to perform the 
#GO enrichment analysis. Then we will take format the nominally significant results
#for input into the REVIGO software. 

#Load the vectors of tested CpG sites generated in 'Define_TestedCpGs.r'
load("testedblood.rda")
load("testedcord.rda")
load("testedbrain.rda")
load("testedlung.rda")

pickcontrols<-list(testedblood,testedcord,testedbrain,testedlung,
		testedblood[which(testedblood%in%testedcord)],
		testedblood[which(testedblood%in%testedbrain)],
		testedblood[which(testedblood%in%testedlung)])

load("cacoindicators.strict.rda") #Generated in 'Define_meQTLtargets.r'
load("taggedcgs.rda") #Generated in 'Define_ASDmeQTLtargets.r'
taggedcgs[[3]]<-as.character(taggedcgs[[3]])
#Ensure that we've included only those 450k CpGs in our list of tested CpGs. 
taggedcgs[[1]]<-taggedcgs[[1]][which(taggedcgs[[1]]%in%testedblood)]
taggedcgs[[2]]<-taggedcgs[[2]][which(taggedcgs[[2]]%in%testedcord)]
taggedcgs[[3]]<-as.character(taggedcgs[[3]][which(taggedcgs[[3]]%in%testedbrain)])

#Here we are setting up the 6 GO analyses we will do, where we compare an element of
#the list 'interest' against the element of the list 'background'. These will be the inputs
#sig.cpg and all.cpg in missMethyl::gometh, respectively. The analyses are:
	#1: Peripheral blood, meQTL targets against background of all tested CpGs
	#2: Cord blood, meQTL targets against background of all tested CpGs
	#3: Fetal brain, meQTL targets against background of all tested CpGs
	#4: Peripheral blood, meQTL targets of ASD-related SNPs and their proxies against background of meQTL targets
	#5: Cord blood, meQTL targets of ASD-related SNPs and their proxies against background of meQTL targets
	#6: Fetal brain, meQTL targets of ASD-related SNPs and their proxies against background of meQTL targets
interest<-list(pickcontrols[[1]][which(cacoindicators.strict[[1]]==1)],
	pickcontrols[[2]][which(cacoindicators.strict[[2]]==1)],
	pickcontrols[[3]][which(cacoindicators.strict[[3]]==1)],
	taggedcgs[[1]],
	taggedcgs[[2]],
	taggedcgs[[3]])
background<-list(pickcontrols[[1]],
	pickcontrols[[2]],
	pickcontrols[[3]],
	pickcontrols[[1]][which(cacoindicators.strict[[1]]==1)],
	pickcontrols[[2]][which(cacoindicators.strict[[2]]==1)],
	pickcontrols[[3]][which(cacoindicators.strict[[3]]==1)])

#Perform GO enrichment, retain results in 'gores' list. 
library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(limma)
library(BiasedUrn)
gores<-list()
pdf("gomethplots.pdf",paper="a4r",width=10)
par(mai=c(1,6,1,1))
for (i in 1:6){
	print(i)
	tissue<-interest[[i]]
	ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
	sigcpgs <- unlist(tissue)
	allcpgs <- background[[i]]
	gst <- gometh(sig.cpg = sigcpgs, all.cpg = allcpgs, plot.bias=FALSE, prior.prob=TRUE)
	table(gst$FDR<0.05)
	res<-topGO(gst,number=Inf)
	gores[[i]]<-res
	res<-res[1:10,]
	res$P.DE<-ifelse(res$P.DE<1E-15,1E-15,res$P.DE)
	barplot((-1*log(res$P.DE,base=10)),horiz=TRUE,xlim=c(0,15),
	names.arg=res$Term,las=1,xlab="-log10(p-value)",cex.names=1)
}
dev.off()
save(gores,file="gores.rda")

#Save gores objects as indiviudal text file with headers that are 'GO ID' and 'p-value' and send to
#REVIGO program. 
load("gores.rda")
tissuenames<-c("PB1","CB1","FB1","PB2","CB2","FB2")
for (i in 1:6){
	tissue<-gores[[i]]
	tissue<-data.frame(rownames(tissue),tissue$P.DE)
	tissue<-tissue[which(tissue[,2]<0.05),]
	colnames(tissue)<-c("GO ID","p-value")
	write.table(tissue,file=paste("GOList_",tissuenames[i],".txt",sep=""),row.names=FALSE,quote=FALSE,col.names=TRUE)
}



#########################################################################################
#Below we read in the REVIGO output to generate several tables/figures that were in the paper. 

#These bar plots comprise the left hand sides of Supplementary Figures 2-4. 
#The right hand side are produced via R code generated from the REVGIO software itself. 
tissuenames<-c("PB1","CB1","FB1","PB2","CB2","FB2")
for (i in 1:6){
	png(filename=paste("GOMETH_",tissuenames[i],"_plot.png",sep=""),width=10, height=6,units="in",res=480,bg="white", type="cairo")
	par(mai=c(1,6,1,1))
	tissue<-read.csv(paste("GOAnnotation/REVIGO_",tissuenames[i],".csv",sep=""),header=TRUE)
	tissue<-tissue[which(tissue$eliminated==0),]
	tissue<-tissue[order(tissue$log10.p.value),]
	res<-tissue[10:1,c("description","log10.p.value")]
	colnames(res)<-c("Term","P.DE")
	res$P.DE<-ifelse(res$P.DE<(-15),-15,res$P.DE)
	barplot((-1*(res$P.DE)),horiz=TRUE,xlim=c(0,10),
	names.arg=wrap.labels(res$Term,50),las=1,xlab="-log10(p-value)",cex.names=1)
	dev.off()
}

#Next we create various tables in the paper for the comparison:
	#meQTL targets of ASD-related SNPs and their proxies against background of meQTL targets
load("gores.rda")
pb<-read.csv("REVIGO_PB2.csv",header=TRUE)
pb<-pb[which(pb$eliminated==0),]
cb<-read.csv("REVIGO_CB2.csv",header=TRUE)
cb<-cb[which(cb$eliminated==0),]
fb<-read.csv("REVIGO_FB2.csv",header=TRUE)
fb<-fb[which(fb$eliminated==0),]

siggo.pb<-gores[[4]][which(gores[[4]]$Ont=="BP"),]
siggo.cb<-gores[[5]][which(gores[[5]]$Ont=="BP"),]
siggo.fb<-gores[[6]][which(gores[[6]]$Ont=="BP"),]

siggo.pb<-siggo.pb[which(rownames(siggo.pb)%in%pb$term_ID),]
siggo.cb<-siggo.cb[which(rownames(siggo.cb)%in%cb$term_ID),]
siggo.fb<-siggo.fb[which(rownames(siggo.fb)%in%fb$term_ID),]

#Generate overlap objects
pb.cb<-as.character(pb$description[which(pb$description%in%cb$description)])
pb.fb<-as.character(pb$description[which(pb$description%in%fb$description)])
cb.fb<-as.character(cb$description[which(cb$description%in%fb$description)])
allthree<-as.character(pb.fb[which(pb.fb%in%cb.fb)])
myterms<-unique(c(pb.cb,pb.fb,cb.fb,allthree))
inPB<-as.numeric(myterms%in%pb$description)
inCB<-as.numeric(myterms%in%cb$description)
inFB<-as.numeric(myterms%in%fb$description)

#Generate rank and scaled rank for individual terms from the terms that
#overlap across 2 or 3 tissues. Do this for each tissue. 
bloodres<-matrix(NA,length(myterms),3) #RANK #Scaled rank #p-value
for(i in 1:length(myterms)){
	if(inPB[i]==1){
		bloodres[i,1]<-which(siggo.pb$Term==myterms[i])
		bloodres[i,2]<-which(siggo.pb$Term==myterms[i])/nrow(siggo.pb)
		bloodres[i,3]<-siggo.pb$P.DE[which(siggo.pb$Term==myterms[i])]
	}
}
bloodres<-data.frame(bloodres)
colnames(bloodres)<-paste0("Blood.",c("Rank","NormRank","pval"))

cordres<-matrix(NA,length(myterms),3) #RANK #Scaled rank #p-value
for(i in 1:length(myterms)){
	if(inCB[i]==1){
		cordres[i,1]<-which(siggo.cb$Term==myterms[i])
		cordres[i,2]<-which(siggo.cb$Term==myterms[i])/nrow(siggo.cb)
		cordres[i,3]<-siggo.cb$P.DE[which(siggo.cb$Term==myterms[i])]
	}
}
cordres<-data.frame(cordres)
colnames(cordres)<-paste0("Cord.",c("Rank","NormRank","pval"))

brainres<-matrix(NA,length(myterms),3) #RANK #Scaled rank #p-value
for(i in 1:length(myterms)){
	if(inFB[i]==1){
		brainres[i,1]<-which(siggo.fb$Term==myterms[i])
		brainres[i,2]<-which(siggo.fb$Term==myterms[i])/nrow(siggo.fb)
		brainres[i,3]<-siggo.fb$P.DE[which(siggo.fb$Term==myterms[i])]
	}
}
brainres<-data.frame(brainres)
colnames(brainres)<-paste0("Brain.",c("Rank","NormRank","pval"))

#Finally, create Table 3 for the manuscript. Order first by what overlap group they are in
#(i.e. all 3 tissues, cord and brain, etc.). Then within those groupings, order by highest
#sum of ranks in each tissue. 
total<-cbind(bloodres,cordres,brainres)
total<-cbind(inPB,inCB,inFB,total)
total$OrderingCol<-0
total$OrderingCol[which(total$inPB==1 & total$inCB==1 & total$inFB==1)]<-1
total$OrderingCol[which(total$inPB==1 & total$inCB==0 & total$inFB==1)]<-2
total$OrderingCol[which(total$inPB==0 & total$inCB==1 & total$inFB==1)]<-3
total$OrderingCol[which(total$inPB==1 & total$inCB==1 & total$inFB==0)]<-4

total$SumRank<-rowSums(cbind(total$Blood.NormRank,total$Cord.NormRank,total$Brain.NormRank),na.rm=TRUE)
total$Term<-myterms

total<-total[order(total$OrderingCol,total$SumRank),]

forpaper<-total[,c("Term","Blood.NormRank","Cord.NormRank","Brain.NormRank")]
forpaper[,2]<-round(forpaper[,2],digits=2)
forpaper[,3]<-round(forpaper[,3],digits=2)
forpaper[,4]<-round(forpaper[,4],digits=2)
#TABLE 3. 
write.csv(forpaper,file="CombinedTissueList.csv",row.names=FALSE,quote=FALSE)

#Write out csv files for post revigo terms for each tissue. 
writepb<-siggo.pb[,c("Term","P.DE")]
writepb$Rank<-seq(1,nrow(writepb))
writepb$ScaledRank<-round(seq(1,nrow(writepb))/nrow(writepb),digits=2)
#Supplemental Table 2
write.csv(writepb,file="SuppTable_Blood.csv",row.names=TRUE,quote=FALSE)

writecb<-siggo.cb[,c("Term","P.DE")]
writecb$Rank<-seq(1,nrow(writecb))
writecb$ScaledRank<-round(seq(1,nrow(writecb))/nrow(writecb),digits=2)
#Supplemental Table 3
write.csv(writecb,file="SuppTable_Cord.csv",row.names=TRUE,quote=FALSE)

writefb<-siggo.fb[,c("Term","P.DE")]
writefb$Rank<-seq(1,nrow(writefb))
writefb$ScaledRank<-round(seq(1,nrow(writefb))/nrow(writefb),digits=2)
#Supplemental Table 4
write.csv(writefb,file="SuppTable_Brain.csv",row.names=TRUE,quote=FALSE)
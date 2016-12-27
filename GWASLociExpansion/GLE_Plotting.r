#In this script we construct the pieces of the GWAS expansion figures seen in Figure 1. 
#There are in total plots that go into each figure.
	#1: Linkage disequilibrium plot
	#2: RefSeq Gene tracker
	#3: Peripheral blood meQTL tracker
	#4: Peripheral blood CpG site tracking line
	#5: Peripheral blood CpG-CpG correlation plot
	#6: Cord blood meQTL tracker
	#7: Cord blood CpG site tracking line
	#8: Cord blood CpG-CpG correlation plot
	#9: Fetal brain meQTL tracker
	#10: Fetal brain CpG site tracking line
	#11: Fetal brain CpG-CpG correlation plot
#These peices can then be arranged on top of one another in Illustrator, PowerPoint, etc.

library(minfi)

#All of the ASD loci and their info
load("/GLE/locustable.rda")
	#collectplus: list of ASD hits
	#adgenelist: lists of what genes are tagged by different asd locibased
	#totalmappings: list of each tissue and what cg positions are tagged by the asd snps there (if any)
	#totalgenes: list of what genes are tagged by the cg sites controlled by the snps. same dimensions as asdgenelist
load("collectplus.rda")
load("asdgenelist.rda")
load("/GLE/totalmappings.rda")
load("/GLE/totalgenes.rda")

#Initialize the gene info so we can have the coordinates for later
load("refFlat.clean.new.order.rda") #Produced in the script /GOBasedEnrichment/AddtlObjects.r

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
bloodbeta<-B #We call this 'bloodbeta' to distinguish between cord, brain, etc. 

#Also load and define beta (% methylation) for cord blood and fetal brain matrix:
cordbeta #Beta matrix for cord blood
brainbeta #Beta matrix for fetal brain data

#Load fetal brain meQTL information
#Load fetal brain meQLT data. 
windowsize<-1000000
fetalbrain<-read.csv("All_Imputed_BonfSignificant_mQTLs.csv",header=TRUE)
fetalbrain<-fetalbrain[which(fetalbrain$SNP_Chr==fetalbrain$DNAm_CHR),]
fetalbrain$Space<-abs(fetalbrain$SNP_BP-fetalbrain$DNAm_BP)
fetalbrain<-fetalbrain[which(fetalbrain$Space<=windowsize),]

#Initialize some things:
plinkex<-"/plink" #Location of plink executable for r^2 estimation. 
genotypedirectory<-"" #Location of genotype files to be used in LD plot (we used 1000G data)
meQTLfilelist<-c("/SEED_meQTLs/","/EARLI_meQTLs/","fetalbrain")
pvaluecutoffs<-c(1.1E-5,2.7E-6,1E-8) #p-value cutoffs for 
library(sqldf)

for (locus in 1:8){ # We will loop through the locustable object that we created in 'GLE_Annotation.r'
	###########################################################################
	##########################PLOT 1: LD Plot##################################
	###########################################################################	
	#Define start and stop coordinates for the total genomic space. 
	chr<-strsplit(locustable$FullLocus,":")[[locus]][[1]]
	chr<-as.numeric(strsplit(chr,"chr")[[1]][[2]])
	begin<-strsplit(locustable$FullLocus,":")[[locus]][[2]]
	End<-begin
	begin<-as.numeric(strsplit(begin,"-")[[1]][[1]])-10
	End<-as.numeric(strsplit(End,"-")[[1]][[2]])+10
	
	#Point to list of SNP names for genotype data being used. We need to find that SNPs in the data itself 
	#that fit into the window space that we have just defined. 
	SNPgeno<-read.table(paste("snplist_",chr,".txt",sep=""),header=FALSE,stringsAsFactors=FALSE)
	snppos<-SNPgeno[,2]
	minpos<-SNPgeno[min(which(SNPgeno[,2]>=begin)),1]
	maxpos<-SNPgeno[max(which(SNPgeno[,2]<=End)),1]
	
	#Read in genotype data (here we use 1000G). 
	seedgeno<-paste(genotypedirectory,"1000G_Recode_",chr,sep="")
	ldout<-paste("ldreport_",locus,sep="") #Nmae of LD output file. 
	#The comparisons argument decides how many SNP to SNP r^2 values that PLINK will generate. 
	#For comparsions = 150 this means that we will get the r^2 info for 150 SNPs downsteam (i.e. greater position)
	#than the SNP. There are 2 options here, one in which comparions is defined as a function of the window size
	#but also a set value. Note that this value can also play a big role in the size of 
	#the LD output file produced. 
	comparisons<-ceiling(1*(max(which(SNPgeno[,2]<=End))-min(which(SNPgeno[,2]>=begin))))
	comparisons<-150
	
	#Use a system command to invoke plink and make the LD output file. 
	system(paste(plinkex," --tfile ",seedgeno," --from ",minpos," --to ",maxpos," --r2 --ld-window ",comparisons," --ld-window-r2 0 --out ",ldout,sep=""))
	
	#Read in LD file
	myLD<-read.table(paste(ldout,".ld",sep=""),header=TRUE)
	LDpos<-sort(unique(c(myLD$BP_A,myLD$BP_B)))
	
	#Here we populate a matrix with the r^2 values; it would be symmetric but we will only be using 1 half (triangle) in the ultimate plot. 
	message("Computing LD matrix")
	LDmatrix<-matrix(NA,length(LDpos),length(LDpos))
	for (dore in 1:(length(LDpos)-1)){
		message("Filling in row ",dore," of total ",length(LDpos))
		nothing<-myLD[which((myLD$BP_A)==LDpos[dore]),]
		LDmatrix[dore,match(nothing$BP_B,LDpos)]<-nothing$R2
	}
	LDmatrix[lower.tri(LDmatrix)]<-t(LDmatrix)[lower.tri(LDmatrix)]
	LDmatrix[upper.tri(LDmatrix)]<-NA
	diag(LDmatrix)<-NA
	
	#Use heatmap.2 to make the the LD block plot. 
	library(gplots)
	png(filename=paste("Locus",locus,"/LDplot_locus",locus,".png",sep=""),width=10, height=4,units="in",res=480,type="cairo",bg="transparent")
	my_palette<-colorRampPalette(c("gray97","red"))
	heatmap.2(LDmatrix,Rowv=FALSE, Colv=FALSE,key=FALSE,dendrogram="none",labRow=FALSE,labCol=FALSE,
		scale="none",trace="none",col=my_palette,density.info="none",margins=c(2,2))
	dev.off()

	###########################################################################
	##########################PLOT 2: RefSeq Gene Tracker######################
	###########################################################################	
	message("Making gene track")
	#Find what genes have start and stop coordinates that live in the defined space.
	chrgenes<-refFlat.clean.new.order[which(refFlat.clean.new.order$CHR==chr),]
	which(as.numeric(chrgenes$txStart)<End & as.numeric(chrgenes$txEnd)>begin)
	chrgenes<-chrgenes[which(as.numeric(chrgenes$txStart)<End & as.numeric(chrgenes$txEnd)>begin),]
	chrgenes$txStart<-ifelse(as.numeric(chrgenes$txStart)<begin,begin,chrgenes$txStart)
	chrgenes$txEnd<-ifelse(as.numeric(chrgenes$txEnd)>End,End,chrgenes$txEnd)
	
	#Define intervals for tracker. We are using a fill plotting device to make these plots
	#and we want to create bins based on position values. tofill is the matrix that will 
	#eventually be used in plotting. We want to place zeroes in the matrix that is defined by
	#position (x-axis) and the gene (y-axis). 
	genespacing<-(1*nrow(chrgenes)+(nrow(chrgenes)+1))
	xrep<-c(seq(begin,End,by=10))	
	xrep[length(xrep)]<-End
	xax<-rep(xrep,genespacing)
	yax<-rep(1:genespacing,each=length(xrep))
	tofill<-matrix(0,length(xrep),genespacing)
	
	#Here is where we populate to fill according to where (x -axis) certain genes (y-axis)
	#live within the space we have defined. 
	startseqs<-seq(2,genespacing,by=2)
	for (michael in 1:nrow(chrgenes)){
		thisgenestart<-max(which(xrep<=chrgenes$txStart[michael]))
		thisgenestop<-min(which(xrep>=chrgenes$txEnd[michael]))
		tofill[thisgenestart:thisgenestop,c(startseqs[michael])]<-1
	}
	fillattach<-c(tofill)
	
	#Format for ggplot2 and define colors. Genes will be denoted in black. 
	toplot<-data.frame(xax,yax,fillattach)
	polycolors<-c("gray97","black")
	
	#Make the plot
	library(ggplot2)
	library(gridExtra)
	library(grid)
	png(filename=paste("/Locus",locus,"/GeneTrack_locus",locus,".png",sep=""),width=10, height=4,units="in",res=480,bg="transparent", type="cairo")
	p5<-qplot(xax,yax,fill=factor(fillattach),data=toplot,geom="tile")+
		theme(panel.background=element_blank(),legend.position="none")+
		scale_fill_manual(values=polycolors)	
	print(p5, vp=viewport(width=unit(6, "inches"),height=unit(2, "inches")))
	dev.off()

	###########################################################################
	######################PLOTS 3,6,9: meQTL tracker###########################
	###########################################################################	
	#First we grab the meQTL info for SNPs that are not ASD-associated. 
	totalcgmaps<-list()
	for (y in 1:3){
		message("Mapping SNP hits for ",y)
		if(y<3){ #First for peripheral blood and cord blood
			thischr<-grep(paste("restotal",chr,"_",sep=""),list.files(meQTLfilelist[y]),fixed=TRUE)
			myres<-lapply(1:length(thischr),function(x){load(paste(meQTLfilelist[y],list.files(meQTLfilelist[y])[thischr[x]],sep="")) 
				ret<-ret[which(ret$pvalue<pvaluecutoffs[y]),] 
				ret$cgpos<-pos[match(ret$gene,rownames(object))]
				#Limit to meQTLs where the SNP and CpG of the meQTL are located within the defined window space. 
				ret<-ret[-which(ret$cgpos<begin | ret$cgpos>End),]
				ret<-ret[-which(as.numeric(as.character((ret$snps)))<begin | as.numeric(as.character((ret$snps)))>End),]
				iteration<-as.list(rep(NA,(length(xrep)-1)))
				if(nrow(ret)>0){
					#Determine where thoe SNP positions lie wrt to the bins created previously. Retain CpG position info. 
					ret$bins<-as.numeric(cut(as.numeric(as.character((ret$snps))),xrep,include.lowest=TRUE))
					iteration<-lapply(1:(length(xrep)-1),function(please){
						if(length(which(ret$bins==please))>0){out<-ret$cgpos[which(ret$bins==please)]
						} else {out<-NA}
						return(out)
					})
				}
				return(iteration)
				}
			)
			#Format and retain this CpG position information. 
			compile<-list()
			for (j in 1:(length(xrep)-1)){
				compile[[j]]<-Filter(Negate(is.na),unlist(lapply(myres,function(x){return(x[[j]])})))		
			}
			mymax<-max(unlist(lapply(compile,length)))
			finished<-do.call("rbind",lapply(compile,function(x){length(x)<-mymax; return(x)}))
			totalcgmaps[[y]]<-finished	
		} else{
		#FIX
			thisbrain<-fetalbrain[which(fetalbrain$SNP_Chr==chr),]
			#Limit to meQTLs where the SNP and CpG of the meQTL are located within the defined window space. 
			thisbrain<-thisbrain[-which(as.numeric(thisbrain$DNAm_BP)<begin | as.numeric(thisbrain$DNAm_BP)>End),]
			thisbrain<-thisbrain[which(thisbrain$SNP_BP>=begin & thisbrain$SNP_BP<=End),]
			#Determine where thoe SNP positions lie wrt to the bins created previously. Retain CpG position info.
			thisbrain$bins<-as.numeric(cut(thisbrain$SNP_BP,xrep,include.lowest=TRUE))
			iteration<-lapply(1:(length(xrep)-1),function(please){
						if(length(which(thisbrain$bins==please))>0){out<-thisbrain$DNAm_BP[which(thisbrain$bins==please)]
						} else {out<-NA}
						return(out)
			})
			#Format and retain this CpG position information.
			mymax<-max(unlist(lapply(iteration,length)))
			finished<-do.call("rbind",lapply(iteration,function(x){length(x)<-mymax; return(x)}))
			totalcgmaps[[y]]<-finished
		}
	}
	#Next we grab the meQTL info for SNPs that are ASD-associated. 
	#We use the collectplus object that we loaded earlier to tell us what the ASD-associated SNPs are. 
	#We places those in the same bins that define our SNP tracker. 
	snpstarter<-strsplit(locustable$Locus,":")[[locus]][[2]]
	snpend<-snpstarter
	snpstarter<-as.numeric(strsplit(snpstarter,"-")[[1]][[1]])
	snpend<-as.numeric(strsplit(snpend,"-")[[1]][[2]])
	snpstarter<-which(collectplus$BP==snpstarter & collectplus$CHR==chr)
	snpend<-which(collectplus$BP==snpend & collectplus$CHR==chr)
	snpstotrack<-collectplus$BP[snpstarter:snpend]
	snpstotrack.bin<-as.numeric(cut(snpstotrack,xrep,include.lowest=TRUE))
	totalcgmaps.asd<-list()
	for (y in 1:3){
		message("Mapping SNP hits for ",y)
		if(y<3){#For peripheral blood and cord blood
			thischr<-grep(paste("restotal",chr,"_",sep=""),list.files(meQTLfilelist[y]),fixed=TRUE)
			#Go theough and retain CpG position information of meQTL targets
			myres<-lapply(1:length(thischr),function(x){load(paste(meQTLfilelist[y],list.files(meQTLfilelist[y])[thischr[x]],sep=""))
				ret<-ret[which(ret$pvalue<pvaluecutoffs[y]),] 
				ret$cgpos<-pos[match(ret$gene,rownames(object))]
				iteration<-as.list(rep(NA,length(snpstotrack)))
				if(nrow(ret)>0){
					iteration<-lapply(snpstotrack,function(please){
						if(length(which(ret$snps==please))>0){out<-ret$cgpos[which(ret$snps==please)]
						} else {out<-NA}
						return(out)
					})
				}
				return(iteration)
				}
			)
			compile<-list()
			for (j in 1:length(snpstotrack)){
				compile[[j]]<-Filter(Negate(is.na),unlist(lapply(myres,function(x){return(x[[j]])})))		
			}
			mymax<-max(unlist(lapply(compile,length)))
			finished<-do.call("rbind",lapply(compile,function(x){length(x)<-mymax; return(x)}))
			totalcgmaps.asd[[y]]<-finished	
		} else{
			#Go through and retain CpG position information of meQTL targets
			thisbrain<-fetalbrain[which(fetalbrain$SNP_Chr==chr),]
			iteration<-lapply(snpstotrack,function(please){
						if(length(which(thisbrain$SNP_BP==please))>0){out<-thisbrain$DNAm_BP[which(thisbrain$SNP_BP==please)]
						} else {out<-NA}
						return(out)
			})
			mymax<-max(unlist(lapply(iteration,length)))
			finished<-do.call("rbind",lapply(iteration,function(x){length(x)<-mymax; return(x)}))
			totalcgmaps.asd[[y]]<-finished
		}
	}	
	#Just as we have made a tracker for SNP positions, make one for CpG positions. 
	#Wehn we draw lines to define meQTLs, the will go from the bin of the SNP tracker that the meQTL is in (bottom of figure)
	#to the bin of the CpG tracker that the meQTL target is in (top of figure). 
	totalcpgpos<-data.frame(pos[which(chrnames==paste("chr",chr,sep=""))])
	colnames(totalcpgpos)<-"OriginalPos"
	totalcpgpos$MatchedBin<-as.numeric(cut(totalcpgpos$OriginalPos,xrep,include.lowest=TRUE))
	
	#Finally we will make the figures themselves depicting the meQTLs. 
	#If the meQTL is not ASD related, it will be colored grey. 
	#If the meQTL is ASD related, it will be colored red. 
	#In total there will be 3 figures (png files) output here. 
	for (y in 1:3){
	png(filename=paste("/Locus",locus,"/meQTLLines_locus",locus,"_tissue_",y,".png",sep=""),width=10, height=4,units="in",res=480,bg="transparent", type="cairo")
	message("Plotting SNP hits for ",y)
		#First we grab the information for the lines for SNPs not related to ASD. 
		regularobject<-totalcgmaps[[y]]
		regularobject<-lapply(1:nrow(regularobject),function(u){matters<-regularobject[u,]
			if(sum(is.na(matters))<length(matters)){matters<-matters[!is.na(matters)]
				matters<-totalcpgpos$MatchedBin[match(matters,totalcpgpos$OriginalPos)]
				out<-cbind(rep(u,length(matters)),matters)} else {out<-cbind(NA,NA)}
			return(out)
			}
		)
		regularobject<-do.call("rbind",regularobject)
		regularobject<-regularobject[which(rowSums(!is.na(regularobject))==2),]
		if(!is.null(nrow(regularobject))){
			regularobject<-cbind(regularobject,rep(0,nrow(regularobject)))} else{regularobject<-c(regularobject,0)}
		if(is.null(nrow(regularobject))){
			temp<-matrix(0,1,3)
			temp[1,]<-regularobject
			regularobject<-temp
		}
		#Next we grab the lines for the SNPs that are related to ASD. 	
		asdobject<-totalcgmaps.asd[[y]]
		rownames(asdobject)<-snpstotrack.bin
		asdobject<-lapply(1:nrow(asdobject),function(u){matters<-asdobject[u,]
			if(sum(is.na(matters))<length(matters)){matters<-matters[!is.na(matters)]
				matters<-totalcpgpos$MatchedBin[match(matters,totalcpgpos$OriginalPos)]
				out<-cbind(rep(as.numeric(rownames(asdobject)[u]),length(matters)),matters)} else {out<-cbind(NA,NA)}
			return(out)
			}
		)
		asdobject<-do.call("rbind",asdobject)
		if(nrow(asdobject)>1){asdobject<-asdobject[which(rowSums(!is.na(asdobject))==2),]}
		if(!is.null(nrow(asdobject))){
			asdobject<-cbind(asdobject,rep(1,nrow(asdobject)))} else{asdobject<-c(asdobject,1)}
		if(is.null(nrow(asdobject))){
			temp<-matrix(0,1,3)
			temp[1,]<-asdobject
			asdobject<-temp
		}
		
		#Finally we use the lines function to make the meQTL lines. Note where we define different
		#colors for different types of lines. 
		forlines<-rbind(regularobject,asdobject)
		maxbin<-as.numeric(cut(max(xrep),xrep),include.lowest=T)
		xrange<-c(0,0,maxbin,maxbin)
		yrange<-c(0,10,0,10)
		bull<-plot(xrange,yrange,pch=NA)
		for (jurassic in 1:nrow(forlines)){
			if (forlines[jurassic,3]==1){
			lines(c(forlines[jurassic,1],forlines[jurassic,2]),c(1,9),col="red",lwd=3)} else{
				lines(c(forlines[jurassic,1],forlines[jurassic,2]),c(1,9),col="gray",lwd=0.2)}
		}
	dev.off()
	}

	###########################################################################
	################PLOTS 4,7,10: CpG position tracking lines##################
	###########################################################################	
	#Next we want to be able to draw lines that go from the CpG (bin) position in the CpG track
	#from the meQTL line figures to the position of the CpG in the CpG-CpG correlation matrix figure. 
	#This is necessary because the correlation matrix figure exists in a site x site space rather
	#than anything having to do with genomic position. So we just draw lines pointing to equally spaced
	#bins of number equal to the total count of CpG sites we considered in this space. It is possible
	#to have some (slight) variation by tissue here despite all studies using the 450k, because
	#different CpGs passed QC and were used in each study (though obviously three is a high degree of overlap). 
	#Note: Technically, a similar figure is necessary to map from the SNPs in the meQTL line figure
	#to position in the LD matrix. However, the density of the SNP data is such that such lines
	#are no useful or necessary for most loci so we omit them. Similar code can be used for that
	#purpose, however. 
	
	#We used the CpG IDs from the beta matrices to extract the right number of CpG sites and their positions. 
	betalist<-list(bloodbeta,cordbeta,brainbeta)
	par(mfrow=c(3,1))
	for (y in 1:3){
	png(filename=paste("/Locus",locus,"/CpGLines_locus",locus,"_tissue_",y,".png",sep=""),width=10, height=4,units="in",res=480,bg="transparent", type="cairo")
		message("CG lines plot for tissue ",y)
		totalcpgpos<-(pos[which(chrnames==paste("chr",chr,sep=""))])
		totalcpgpos<-totalcpgpos[which(totalcpgpos>=begin & totalcpgpos <= End)]
		tissue<-betalist[[y]]
		tissue<-tissue[which(rowSums(is.na(tissue))!=ncol(tissue)),]
		tissuepos<-pos[match(rownames(tissue),rownames(object))]
		tissue<-tissue[match(totalcpgpos,tissuepos),]
		tissuepos<-tissuepos[match(totalcpgpos,tissuepos)]
		
		#Bottom tracker is the CpG position and should emerge from the positions on the CpG track. 
		bottomtracker<-as.numeric(cut(tissuepos,xrep,include.lowest=TRUE))
		#toptracker is where the lines go to and are meant to line up with the correlation matrix. 
		#This should be simple equidistant locations, as explained above. 
		toptracker<-(seq(1:length(tissuepos))/length(tissuepos))
		toptracker<-round(toptracker*length(xrep))
		xrange<-c(0,0,maxbin,maxbin)
		yrange<-c(0,10,0,10)
		bull<-plot(xrange,yrange,pch=NA)
		#Finally we draw the lines. 
		for (jurassic in 1:length(toptracker)){
			lines(c(bottomtracker[jurassic],toptracker[jurassic]),c(1,9),col="black",lwd=1)
		}
	dev.off()
	}
	
	###########################################################################
	################PLOTS 5,8,11: CpG-CpG correlation plot#####################
	###########################################################################	
	#These plots are constructed very similarly to the LD plots (Plot 1). Instead of
	#deriving a correlation metric from PLINK, we have to calculate the correlation of
	#CpG sites. We do this and directly populate the output in matrix form that can 
	#be readily used in a heatmap.2 context. 
	
	library(gplots)
	for (y in 1:3){
		#First collect position information to know how big our initialized matrix in which
		#we put our correlation infromation should be. 
		png(filename=paste("/Locus",locus,"/CGcorr_locus",locus,"_tissue_",y,".png",sep=""),width=10, height=4,units="in",res=480,bg="transparent", type="cairo")
		message("CG corr plot for tissue ",y)
			totalcpgpos<-(pos[which(chrnames==paste("chr",chr,sep=""))])
			totalcpgpos<-totalcpgpos[which(totalcpgpos>=begin & totalcpgpos <= End)]
			tissue<-betalist[[y]]
			tissue<-tissue[which(rowSums(is.na(tissue))!=ncol(tissue)),]
			tissuepos<-pos[match(rownames(tissue),rownames(object))]
			tissue<-tissue[match(totalcpgpos,tissuepos),]
			tissuepos<-tissuepos[match(totalcpgpos,tissuepos)]
		
		#Initialize matrix and populate it. Again it is symmetric but we only need only half (triangle)
		#for plotting purposes. 
		CPGmatrix<-matrix(NA,length(tissuepos),length(tissuepos))
		for (dore in 1:(length(tissuepos)-1)){
			for (joan in dore:length(tissuepos)){
				if(abs(dore-joan)<=100){
				CPGmatrix[dore,joan]<-abs(cor(tissue[dore,],tissue[joan,]))}
			}
		}
		diag(CPGmatrix)<-NA
		print(range(CPGmatrix,na.rm=TRUE))
		
		#Finally, construct the heatmap. 
		my_palette<-colorRampPalette(c("gray97","red"))
		mybreaks<-seq(0,1,by=0.01)
		heatmap.2(CPGmatrix,Rowv=NULL, Colv=NULL,key=TRUE,dendrogram="none",labRow=FALSE,labCol=FALSE,
		scale="none",trace="none",col=my_palette,breaks=mybreaks,density.info="none",margins=c(2,2))
		dev.off()
	}
	
}
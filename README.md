# ASDmeQTL_Manuscript
Scripts for manuscript "Cross tissue integration of genetic and epigenetic data offers insight into autism spectrum disorder"

Now available here: http://biorxiv.org/content/early/2016/12/06/091330

Directory descriptions:
	meQTLParameters: Implement data-driven procedure to select meQTL query parameters 
		(window size, MAF threshold, DNAM sd threshold) to ensure adequate power for SNP & DNAm data
		
	meQTLScreening&FDR: Perform meQTL query and implement procedure to control meQTL FDR 
	
	SNPBasedEnrichment: Perform analysis to examine if ASD-related SNPs are enriched for meQTLs, across tissue type
	
	GOBasedEnrichment: Examine which biological pathways are over-represented by meQTL targets of 
		ASD SNPs and meQTL targets generally, across tissue type
		
	GWASLociExpansion: Demonstrate utility of meQTL target information to enhance regions implicated 
		by GWAS SNP positions
		
	RegulatoryFeatureEnrichment: Determine which DNaseI HS sites, chromatin marks, and TFBSs that 
		meQTL targets of pyschiatric disorder-related SNPs, and meQTL targets generally, preferentially overlap


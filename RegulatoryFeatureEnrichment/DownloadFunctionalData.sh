###################################
### Download the data from UCSC ###
###################################
## TFBS data
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV2.bed.gz ./TFBS.gz

## DNaseIHS data
#CD14_monocytes_Pk_DNaseIHS
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeOpenChromDnase/wgEncodeOpenChromDnaseMonocd14Pk.narrowPeak.gz ./CD14_monocytes_Pk_DNaseIHS.gz

#CD4_Pk_DNaseIHS
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/wgEncodeUwDnaseCd4naivewb11970640PkRep1.narrowPeak.gz ./CD4_Pk_DNaseIHS.gz

#Cerebellum_Pk_DNaseIHS
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeOpenChromDnase/wgEncodeOpenChromDnaseCerebellumocPk.narrowPeak.gz ./Cerebellum_Pk_DNaseIHS.gz

#Cerebrum_Pk_DNaseIHS
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeOpenChromDnase/wgEncodeOpenChromDnaseCerebrumfrontalocPk.narrowPeak.gz ./Cerebrum_Pk_DNaseIHS.gz

#FrCortex_Pk_DNaseIHS
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeOpenChromDnase/wgEncodeOpenChromDnaseFrontalcortexocPk.narrowPeak.gz ./FrCortex_Pk_DNaseIHS.gz

#IMR90_Pk_DNaseIHS
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeOpenChromDnase/wgEncodeOpenChromDnaseImr90Pk.narrowPeak.gz ./IMR90_Pk_DNaseIHS.gz

###############################################
### Download from Epigenome Roadmap Project ###
###############################################
### For samples that rmare of European Ancestry (TC010, TC015), if there was data from more than one individual available, both were downloaded. Data were then combined. Data from non-European Ancestry background (ie TC012, TC016) were not included for analysis
### For fetal samples, data from detal days 120 and 122 were included (Broad samples, rather than UCSF).

#H3K27ac_Blood
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1127nnn/GSM1127145/suppl/GSM1127145_UCSF-UBC.Peripheral_Blood_Mononuclear_Primary_Cells.H3K27ac.TC015.bed.gz -O H3K27ac_Blood.gz

#H3K27ac_Lung
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM906nnn/GSM906395/suppl/GSM906395_UCSD.Lung.H3K27ac.STL002.bed.gz -O H3K27ac_Lung_STL002.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1013nnn/GSM1013123/suppl/GSM1013123_UCSD.Lung.H3K27ac.STL001.bed.gz -O H3K27ac_Lung_STL001.gz

#H3K9ac_Blood
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM613nnn/GSM613879/suppl/GSM613879_UCSF-UBC.Peripheral_Blood_Mononuclear_Primary_Cells.H3K9ac.TC010.HS2623.bed.gz -O H3K9ac_Blood.gz

#H3K27me3_Blood
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM613nnn/GSM613877/suppl/GSM613877_UCSF-UBC.Peripheral_Blood_Mononuclear_Primary_Cells.H3K27me3.TC010.HS2621.bed.gz -O H3K27me3_Blood_TC010.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1127nnn/GSM1127130/suppl/GSM1127130_UCSF-UBC.Peripheral_Blood_Mononuclear_Primary_Cells.H3K27me3.TC015.bed.gz -O H3K27me3_Blood_TC015.gz

#H3K27me3_FetalBrain
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM621nnn/GSM621393/suppl/GSM621393_BI.Fetal_Brain.H3K27me3.UW_H-22510.bed.gz -O H3K27me3_FetalBrain_GSM621393.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM916nnn/GSM916061/suppl/GSM916061_BI.Fetal_Brain.H3K27me3.UW_H22676.bed.gz -O H3K27me3_FetalBrain_GSM916061.gz


#H3K36me3_Blood
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM613nnn/GSM613880/suppl/GSM613880_UCSF-UBC.Peripheral_Blood_Mononuclear_Primary_Cells.H3K36me3.TC010.HS2624.bed.gz -O H3K36me3_Blood_TC010.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1127nnn/GSM1127131/suppl/GSM1127131_UCSF-UBC.Peripheral_Blood_Mononuclear_Primary_Cells.H3K36me3.TC015.bed.gz -O /H3K36me3_Blood_TC015.gz

#H3K36me3_FetalBrain
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM621nnn/GSM621410/suppl/GSM621410_BI.Fetal_Brain.H3K36me3.UW_H-22510.bed.gz -O H3K36me3_FetalBrain.gz

#H3K36me3_Lung
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM956nnn/GSM956014/suppl/GSM956014_UCSD.Lung.H3K36me3.STL002.bed.gz -O H3K36me3_Lung_STL002.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1059nnn/GSM1059437/suppl/GSM1059437_UCSD.Lung.H3K36me3.STL001.bed.gz -O H3K36me3_Lung_STL001.gz

#H3K4me1_Blood
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1127nnn/GSM1127143/suppl/GSM1127143_UCSF-UBC.Peripheral_Blood_Mononuclear_Primary_Cells.H3K4me1.TC015.bed.gz -O H3K4me1_Blood.gz

#H3K4me1_FetalBrain
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM706nnn/GSM706850/suppl/GSM706850_BI.Fetal_Brain.H3K4me1.UW_H22676.bed.gz -O H3K4me1_FetalBrain.gz

#H3K4me1_Lung
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM910nnn/GSM910572/suppl/GSM910572_UCSD.Lung.H3K4me1.STL002.bed.gz -O H3K4me1_Lung_STL002.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1059nnn/GSM1059443/suppl/GSM1059443_UCSD.Lung.H3K4me1.STlhL001.bed.gz -O H3K4me1_Lung_STL001.gz

#H3K4me3_Blood
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1127nnn/GSM1127126/suppl/GSM1127126_UCSF-UBC.Peripheral_Blood_Mononuclear_Primary_Cells.H3K4me3.TC015.bed.gz -O H3K4me3_Blood.gz

#H3K4me3_FetalBrain
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM621nnn/GSM621457/suppl/GSM621457_BI.Fetal_Brain.H3K4me3.UW_H-22510.bed.gz -O H3K4me3_FetalBrain

#H3K4me3_Lung
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM915nnn/GSM915336/suppl/GSM915336_UCSD.Lung.H3K4me3.STL002.bed.gz -O H3K4me3_Lung.gz

#H3K9me3_Blood
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM613nnn/GSM613878/suppl/GSM613878_UCSF-UBC.Peripheral_Blood_Mononuclear_Primary_Cells.H3K9me3.TC010.HS2622.bed.gz -O H3K9me3_Blood_TC010.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1127nnn/GSM1127133/suppl/GSM1127133_UCSF-UBC.Peripheral_Blood_Mononuclear_Primary_Cells.H3K9me3.TC015.bed.gz -O H3K9me3_Blood_TC015.gz

#H3K9me3_FetalBrain
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM621nnn/GSM621427/suppl/GSM621427_BI.Fetal_Brain.H3K9me3.UW_H-22510.bed.gz -O H3K9me3_FetalBrain_GSM621427.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM916nnn/GSM916054/suppl/GSM916054_BI.Fetal_Brain.H3K9me3.UW_H22676.bed.gz -O H3K9me3_FetalBrain_GSM916054.gz

#H3K9me3_Lung
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM906nnn/GSM906411/suppl/GSM906411_UCSD.Lung.H3K9me3.STL002.bed.gz -O H3K9me3_Lung_STL002.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1120nnn/GSM1120355/suppl/GSM1120355_UCSD.Lung.H3K9me3.STL001.bed.gz -O H3K9me3_Lung_STL001.gz

####################
### Gunzip files ###
####################
gunzip *.gz

###################################################################
### Combine data where more than one individual is contributing ###
###################################################################
## Note this takes some time, so you can break these up into separate jobs if you want to speed things up

#H3K27ac_Lung
bedtools intersect -a H3K27ac_Lung_STL002 -b H3K27ac_Lung_STL001 > H3K27ac_Lung
rm H3K27ac_Lung_STL002
rm H3K27ac_Lung_STL001

#H3K36me3_Lung
bedtools intersect -a H3K36me3_Lung_STL002 -b H3K36me3_Lung_STL001 > H3K36me3_Lung
rm H3K36me3_Lung_STL002
rm H3K36me3_Lung_STL001

#H3K27me3_Blood
bedtools intersect -a H3K27me3_Blood_TC010 -b H3K27me3_Blood_TC015 > H3K27me3_Blood
rm H3K27me3_Blood_TC010
rm H3K27me3_Blood_TC015

#H3K27me3_FetalBrain
bedtools intersect -a H3K27me3_FetalBrain_GSM621393 -b H3K27me3_FetalBrain_GSM916061 > H3K27me3_FetalBrain
rm H3K27me3_FetalBrain_GSM621393
rm H3K27me3_FetalBrain_GSM916061

#H3K36me3_Blood
bedtools intersect -a  H3K36me3_Blood_TC010 -b H3K36me3_Blood_TC015 > H3K36me3_Blood
rm H3K36me3_Blood_TC010
rm H3K36me3_Blood_TC015

#H3K4me1_Lung
bedtools intersect -a H3K4me1_Lung_STL002 -b H3K4me1_Lung_STL001 > H3K4me1_Lung
rm H3K4me1_Lung_STL002
rm H3K4me1_Lung_STL001

#H3K9me3_Blood
bedtools intersect -a H3K9me3_Blood_TC010 -b H3K9me3_Blood_TC015 > H3K9me3_Blood
rm H3K9me3_Blood_TC010
rm H3K9me3_Blood_TC015

#H3K9me3_FetalBrain
bedtools intersect -a H3K9me3_FetalBrain_GSM621427 -b H3K9me3_FetalBrain_GSM916054  > H3K9me3_FetalBrain
rm H3K9me3_FetalBrain_GSM621427
rm H3K9me3_FetalBrain_GSM916054

#H3K9me3_Lung
bedtools intersect -a H3K9me3_Lung_STL002 -b H3K9me3_Lung_STL001 > H3K9me3_Lung
rm H3K9me3_Lung_STL002
rm H3K9me3_Lung_STL001




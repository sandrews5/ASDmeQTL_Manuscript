PLINK= #Location of plink executable. 
ROOT_DIR= #Location of SNP genotypes. 
OUT_DIR= #Where to put overlap GENOTYPE FILES
SUB_DIR= #Location of SNP lists to be pruned genreated in 'PGCoverlap.r'

for CHR in {1..22}
do	
${PLINK} --tfile ${ROOT_DIR}SEED_MethGeno_${CHR} --extract ${SUB_DIR}snplist_${CHR}.txt --recode transpose --out ${OUT_DIR}SEED_Overlap_${CHR}
done

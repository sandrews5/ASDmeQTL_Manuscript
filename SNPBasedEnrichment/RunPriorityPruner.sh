CHR=$SGE_TASK_ID;
PRUNER=PriorityPruner.jar #Location of PriorityPruner executable
PRUNE_DIR= #Where input tables live and where pruned genotype files will be placed. 
OV_DIR= #Where the overlap genotype files (produced in 'PriorityPruner_MakeGenotypes.sh') live

java -jar ${PRUNER} --tfile ${OV_DIR}SEED_Overlap_${CHR} --snp_table ${PRUNE_DIR}inputable_ASD_${CHR}.txt \
	--r2 0.7 --out ${PRUNE_DIR}SEED_Pruned_${CHR}

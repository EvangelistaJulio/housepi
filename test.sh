#!/bin/sh
#SBATCH -J test
#SBATCH -p mmmi
#SBATCH --mem-per-cpu=2000
#SBATCH --cpus-per-task=1
#SBATCH -e ./log/test_%A.err
#SBATCH -o ./log/test_%A.log

module load gcc/11.3.0

WD="/pasteur/zeus/projets/p01/MMMICovid/fcamagli/household-epidemics"
SET="${WD}/default_setup.csv"
DAT="${WD}/data"
OUT="${WD}/output"
EXE="${WD}/program.exe"

mkdir -p ${OUT}

srun ${EXE} ${SET} ${DAT}/database.txt ${DAT}/community_infection.txt ${DAT}/parameters.txt ${OUT}/test.txt ${OUT}/augm.txt

exit 0

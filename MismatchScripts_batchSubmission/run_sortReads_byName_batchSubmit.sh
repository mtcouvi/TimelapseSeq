#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-03:00
#SBATCH -p short
#SBATCH --mem=50G
#SBATCH --mail-type=END
#SBATCH --mail-user=mtcouvi@gmail.com

module load gcc/6.2.0
module load R/4.0.1

prefix=$1
echo $prefix

path=`pwd`/${prefix}

/n/groups/churchman/mc348/TimelapseSeq/Scripts/sortReads_byName.R reads.MM.bed $path


############# For quick batch submission #############
# Libs='G1 G2 G3 G4 G5 H1 H2 H3 H4 H5'
# dirName='MMinformed4_Top100_withDups'
# 
# for lib in $Libs
# do
# sbatch /n/groups/churchman/mc348/TimelapseSeq/Scripts/run_sortReads_byName_batchSubmit.sh ${lib}_${dirName}
# done

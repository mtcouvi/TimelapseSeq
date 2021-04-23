#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-03:00
#SBATCH -p short
#SBATCH --mem=50G
#SBATCH --mail-type=END
#SBATCH --mail-user=mtcouvi@gmail.com

### LAST UPDATED:	02/15/2020 by Brendan Smalec

### USE: 			This script will concatenate the mismatches from paired reads into 
###					one line (total mismatches per fragment) while also removing 
###					mismatches that are duplicated on both reads (aka, they are within
###					the overlapping region). 

###					The MM_per_fragment script keeps the information for read #1 of the
###					pair (coordinates, sequence *edited 10/2020 by MTC to use full fragment as the sequence, etc) while updating the following fields
###					to reflect the mismatch information for both read #1 and read #2:
###					mismatch read nt, mismatch genome, position, mismatch genome nt, 
###					total # of mismatches, # of T>C mismatches, and whether or not the 
###					fragment contains any T>C mismatches (binary call). However this 
###					does NOT update the mismatch read position or the length of the read,
###					so if these parameters will be used for subsequent analysis, the 
###					script will need to be adapted. 

### REQUIREMENTS:	reads.MM.bed_sort.bed (in the same directory)
###					MM_per_fragment.R (R script, in the same directory)

# load required modules
module load gcc/6.2.0
module load samtools/1.9
module load bedtools/2.27.1
module load R/4.0.1

path=`pwd`
prefix=$1

# break into smaller files to run MM_per_fragment in parallel
echo "Splitting up paired reads"
split --lines=20000 ${prefix}/reads.MM.bed_sort.bed ${prefix}/tmp.
echo "Done"


# submit all tmp files for analysis in parallel
echo "Submitting parallel jobs"
temp_files=`ls ${prefix}/tmp*`

for x in $temp_files
do
		
		fname=${x/${prefix}\//}


	
		OUT="${x}_MMperF.sh" 		#THIS is the name of the script to be written. Include var x.
		echo "#!/bin/bash" >> $OUT
		echo "#SBATCH -n 1                              # Request one core" >> $OUT
		echo "#SBATCH -t 480" >> $OUT
		echo "#SBATCH -p short                           # Partition to run in" >> $OUT
		#echo "#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL" >> $OUT
		#echo "#SBATCH --mail-user=x@g.harvard.edu   # Email to which notifications will be sent" >> $OUT
		echo "" >> $OUT				
		echo "module load gcc/6.2.0" >> $OUT
		echo "module load R/4.0.1" >> $OUT
		echo "mkdir ${x}_analysis" >> $OUT
# 		echo "cd ${x}_analysis" >> $OUT
		echo "" >> $OUT
		echo "/n/groups/churchman/mc348/TimelapseSeq/Scripts/MM_per_fragment_20-02-28_MC_batchSubmit.R ${x} $path $fname" >> $OUT
		echo "echo $PWD/${x}" >> $OUT
# 		echo "ls" >> $OUT
		echo "mv ${x}_analysis/${fname}_fragments.bed ${path}/${prefix}/." >> $OUT
# 		echo "cd .." >> $OUT
		echo "rm -r ${x}_analysis" >> $OUT
		#echo "Rscript --vanilla ~/scripts/o2_grnaSearch2.sh ../../../fastaToSearch/seqtosearch_${i}_split${n}.fa" >> $OUT
		#echo "date '+%A %W %Y %X'" >> $OUT
		echo "" >> $OUT
		chmod u+x $OUT						# MUST change permission, otherwise will not submit
		sbatch $OUT							# Submit the script 
	
done

# when this is done running, test that all tmp files have been analyzed by running:
# ./test_MM_per_F.sh
# like with test_submitMM.sh

# once all of these have run, then run the following:
# cat tmp*fragments.bed > reads_MM.fragments.bed

# sort
#echo "Sorting fragments"
#bedtools sort -i reads_MM.fragments.bed > reads_MM.fragments.sort.bed
#echo "Done"
# sbatch -p short -t 0-12:00 --mem=50G --wrap="bedtools sort -i reads_MM.fragments.bed > reads_MM.fragments.sort.bed"
# rm reads_MM.fragments.bed

############# For quick batch submission #############
# Libs='G1 G2 G3 G4 G5 H1 H2 H3 H4 H5' # 'G1 G2 G3 G4 G5 H1 H2 H3 H4 H5'
# dirName='MMinformed4_Top100_withDups'
# 
# for lib in $Libs
# do
# sbatch -e logs/run_MM_per_frag_${lib}_${dirName}.err -o logs/run_MM_per_frag_${lib}_${dirName}.log /n/groups/churchman/mc348/TimelapseSeq/Scripts/run_MM_per_fragment_20_02_28_MC_batchSubmit.sh ${lib}_${dirName}
# done

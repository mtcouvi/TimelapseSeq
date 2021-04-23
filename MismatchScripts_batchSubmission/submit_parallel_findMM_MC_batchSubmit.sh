#!/bin/bash
#SBATCH -n 1 # Request one core
#SBATCH -t 60                             
#SBATCH -p short     

### LAST UPDATED:	10/23/2019 by Brendan Smalec
###					3/15/2021 by Mary Couvillion

### USE:			This script runs findMismatches_complete.R written by me and Kate
###					to output the read position, read nt, genome position, and genome nt
###					of each mismatch within a read. It also calculates the total number of
###					mismatches per read, the number of T>C mismatches within a read,
###					whether or not the read contains any T>C mismatches (binary call), and
###					finally the length of the read. 
###					Note that there has been several iterations of findMismatches_complete
###					and although we now convert reads to pseudo-bed format before running,
###					we used to run it with .sam files. Some of the comments and file 
###					naming may still refer to .sam format, however all code has been 
###					changed to run on .bed format. 

### REQUIREMENTS:	Processed alignments, split up into smaller files with prefix tmp.
###						which should be in the same directory
###					findMismatches_complete_19_10_23.R (R script, in same directory)
###					test_submitMM.sh (script, in the same directory)

prefix=$1
genome=$2


if [ "${genome}" = "Hela" ]
then
fastapath="/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge.fasta"
elif [ "${genome}" = "K562" ]
then
fastapath="/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/K562_ensGRCh38_dm6_ercc_cat.fasta"
elif [ "${genome}" = "mouse" ]
then
fastapath="/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/mouseNIH3T3_mm10_dm6_ercc_cat.fasta"
fi


# create a list of all alignment files             
temp_files=`ls ${prefix}/tmp*`
path=`pwd`
# submit a job to run findMismatches_complete.R on each tmp alignment file
# many lines commented out here are from the script i copied this from and don't want
# to lose the commands, although they are not used here

for x in $temp_files

do
	
		OUT="${x}_mm.sh" 		#THIS is the name of the script to be written. Include var x.
		echo "#!/bin/bash" >> $OUT
		echo "#SBATCH -n 1                              # Request one core" >> $OUT
		echo "#SBATCH -t 240" >> $OUT
		echo "#SBATCH -p short                           # Partition to run in" >> $OUT
# 		echo "#SBATCH -e ${path}/${prefix}/slurm-%j.err" >> $OUT
# 		echo "#SBATCH -o ${path}/${prefix}/slurm-%j.out" >> $OUT
		echo "" >> $OUT				
		echo "module load gcc/6.2.0" >> $OUT
		echo "module load bedtools/2.27.1" >> $OUT
		echo "echo ${x}" >> $OUT
		echo "mkdir ${x}_analysis" >> $OUT
		echo "" >> $OUT
		echo "/n/groups/churchman/mc348/TimelapseSeq/Scripts/findMismatches_complete_20_02_28_MC_variableGenome.R $x $path $fastapath" >> $OUT
		echo "rm -r ${x}_analysis" >> $OUT
		echo "" >> $OUT
		chmod u+x $OUT						# MUST change permission, otherwise will not submit
		sbatch $OUT							# Submit the script 
	
done

# when this is done running, test that all tmp files have been analyzed by running:
# ./test_submitMM.sh 
# this will output a list of tmp files with no corresponding MM files - if there are any,
# replace that list above instead of `ls tmp*` and run again (sometimes adding a longer
# run time lets them all run, but sometimes they just don't finish running during the 
# first batch and I'm not sure why)

# once all tmp files have been analyzed, cat them together:
# cat tmp*MM.bed > reads.MM.bed

############# For quick batch submission #############
# Libs='G1 G2 G3 G4 G5 H1 H2 H3 H4 H5'
# dirName='MMinformed4_Top100_withDups'
# genome='mouse'
# 
# for lib in $Libs
# do
# sbatch -e logs/submit_parallel_${lib}_${dirName}.err -o logs/submit_parallel_${lib}_${dirName}.log /n/groups/churchman/mc348/TimelapseSeq/Scripts/submit_parallel_findMM_MC_batchSubmit.sh ${lib}_${dirName} $genome
# done

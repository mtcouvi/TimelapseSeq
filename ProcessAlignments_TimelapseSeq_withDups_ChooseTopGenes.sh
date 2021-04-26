#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-12:00
#SBATCH -p short
#SBATCH --mem=50G
#SBATCH --mail-type=END
#SBATCH --mail-user=mtcouvi@gmail.com

###

####SBATCH -c 4
####SBATCH -t 0-24:00
####SBATCH -p medium
####SBATCH --mem=50G
####SBATCH --mail-type=END
####SBATCH --mail-user=mtcouvi@gmail.com

### For small, in silico datasets:
### SBATCH -c 4
### SBATCH -t 0-01:00
### SBATCH -p short
### SBATCH --mem=10G
### SBATCH --mail-type=END
### SBATCH --mail-user=mtcouvi@gmail.com

### Modified from Brendan's processAlignments.sh script

### USE: 			This script will remove secondary alignments and duplicates from the 
###					STAR alignment, remove reads that map to >4 locations and paired
###					reads where the mate did not map properly. The CIGAR string of the
###					remaining reads are then processed so that the mismatches within 
###					each read can be called (soft-clipped bases are removed, and 
###					mismatches are distinguished from matches). Finally, reads are 
###					converted from sam to a pseudo-bed format for quicker downstream
###					analysis.

### REQUIREMENTS:	'Aligned.out.sam' (STAR output, in the same directory)
###  				modifyBed.R (script, in the same directory)
###					path to genome fasta, set below:

module load picard/2.8.0

libName=$1
Expt=$2
geneNum=$3
data=$4

# Set path to fasta
if [ "${data}" = "Hela" ]
then
fastapath='/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge.fasta'
elif [ "${data}" = "mouse" ]
then
fastapath='/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/mouseNIH3T3_mm10_dm6_ercc_cat.fasta'
fi



# Downsample bam

# Straight downsample
# samtools view -@ 3 -bs 42.3 ${libName}_${MapMethod}_Aligned.sortedByCoord.noSpike.bam > ${libName}_${MapMethod}_Aligned.sortedByCoord.subsamp.bam

# Top genes
samtools view -@ 3 -b -o ${libName}_Top${geneNum}.bam -L Top${geneNum}_${Expt}.bed ${libName}_Aligned.sortedByCoord.m_noSpike.bam
# 
# # sort by read name (required for fixmate)
samtools sort -@ 3 -n -o ${libName}.nameSort.bam ${libName}_Top${geneNum}.bam  ##Aligned.sortedByCoord.subsamp ${libName}-tot.sort.human.ds.bam
# 
# # Remove reads whose pairs did not make it through the bed filter (this doesn't get rid of the ones that have multiple alignments!)
samtools view -h -@ 3 ${libName}.nameSort.bam | rev | uniq -D -f16 | rev | samtools view -h -b -S > ${libName}.nameSort.paired.bam


# add fixmate m field (required for markdup)
echo add fixmate
/n/groups/churchman/bms36/programs/samtools-1.10/bin/samtools fixmate -@ 3 -m ${libName}.nameSort.paired.bam ${libName}.nameSort.fm.bam
rm ${libName}.nameSort.bam
echo Done

# sort by coordinate (required for markdup)
samtools sort -@ 3 -o ${libName}_2temp.bam ${libName}.nameSort.fm.bam

echo Sorted input
samtools view -@ 3 -c ${libName}_2temp.bam 

# identify PCR and optical duplicates
# -d The optical duplicate distance. Suggested settings of 100 for HiSeq style platforms or about 2500 for NovaSeq ones. Default is 0 to not look for optical duplicates. When set, duplicate reads are tagged with dt:Z:SQ for optical duplicates and dt:Z:LB otherwise. Calculation of distance depends on coordinate data embedded in the read names produced by the Illumina sequencing machines.
# use 100 for Next-seq
# echo Identifying duplicates
# /n/groups/churchman/bms36/programs/samtools-1.10/bin/samtools markdup -s -d 100 -@ 3 ${libName}_${MapMethod}_1temp.bam ${libName}_${MapMethod}_2temp.bam
# echo Done
# rm ${libName}_${MapMethod}_1temp.bam

# report read statistics
echo Running samtools flagstat
samtools flagstat -@ 3 ${libName}_2temp.bam
echo Done

# remove: PCR and optical duplicates, secondary alignments, and supplementary alignments

# not primary
samtools view -@ 3 -F 0x100 -o ${libName}_4temp.bam ${libName}_2temp.bam
rm ${libName}_2temp.bam

echo nonprimary removed
samtools view -@ 3 -c ${libName}_4temp.bam



# PCR or optical duplicate
# samtools view -@ 3 -F 0x400 -o ${libName}_${MapMethod}_4temp.bam ${libName}_${MapMethod}_3temp.bam
# rm ${libName}_${MapMethod}_3temp.bam

# echo PCR/optical duplicates removed
# samtools view -@ 3 -c ${libName}_${MapMethod}_4temp.bam

# supplementary alignment
samtools view -@ 3 -F 0x800 -o ${libName}_5temp.bam ${libName}_4temp.bam
rm ${libName}_4temp.bam

echo supplementary alignments removed
samtools view -@ 3 -c ${libName}_5temp.bam

# reads with unmapped mate
samtools view -@ 3 -F 0x8 -o ${libName}_6temp.bam ${libName}_5temp.bam
rm ${libName}_5temp.bam

echo reads with unmapped mate removed
samtools view -@ 3 -c ${libName}_6temp.bam


# remove reads that map to more than 4 positions, even if it is the primary alignment
samtools view -@ 3 -q 1 -o ${libName}_filtered_temp.bam ${libName}_6temp.bam
rm ${libName}_6temp.bam

echo reads with more than 4 mapping locations removed
samtools view -@ 3 -c ${libName}_filtered_temp.bam 

# Remove unpaired reads again and should get rid of the rest because they should now be unique (note: out of 10 files with 2-5 mil reads, 2 reads got through. See below for command to manually remove them)
# sort by read name (required for uniq)
samtools sort -@ 3 -n -o ${libName}_filtered_temp.nameSort.bam ${libName}_filtered_temp.bam  
remove unique qnames
samtools view -h -@ 3 ${libName}_filtered_temp.nameSort.bam | rev | uniq -D -f19 | rev | samtools view -h -b -S > ${libName}_filtered_temp.paired.nameSort.bam
## To manually remove specific reads:
## samtools view -h -@ 3 ${libName}_filtered_temp.nameSort.bam | rev | uniq -D -f19 | rev | grep -v 'A00794:219:HNH7MDRXX:2:2169:10628:14043\|A00794:219:HNH7MDRXX:2:2169:10628:23281' | samtools view -h -b -S > ${libName}_filtered_temp.paired.nameSort.bam

# sort by coordinate (required for biostar84452.jar)
samtools sort -@ 3 -o ${libName}_filtered_temp.paired.bam ${libName}_filtered_temp.paired.nameSort.bam


######## Downsample to numReads reads ######
# Get readcount in bam
# let "rc=`samtools view -c ${libName}_${MapMethod}_filtered_temp.bam`/2"
# ratio=`echo 'scale=3;'$numReads'/'$rc | bc` # need the bc to do floating-point arithmetic
# # Downsample resulting bam so subsequent steps don't take so long
# samtools view -@ 3 -bs 42${ratio} ${libName}_${MapMethod}_filtered_temp.bam > ${libName}_${MapMethod}_filtered_temp_ds.bam
#################################################

# convert from .bam to .sam
echo converting bam to sam
samtools view -@ 3 -h -O SAM -o ${libName}_filtered_temp.sam ${libName}_filtered_temp.paired.bam
echo Done

# remove soft clipped bases from reads
echo removing soft clipped bases
java -jar /n/groups/churchman/bms36/programs/jvarkit/dist/biostar84452.jar -o ${libName}_filtered.noSoft_temp.sam ${libName}_filtered_temp.sam
rm ${libName}_filtered_temp.sam
echo Done
############ Before running samfixcigar, need to create fasta index and sequence dictionary
############ Load picard: module load picard/2.8.0
############ First, fix fasta so that all lines are the same length (need extra memory for this, I used 50G): java -jar $PICARD/picard-2.8.0.jar NormalizeFasta I=/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/BrendanReads/hg38_dm6_ercc_cat.fasta O=/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/BrendanReads/hg38_dm6_ercc_cat.formatted.fasta LINE_LENGTH=60
############ Next, make fasta index: samtools faidx /n/groups/churchman/mc348/TimelapseSeq/SeqFiles/BrendanReads/hg38_dm6_ercc_cat.formatted.fasta
############ Make picard sequence dictionary: java -jar $PICARD/picard-2.8.0.jar CreateSequenceDictionary R=/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/BrendanReads/hg38_dm6_ercc_cat.formatted.fasta O=/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/BrendanReads/hg38_dm6_ercc_cat.formatted.dict
# converting cigar string to discriminate between mismatches and matches
echo Converting cigar string
java -jar /n/groups/churchman/bms36/programs/jvarkit/dist/samfixcigar.jar -r $fastapath ${libName}_filtered.noSoft_temp.sam > ${libName}_filtered.noSoft.fixCigar_temp.sam
rm ${libName}_filtered.noSoft_temp.sam # /n/groups/churchman/mc348/TimelapseSeq/SeqFiles/BrendanReads/hg38_dm6_ercc_cat.formatted.fasta
echo Done

# get rid of sam header
echo Removing sam header
samtools view -@ 3 -O SAM -o ${libName}_filtered.noSoft.fixCigar.noHead_temp.sam ${libName}_filtered.noSoft.fixCigar_temp.sam
echo Done

# convert to bam
samtools view -@ 3 -b -o ${libName}_filtered.noSoft.fixCigar_temp.bam ${libName}_filtered.noSoft.fixCigar_temp.sam
rm ${libName}_filtered.noSoft.fixCigar_temp.sam

# convert to modified bed format
bedtools bamtobed -cigar -i ${libName}_filtered.noSoft.fixCigar_temp.bam > ${libName}_reads.bed
rm ${libName}_filtered.noSoft.fixCigar_temp.bam

echo bed lines:
wc -l ${libName}_reads.bed

# run modifyBed.R script to edit bed file so MM analysis can be run on it *Need more than 25G memory for this step
Rscript /n/groups/churchman/mc348/TimelapseSeq/Scripts/modifyBed.R ${libName}_reads.bed ${libName}_filtered.noSoft.fixCigar.noHead_temp.sam
rm ${libName}_reads.bed

mkdir ${libName}_Top${geneNum}_withDups

# split up file to run MM analysis
split --lines=10000 ${libName}_reads.bedwithC.bed ${libName}_Top${geneNum}_withDups/tmp.



########## For quick batch submission ##########
# Libs='G1 G2 G3 G4 G5'
# Libs='H1 H2 H3 H4 H5'
# MapMethod='MMinformed4'
# Expt='mouse_H'
# geneNum='100'
# data='mouse'

# for lib in $Libs
# do
# sbatch -e logs/ProcAlmts_${lib}_${MapMethod}.err -o logs/ProcAlmts_${lib}_${MapMethod}.log ../Scripts/ProcessAlignments_TimelapseSeq_withDups_ChooseTopGenes.sh ${lib}_${MapMethod} $Expt $geneNum $data
# done



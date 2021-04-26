#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-06:00
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

libName=$1
MapMethod=$2
MTstatus=$3
# numReads=$3
data=$4

# Set path to fasta
if [ "${data}" = "Hela" ]
then
fastapath='/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge.fasta'
elif [ "${data}" = "mouse" ]
then
fastapath='/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/mouseNIH3T3_mm10_dm6_ercc_cat.fasta'
fi


module load picard/2.8.0



# Downsample bam so this doesn't take so long
# samtools view -@ 3 -bs 42.3 ${libName}_${MapMethod}_Aligned.sortedByCoord.noSpike.bam > ${libName}_${MapMethod}_Aligned.sortedByCoord.subsamp.bam


# sort by read name (required for fixmate)
samtools sort -@ 3 -n -o ${libName}_${MapMethod}.nameSort.bam ${libName}_${MapMethod}_${MTstatus}.bam #Aligned.sortedByCoord.subsamp

# Remove reads whose pairs did not make it through the bed filter (this doesn't get rid of the ones that have multiple alignments!)
samtools view -h -@ 3 ${libName}_${MapMethod}.nameSort.bam | rev | uniq -D -f16 | rev | samtools view -h -b -S > ${libName}_${MapMethod}.nameSort.paired.bam
echo before first unpaired removal:
samtools view -@ 3 -c ${libName}_${MapMethod}.nameSort.bam
echo after first unparied removal:
samtools view -@ 3 -c ${libName}_${MapMethod}.nameSort.paired.bam

# add fixmate m field (required for markdup)
/n/groups/churchman/bms36/programs/samtools-1.10/bin/samtools fixmate -@ 3 -m ${libName}_${MapMethod}.nameSort.paired.bam ${libName}_${MapMethod}.nameSort.fm.bam
rm ${libName}_${MapMethod}.nameSort.bam
rm ${libName}_${MapMethod}.nameSort.paired.bam

# sort by coordinate (required for markdup)
samtools sort -@ 3 -o ${libName}_${MapMethod}_2temp.bam ${libName}_${MapMethod}.nameSort.fm.bam

echo Sorted input
samtools view -@ 3 -c ${libName}_${MapMethod}_2temp.bam 

# identify PCR and optical duplicates
# -d The optical duplicate distance. Suggested settings of 100 for HiSeq style platforms or about 2500 for NovaSeq ones. Default is 0 to not look for optical duplicates. When set, duplicate reads are tagged with dt:Z:SQ for optical duplicates and dt:Z:LB otherwise. Calculation of distance depends on coordinate data embedded in the read names produced by the Illumina sequencing machines.
# use 100 for Next-seq
# echo Identifying duplicates
# /n/groups/churchman/bms36/programs/samtools-1.10/bin/samtools markdup -s -d 100 -@ 3 ${libName}_${MapMethod}_1temp.bam ${libName}_${MapMethod}_2temp.bam
# echo Done
# rm ${libName}_${MapMethod}_1temp.bam

# report read statistics
echo Running samtools flagstat
samtools flagstat -@ 3 ${libName}_${MapMethod}_2temp.bam
echo Done

# remove: PCR and optical duplicates, secondary alignments, and supplementary alignments
# not primary
samtools view -@ 3 -F 0x100 -o ${libName}_${MapMethod}_4temp.bam ${libName}_${MapMethod}_2temp.bam
rm ${libName}_${MapMethod}_2temp.bam

echo nonprimary removed
samtools view -@ 3 -c ${libName}_${MapMethod}_4temp.bam

# PCR or optical duplicate
# samtools view -@ 3 -F 0x400 -o ${libName}_${MapMethod}_4temp.bam ${libName}_${MapMethod}_3temp.bam
# rm ${libName}_${MapMethod}_3temp.bam

# echo PCR/optical duplicates removed
# samtools view -@ 3 -c ${libName}_${MapMethod}_4temp.bam

# supplementary alignment
samtools view -@ 3 -F 0x800 -o ${libName}_${MapMethod}_5temp.bam ${libName}_${MapMethod}_4temp.bam
rm ${libName}_${MapMethod}_4temp.bam

echo supplementary alignments removed
samtools view -@ 3 -c ${libName}_${MapMethod}_5temp.bam

# reads with unmapped mate
samtools view -@ 3 -F 0x8 -o ${libName}_${MapMethod}_6temp.bam ${libName}_${MapMethod}_5temp.bam
rm ${libName}_${MapMethod}_5temp.bam

echo reads with unmapped mates removed by samtools
samtools view -@ 3 -c ${libName}_${MapMethod}_6temp.bam

# remove reads that map to more than 4 positions, even if it is the primary alignment
samtools view -@ 3 -q 1 -o ${libName}_${MapMethod}_filtered_temp.bam ${libName}_${MapMethod}_6temp.bam
rm ${libName}_${MapMethod}_6temp.bam

echo reads with more than 4 mapping locations removed
samtools view -@ 3 -c ${libName}_${MapMethod}_filtered_temp.bam 

# # Remove unpaired reads again and should get rid of the rest because they should now be unique (note: out of 10 files with 2-5 mil reads, 2 reads got through. Manually remove them)
# # sort by read name (required for uniq)
samtools sort -@ 3 -n -o ${libName}_${MapMethod}_filtered_temp.nameSort.bam ${libName}_${MapMethod}_filtered_temp.bam  
# remove unique qnames
samtools view -h -@ 3 ${libName}_${MapMethod}_filtered_temp.nameSort.bam | rev | uniq -D -f19 | rev | samtools view -h -b -S > ${libName}_${MapMethod}_filtered_temp.paired.nameSort.bam
echo before second unpaired removal:
samtools view -@ 3 -c ${libName}_${MapMethod}_filtered_temp.nameSort.bam
echo after second unpaired removal:
samtools view -@ 3 -c ${libName}_${MapMethod}_filtered_temp.paired.nameSort.bam

# for removing specific reads too (check logs to make sure no reads escaped, then if so do this)
# samtools view -h -@ 3 ${libName}_filtered_temp.nameSort.bam | rev | uniq -D -f19 | rev | grep -v 'A00794:219:HNH7MDRXX:2:2169:10628:14043\|A00794:219:HNH7MDRXX:2:2169:10628:23281' | samtools view -h -b -S > ${libName}_filtered_temp.paired.nameSort.bam

# sort by coordinate (required for biostar84452.jar)
samtools sort -@ 3 -o ${libName}_${MapMethod}_filtered_temp.paired.bam ${libName}_${MapMethod}_filtered_temp.paired.nameSort.bam

######## Downsample to numReads reads ######
# Get readcount in bam
# let "rc=`samtools view -c ${libName}_${MapMethod}_filtered_temp.bam`/2"
# let "ratio=100*$numReads/$rc"
# # Downsample resulting bam so subsequent steps don't take so long
# samtools view -@ 3 -bs 42.${ratio} ${libName}_${MapMethod}_filtered_temp.bam > ${libName}_${MapMethod}_filtered_temp_ds.bam

# convert from .bam to .sam
samtools view -@ 3 -h -O SAM -o ${libName}_${MapMethod}_filtered_temp.sam ${libName}_${MapMethod}_filtered_temp.paired.bam


# remove soft clipped bases from reads
echo removing soft-clipped bases
java -jar /n/groups/churchman/bms36/programs/jvarkit/dist/biostar84452.jar -o ${libName}_${MapMethod}_filtered.noSoft_temp.sam ${libName}_${MapMethod}_filtered_temp.sam
echo done removing soft-clipped bases

rm ${libName}_${MapMethod}_filtered_temp.sam

############ Before running samfixcigar, need to create fasta index and sequence dictionary
############ Load picard: module load picard/2.8.0
############ First, fix fasta so that all lines are the same length (need extra memory for this, I used 50G): java -jar $PICARD/picard-2.8.0.jar NormalizeFasta I=GRCh38_ncRNAs_ERCC_merge.fasta O=GRCh38_ncRNAs_ERCC_merge_linelengthfixed.fasta LINE_LENGTH=60
############ Next, make fasta index: samtools faidx ensGRCh38_h_MTsnp_ncRNAs_allERCC_merge.formatted.fasta
############ Make picard sequence dictionary: java -jar $PICARD/picard-2.8.0.jar CreateSequenceDictionary R=ensGRCh38_h_MTsnp_ncRNAs_allERCC_merge.formatted.fasta O=ensGRCh38_h_MTsnp_ncRNAs_allERCC_merge.formatted.dict
# converting cigar string to discriminate between mismatches and matches
java -jar /n/groups/churchman/bms36/programs/jvarkit/dist/samfixcigar.jar -r $fastapath ${libName}_${MapMethod}_filtered.noSoft_temp.sam > ${libName}_${MapMethod}_filtered.noSoft.fixCigar_temp.sam
rm ${libName}_${MapMethod}_filtered.noSoft_temp.sam

# get rid of sam header
samtools view -@ 3 -O SAM -o ${libName}_${MapMethod}_filtered.noSoft.fixCigar.noHead_temp.sam ${libName}_${MapMethod}_filtered.noSoft.fixCigar_temp.sam

# convert to bam
samtools view -@ 3 -b -o ${libName}_${MapMethod}_filtered.noSoft.fixCigar_temp.bam ${libName}_${MapMethod}_filtered.noSoft.fixCigar_temp.sam
rm ${libName}_${MapMethod}_filtered.noSoft.fixCigar_temp.sam

# convert to modified bed format
bedtools bamtobed -cigar -i ${libName}_${MapMethod}_filtered.noSoft.fixCigar_temp.bam > ${libName}_${MapMethod}_reads.bed
rm ${libName}_${MapMethod}_filtered.noSoft.fixCigar_temp.bam

echo bed lines
wc -l ${libName}_${MapMethod}_reads.bed

# run modifyBed.R script to edit bed file so MM analysis can be run on it *Need more than 25G memory for this step
~/R-3.5.1/library/modifyBed.R ${libName}_${MapMethod}_reads.bed ${libName}_${MapMethod}_filtered.noSoft.fixCigar.noHead_temp.sam
rm ${libName}_${MapMethod}_reads.bed

mkdir ${libName}_${MTstatus}_${MapMethod}_withDups

# split up file to run MM analysis
split --lines=10000 ${libName}_${MapMethod}_reads.bedwithC.bed ${libName}_${MTstatus}_${MapMethod}_withDups/tmp.



########## For quick batch submission ##########
# Libs='TL_0m TL_15m_10 TL_15m_25 TL_15m_50 TL_240m_10 TL_240m_25 TL_240m_50 TL_0m_noC TL_240m_10_noC TL_240m_25_noC TL_240m_50_noC'
# Libs="TL3_0m TL3_7m TL3_15m TL3_30m TL3_45m TL3_60m TL3_90m TL3_120m TL3_240m TL3_0m_noC TL3_120m_noC TL3_240m_noC"
# Libs="TL4_0m TL4_7m TL4_15m TL4_30m TL4_45m TL4_60m TL4_90m TL4_120m TL4_240m"
# Libs='MUT-TL3-0m_TC000349 MUT-TL3-0m_TC000502 MUT-TL3-0m_TC000605 MUT-TL3-0m_TC000745 MUT-TL3-0m_TC00091'
# Libs='MUT-TL3_0m_TC0.000 MUT-TL3_0m_TC0.00349 MUT-TL3_0m_TC0.00502 MUT-TL3_0m_TC0.00605 MUT-TL3_0m_TC0.00745 MUT-TL3_0m_TC0.0091'
# Exp="TL3" #TL5 TL4
# # Libs="0m_tot 15m_tot 30m_tot 60m_tot 0m_poly 15m_poly 30m_poly 60m_poly"
# Libs="0m 7m 15m 30m 45m 60m 90m 120m 240m"
# MapMethod='t5MTMMinformed6'
# MTstatus='MT' # MT MTnorRNA
# data="Hela"
# 
# for lib in $Libs
# do
# sbatch -e logs/ProcAlmts_manRemUnp_${lib}.err -o logs/ProcAlmts_manRemUnp_${lib}.log /n/groups/churchman/mc348/TimelapseSeq/Scripts/ProcessAlignments_TimelapseSeq_MT_withDups.sh ${Exp}_${lib} $MapMethod $MTstatus $data
# done


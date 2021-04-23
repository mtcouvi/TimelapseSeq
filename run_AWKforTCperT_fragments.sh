#!/bin/bash

#SBATCH -c 1
#SBATCH -t 70
#SBATCH -p short
#SBATCH --mem=5G
#SBATCH -o logs/run_AWKforTCperT_%j.out
#SBATCH -e logs/run_AWKforTCperT_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<email>


#### Paired-end SAM strand specification: https://biobeat.wordpress.com/2013/04/29/directional-rna-seq-part-1-extract-strand-information-from-sam-file/
#### Top strand: 99 = R1+,  147 = R2- (primary), 355 = R1+, 403 = R2- (not primary)
#### Bottom strand: 83 = R1-, 163 = R2+, (primary), 339 = R1-, 419 = R2+ (not primary)


### USE: 		Run from one directory up from directories with 
###				reads_MM.fragments.sort.bed file 

######################## For sbatch submission #################

# Exp='TL3'
# Libs="0m 7m 15m 30m 45m 60m 90m 120m 240m"
# MapMethod='MT_t5MTMMinformed6_withDups'
# for lib in $Libs
# do
# sbatch /n/groups/churchman/mc348/TimelapseSeq/Scripts/run_AWKforTCperT_fragments.sh ${Exp}_${lib} $MapMethod
# done
#################################################################



LibName=$1
MapMethod=$2


# All reads
# awk -F'\t' 'BEGIN{OFS="\t"; print "chr\tstart\tend\tstrand\tseq\tfrag_length\ttot_mismatches\tTC_mismatches";} {print $1"\t"$2"\t"$3"\t"$6"\t"$9"\t"$21"\t"$18"\t"$19}' ${LibName}_${MapMethod}/reads_MM.fragments.sort.bed > ${LibName}_${MapMethod}/frag_MMfrequency_All.txt

# And for all mito reads
awk -F'\t' 'BEGIN{OFS="\t"; print "chr\tstart\tend\tstrand\tseq\tfrag_length\ttot_mismatches\tTC_mismatches";} {if ($1 == "h_MT") print $1"\t"$2"\t"$3"\t"$6"\t"$9"\t"$21"\t"$18"\t"$19}' ${LibName}_${MapMethod}/reads_MM.fragments.sort.bed > ${LibName}_${MapMethod}/frag_MMfrequency_MTall.txt


# Mito rRNAs only
awk -F'\t' 'BEGIN{OFS="\t"; print "chr\tstart\tend\tstrand\tseq\tfrag_length\ttot_mismatches\tTC_mismatches";} {if ($1 == "h_MT" && $3 < 3305) print $1"\t"$2"\t"$3"\t"$6"\t"$9"\t"$21"\t"$18"\t"$19}' ${LibName}_${MapMethod}/reads_MM.fragments.sort.bed > ${LibName}_${MapMethod}/frag_MMfrequency_MTrRNA.txt
 
# Mito reads, skip rRNAs
awk -F'\t' 'BEGIN{OFS="\t"; print "chr\tstart\tend\tstrand\tseq\tfrag_length\ttot_mismatches\tTC_mismatches";} {if ($1 == "h_MT" && $3 > 3305) print $1"\t"$2"\t"$3"\t"$6"\t"$9"\t"$21"\t"$18"\t"$19}' ${LibName}_${MapMethod}/reads_MM.fragments.sort.bed > ${LibName}_${MapMethod}/frag_MMfrequency_MTnorRNA.txt





#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-02:00
#SBATCH -p short
#SBATCH --mem=5G
#SBATCH --mail-type=END
#SBATCH --mail-user=mtcouvi@gmail.com
###


LibName=$1
MapMethod=$2

echo $LibName


# Index
# samtools index -@ 3 ${LibName}_${MapMethod}_Aligned.sortedByCoord.out.bam

# Extract mito reads 
samtools view -@ 3 -b ${LibName}_${MapMethod}_Aligned.sortedByCoord.noSpike.bam h_MT:3305 -o ${LibName}_${MapMethod}_MTnorRNA.bam
samtools view -@ 3 -b ${LibName}_${MapMethod}_Aligned.sortedByCoord.noSpike.bam h_MT -o ${LibName}_${MapMethod}_MT.bam


# Index
samtools index -@ 3 ${LibName}_${MapMethod}_MTnorRNA.bam
samtools index -@ 3 ${LibName}_${MapMethod}_MT.bam

# Optional for viewing bams:
# Downsample
# samtools view -@ 3 -bs 42.3 ${LibName}_${MapMethod}_MT.bam > ${LibName}_${MapMethod}_MT.ds.bam
# Split into + and - files
# samtools view -@ 3 -h ${LibName}_${MapMethod}_MT.ds.bam | awk '{if ($1 ~ /^@/ || $2 == "99" || $2 == "355" || $2 == "147" || $2 == "403") print $0}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_MT.ds.P.bam
# samtools view -@ 3 -h ${LibName}_${MapMethod}_MT.ds.bam | awk '{if ($1 ~ /^@/ || $2 == "83" || $2 == "339" || $2 == "163" || $2 == "419") print $0}' | samtools view -@ 3 -h -b -o ${LibName}_${MapMethod}_MT.ds.M.bam
# Index
# samtools index ${LibName}_${MapMethod}_MT.ds.P.bam
# samtools index ${LibName}_${MapMethod}_MT.ds.M.bam

########## Make bamlist AFTER running above on all samples ##########

MapMethod="t5MTMMinformed6"
Exp="TL4"
MTstatus="MT" # MT MTnorRNA

bglib="0m"
no4sUlib="no4sU"
# # # # Copy 0m sample to make 4sU sample for comparison in GS
cp ${Exp}_${bglib}_${MapMethod}_${MTstatus}.bam ${Exp}_${no4sUlib}_${MapMethod}_${MTstatus}.bam
cp ${Exp}_${bglib}_${MapMethod}_${MTstatus}.bam.bai ${Exp}_${no4sUlib}_${MapMethod}_${MTstatus}.bam.bai



# # # # Make bamlist
Exp="TL4" # TL5_poly
touch ${Exp}_${MTstatus}_${MapMethod}.bamlist
# Libs='TL_no4sU TL_0m TL_15m_10 TL_15m_25 TL_15m_50 TL_240m_10 TL_240m_25 TL_240m_50 TL_0m_noC TL_240m_10_noC TL_240m_25_noC TL_240m_50_noC'
# Libs="TL5_no4sU_tot TL5_0m_tot TL5_15m_tot TL5_30m_tot TL5_60m_tot"
# Libs="TL5_no4sU_poly TL5_0m_poly TL5_15m_poly TL5_30m_poly TL5_60m_poly"
Libs="no4sU 0m 7m 15m 30m 45m 60m 90m 120m 240m"
for lib in $Libs
do
echo "${Exp}_${lib}_${MapMethod}_${MTstatus}.bam" >> ${Exp}_${MTstatus}_${MapMethod}.bamlist
done 
# 
# touch Timelapse1_MTall_${MapMethod}.bamlist
# Libs='Timelapse_no4sU Timelapse_0m_CAL Timelapse_15m_CAL Timelapse_30m_CAL Timelapse_60m_CAL Timelapse_120m_CAL Timelapse_240m_CAL Timelapse_noChem240m_CAL'
# for lib in $Libs
# do
# echo "${lib}_${MapMethod}_MTall.bam" >> Timelapse1_MTall_${MapMethod}.bamlist
# done 

# Experiment="K1"
# MapMethod="MMinformed4"
# # 
# # # # Copy 0m sample to make 4sU sample for comparison in GS
# # cp ${Experiment}_0m_${MapMethod}_MT.bam ${Experiment}_no4sU_${MapMethod}_MT.bam
# # cp ${Experiment}_0m_${MapMethod}_MT.bam.bai ${Experiment}_no4sU_${MapMethod}_MT.bam.bai
# cp ${Experiment}_${MapMethod}_MT.bam ${Experiment}_no4sU_${MapMethod}_MT.bam
# cp ${Experiment}_${MapMethod}_MT.bam.bai ${Experiment}_no4sU_${MapMethod}_MT.bam.bai
# cp ${Experiment}_${MapMethod}_MTnorRNA.bam ${Experiment}_no4sU_${MapMethod}_MTnorRNA.bam
# cp ${Experiment}_${MapMethod}_MTnorRNA.bam.bai ${Experiment}_no4sU_${MapMethod}_MTnorRNA.bam.bai
# 
# # # # Make bamlist
# Experiment="K562_J"
# MapMethod="MMinformed4"
# 
# touch ${Experiment}_MT_${MapMethod}.bamlist
# Libs='J1_no4sU J1 J2 J3 J4 J5'
# # Libs='K1_no4sU K1 K2 K3 K4 K5'
# for lib in $Libs
# do
# echo "${lib}_${MapMethod}_MT.bam" >> ${Experiment}_MT_${MapMethod}.bamlist
# done 

# Have to change names of the files so they don't have . in the number
# Libs='MUT-TL3_0m'
# MapMethod='t5MTMMinformed6_MT'
# rates=("0.000" "0.00349" "0.00502" "0.00605" "0.00745" "0.0091")
# newrates=("0000" "000349" "000502" "000605" "000745" "00091")
# fracNew=0.1
# for LibName in $Libs
# do
# for i in {0..5}
# do
# mv ${LibName}_TC${rates[i]}-${fracNew}_${MapMethod}_MT.bam ${LibName}_TC${newrates[i]}-${fracNew}_${MapMethod}_MT.bam
# mv ${LibName}_TC${rates[i]}-${fracNew}_${MapMethod}_MT.bam.bai ${LibName}_TC${newrates[i]}-${fracNew}_${MapMethod}_MT.bam.bai
# done
# done

# # 
# # # # Copy 0m sample to make 4sU sample for comparison in GS
# cp MUT-TL3_0m_TC0000-${fracNew}_${MapMethod}_MT.bam MUT-TL3_0m_TCno4sU-${fracNew}_${MapMethod}_MT.bam
# cp MUT-TL3_0m_TC0000-${fracNew}_${MapMethod}_MT.bam.bai MUT-TL3_0m_TCno4sU-${fracNew}_${MapMethod}_MT.bam.bai
# 
# # # Make bamlist
# touch MUT-TL3_0m-${fracNew}_${MapMethod}.bamlist
# Libs='MUT-TL3_0m_TCno4sU MUT-TL3_0m_TC0000 MUT-TL3_0m_TC000349 MUT-TL3_0m_TC000502 MUT-TL3_0m_TC000605 MUT-TL3_0m_TC000745 MUT-TL3_0m_TC00091'
# for lib in $Libs
# do
# echo "${lib}-${fracNew}_${MapMethod}_MT.bam" >> MUT-TL3_0m-${fracNew}_${MapMethod}.bamlist
# done 
#  
# # 
# # ########## For quick batch submission ##########
# # Libs='TL_0m TL_15m_10 TL_15m_25 TL_15m_50 TL_240m_10 TL_240m_25 TL_240m_50 TL_0m_noC TL_240m_10_noC TL_240m_25_noC TL_240m_50_noC'
# Libs='TL3_0m TL3_7m TL3_15m TL3_30m TL3_45m TL3_60m TL3_90m TL3_120m TL3_240m TL3_0m_noC TL3_120m_noC TL3_240m_noC'
# Libs='TL4_0m TL4_7m TL4_15m TL4_30m TL4_45m TL4_60m TL4_90m TL4_120m TL4_240m'
# Libs='MUT-TL3_0m_TC0.000 MUT-TL3_0m_TC0.00349 MUT-TL3_0m_TC0.00502 MUT-TL3_0m_TC0.00605 MUT-TL3_0m_TC0.00745 MUT-TL3_0m_TC0.0091'
# Libs='SimonK562_none_1 SimonK562_none_2 SimonK562_s6G_1 SimonK562_s6G_2 SimonK562_s4U_1 SimonK562_s4U_2'
# Libs='J1 J2 J3 J4 J5 K1 K2 K3 K4 K5'
# Libs='0m_tot 15m_tot 30m_tot 60m_tot 0m_poly 15m_poly 30m_poly 60m_poly'
Libs="0m 7m 15m 30m 45m 60m 90m 120m 240m"
MapMethod='t5MTMMinformed6'
Expt='TL4'
for lib in $Libs
do
sbatch -e logs/ExtractMito_${lib}.err -o logs/ExtractMito_${lib}.log ../Scripts/ExtractMitoReadsForGRANDSLAM.sh ${Expt}_${lib} $MapMethod
done


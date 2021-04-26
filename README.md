# TimelapseSeq Mismatch calling

### Preprocessing steps:  
    a. Trim adaptors: TrimAdapter_Timelapse.sh  
    b. Align: MapForGRANDSLAM.sh  
    c. Filter reads:  
        opt1. Mito-mapping reads: ExtractMitoReadsForGRANDSLAM.sh  
        opt2. Top genes:  
            c2a. Combine readCounts files (from STAR)
            c2b. AddGeneName_RPKMS_customGTF.R
            c2c. ChooseTopGenesSummingNreadsANDSelectGenesInBed.R  
                If using a new genome, will need to make the bed file for this script with GENCODE_gtf2bed.py
         


### Running Mismatch scripts
1. ProcessAlignments_TimelapseSeq_MT_withDups.sh libName MapMethod MTstatus data   
or  
ProcessAlignments_TimelapseSeq_ChooseTopGenes_withDups.sh libName MapMethod geneNum data

    ** Be sure fasta path in findMismatches_complete_20_02_28_MC.R matches reference used for mapping: this info is set on command line with <data> ***

2. Call mismatches

    From directory one up from tmp files
    
    example)
    ```
    Exp="TL4"
    Libs="0m 7m 15m 30m 45m 60m 90m"
    dirName='MT_t5MTMMinformed6_withDups'
    genome="Hela"
    
    for lib in $Libs
    do  
    sbatch -e logs/submit_parallel_${lib}_${dirName}.err -o logs/submit_parallel_${lib}_${dirName}.log /n/groups/churchman/mc348/TimelapseSeq/Scripts/submit_parallel_findMM_MC_batchSubmit.sh ${Exp}_${lib}_${dirName} $genome
    done  
    ```
    ```
    ls -l slurm* #to check size of all slurms files
    rm slurm*
    ```
3. Check to make sure all files were created:
From directory one up from tmp files
 
    ```
    for lib in $Libs
    do
    /n/groups/churchman/mc348/TimelapseSeq/Scripts/test_submitMM_batchSubmit.sh ${Exp}_${lib}_${dirName}
    done
    ```
    
4. Concatenate files
    ```
    for lib in $Libs
    do
    sbatch -p short -t 0-01:00 --wrap="cat ${Exp}_${lib}_${dirName}/tmp*MM.bed > ${Exp}_${lib}_${dirName}/reads.MM.bed" 
    done
    ```
    ```
    rm slurm*
    rm *_${dirName}/tmp*
    ```
6. Sort
    ```
    for lib in $Libs
    do
    sbatch /n/groups/churchman/mc348/TimelapseSeq/Scripts/run_sortReads_byName_batchSubmit.sh ${Exp}_${lib}_${dirName}
    done
    ```
    ```
    # check slurms
    rm slurm*
    # remove unsorted files
    rm *_${dirName}/reads.MM.bed
    ```
    
7. Combine files into fragments
    ```
    for lib in $Libs
    do
    sbatch -e logs/run_MM_per_frag_${Exp}_${lib}_${dirName}.err -o logs/run_MM_per_frag_${Exp}_${lib}_${dirName}.log /n/groups/churchman/mc348/TimelapseSeq/Scripts/run_MM_per_fragment_20_02_28_MC_batchSubmit.sh ${Exp}_${lib}_${dirName}
    done
    ```
    ```
    check slurms 
    rm slurm* 
    ```
    
8. Check that all files were created
    ```
    for lib in $Libs
    do
    /n/groups/churchman/mc348/TimelapseSeq/Scripts/test_MM_per_F.sh ${Exp}_${lib}_${dirName}
    done
    ```

9. Concatenate files
    ```
    for lib in $Libs
    do
    sbatch -p short -t 0-01:00 --wrap="cat ${Exp}_${lib}_${dirName}/tmp*fragments.bed > ${Exp}_${lib}_${dirName}/reads_MM.fragments.bed" 
    done
    ```
    ```
    #check slurms
    rm slurm*
    rm *_${dirName}/tmp*
    ```
    
10. Sort
    ```
    for lib in $Libs
    do
    sbatch -p short -t 0-03:00 --mem=50G --wrap="bedtools sort -i ${Exp}_${lib}_${dirName}/reads_MM.fragments.bed > ${Exp}_${lib}_${dirName}/reads_MM.fragments.sort.bed"
    done
    ```

11. Remove unsorted file
   
    ```
    rm *_${dirName}/reads_MM.fragments.bed
    ```

12. Shorten reads_MM.fragments.sort.bed (remove columns and select only reads of interest)
    need to modify this script for reads of interest
    ```
    Exp='TL3'
    Libs="0m 7m 15m 30m 45m 60m 90m 120m 240m"
    MapMethod='MT_t5MTMMinformed6_withDups'
    
    for lib in $Libs
    do
    sbatch /n/groups/churchman/mc348/TimelapseSeq/Scripts/run_AWKforTCperT_fragments.sh ${Exp}_${lib} $MapMethod
    done
    ```
13. R scripts to plot MM frequencies
    MismatchFrequencyAveTsAndNumReadsWithMM.R




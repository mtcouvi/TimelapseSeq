
library(data.table)
library(rlist)

pseudo = 'exclude' # include exclude
ncRNA = 'exclude'
n = 100

Experiment = 'Brendan_mouseNIH3T3' #'Hela_2020_09' 'Brendan_K562_InSilico'
Expt = 'mouse_H'
# fileName1 = paste0(Expt,'_MMinformed4_ReadsPerGene') #'TL3_MMinformed4_ReadsPerGene' J1K1_MMinformed5snpCor_ReadsPerGene
path = paste0('/Users/Mary/Desktop/Data/TimelapseSeq/',Experiment,'/RPK/')
inputBedFile = 'mouseNIH3T3_mm10_MTmod.bed' # hg38_Actually_hg37_J1K1snpIndelMasked.bed


if (ncRNA == 'include') {
DT <- data.table(read.table(paste0(path, fileName, '_reads_noPseudo.txt'), sep = '\t', quote = '', header=TRUE, stringsAsFactors = FALSE))
}
if (ncRNA == 'exclude') {
DT <- data.table(read.table(paste0(path, fileName, '_reads_noPseudo_ncRNA.txt'), sep = '\t', quote = '', header=TRUE, stringsAsFactors = FALSE))
}

numSamps = ncol(DT)-3
samples = colnames(DT)[2:(numSamps+1)]
for (i in c(1:numSamps)) {
as.numeric(DT[[samples[i]]])}


# Remove mito genes
if (Experiment == 'Brendan_K562_InSilico') {DT <- DT[!grep('^MT-', GeneName)]}
if (Experiment == 'Brendan_mouseNIH3T3') {DT <- DT[!grep('^mt-', GeneName)]}


# Sort in descending order
DTsort1 <- DT[order(-DT[[samples[1]]])]
DTsort5 <- DT[order(-DT[[samples[5]]])]

TotReadCount1 <- sum(head(DTsort, n)[[samples[1]]]) # This is now 1594604 without mito genes or ncRNAs
TotReadCount5 <- sum(head(DTsort, n)[[samples[5]]]) # This is now 1594604 without mito genes or ncRNAs


# Make a bed file from the list with coordinates, to use to pull the genes out of a bam file
# Get the genes
TopGenes <- head(DTsort1, n)$GeneName

# Open bed file with gene coordinates
# First have to make a bed from the modified gtf (use GENCODE_gtf2bed.py)
bedDT <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/TimelapseSeq/SeqFiles/', inputBedFile), col.names=c('chr','start','end','gene_name','null','strand')))

# Pull top genes out of bed
TopBED <- bedDT[bedDT$gene_name %in% TopGenes]

# Get coordinates for longest transcript, +/- 1 kb
TopBED[, `:=` (newStart = min(start)-1000, newEnd = max(end)+1000),  by = gene_name]
# Just keep one entry for each transcript (all have same newStart and newEnd)
TopBED <- unique(TopBED, by='gene_name')
# Rewrite table with newStart and newEnd
newTopBED <- data.table(TopBED$chr, TopBED$newStart, TopBED$newEnd, TopBED$gene_name, TopBED$null, TopBED$strand)

# Write table
write.table(newTopBED, file=paste0(path,'/Top',n,'_', Expt,'.bed'), row.names=FALSE, col.names=FALSE, sep=("\t"), quote=FALSE)

# source('/Users/Mary/Desktop/Data/TimelapseSeq/Scripts/ChooseTopGenesSummingNreadsANDSelectGenesInBed.R')

########### NOTES ##############
################################

# Excluding RNA genes, top 300 (1.4%) of genes give about 3 million reads (8 are mito mRNAs):
# sum(head(DTsortJ1, 300)$J1) 
# 3081907
# sum(head(DTsortJ1, 300)$K1)
# 3051485
# sum(head(DTsortK1, 300)$K1) 
# 3183001
# sum(head(DTsortK1, 300)$J1)
# 2956594

# Including RNA genes, top 150 (0.6%) of genes give about 3 million reads (9 are mito genes):
# sum(head(DTsortJ1, 150)$J1) 
# 3148269
# sum(head(DTsortJ1, 150)$K1)
# 3057143
# sum(head(DTsortK1, 150)$K1) 
# 3146893
# sum(head(DTsortK1, 150)$J1)
# 3043368
 
# >>> Talking to Erik decided to go with coding and nuclear, and to be more stringent: top 100 genes

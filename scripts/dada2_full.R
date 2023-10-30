# Load libraries
library("BiocManager")
library("dada2")
library("phyloseq")
library("BiocGenerics")
library("S4Vectors")
library("IRanges")
library("XVector")
library("Biostrings")
library("permute")
library("vegan")
library("GenomeInfoDb")
library("GenomicRanges")
library("matrixStats")
library("MatrixGenerics")
library("Biobase")
library("SummarizedExperiment")
library("DESeq2")
library('devtools')
library("cluster")
library("pairwiseAdonis")
library("biomformat")
library("argparser")

# Currently we only have run 1 (run 2 will come after christmas)
new_day <- gsub("-", "", as.character(Sys.Date()))

# Create parser
# parser <- arg_parser("Script for dada2 processing")
# 
# # Add command line arguments
# parser <- add_argument(parser, "--path", help = "Path to project")
# parser <- add_argument(parser, "--path_fastq", help = "Path to fastq_trimmed")
# parser <- add_argument(parser, "--pathQC", help = "Path to QC directory")
# 
# args <- parse_args(parser)
# 
# # Load files from parser
# path <- args$path
# path_fastq <- args$path_fastq
# pathQC <- args$pathQC


# Load files example:
path <- "~/workspace/MPI/SG_MemoryRapamycin/"
path_fastq <- "~/workspace/MPI/SG_MemoryRapamycin/data/fastq_trimmed/"
pathQC <- "~/workspace/MPI/SG_MemoryRapamycin/analysis/dada_quality/"

# forward and reverse fastq filenames have format: sample-name_xx_xx_R1_001.fastq
fnFs <- sort(list.files(path_fastq, pattern="_R1_001_trimm.fastq", full.names = TRUE))
fnRs <- sort(list.files(path_fastq, pattern="_R2_001_trimm.fastq", full.names = TRUE))
# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
names(fnFs) <- sample.names
names(fnFs) <- sample.names
# Get the filtered ones

filtFs <- file.path(path, "data/fastq_filtered/", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "data/fastq_filtered/", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)



getN <- function(x) sum(getUniques(x))
track <- cbind(out)
#track <- cbind(sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
#colnames(track) <- c("merged", "nonchim")
colnames(track) <- c("input", "filtered")
rownames(track) <- sample.names
write.table(track, paste0(pathQC,"QC_trackReads_input.csv"), sep = ";", quote = F)


errF <- learnErrors(filtFs, multithread = TRUE)
saveRDS(errF, paste0(pathQC, "errF-", new_day, ".rds"))
errR <- learnErrors(filtRs, multithread=TRUE)
saveRDS(errF, paste0(pathQC, "errR-", new_day, ".rds"))


# Check that all files exist and dont evaluate the ones that dont
exists <- file.exists(filtFs)
sample.names <- sample.names[exists]

# Dereplicate FASTQ files to speed up computation
derepFs = derepFastq(filtFs, verbose = TRUE)
derepRs = derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) = sample.names
names(derepRs) = sample.names

# Apply core sequence-variant inference algorithm 
dadaRs = dada(derepRs, err=errR, multithread=TRUE, pool = "pseudo")
dadaFs = dada(derepFs, err=errF, multithread=TRUE, pool = "pseudo")

# merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Make sequence table for analysis
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, paste0(pathQC, "seqtab-", new_day,".rds"))
# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
table(nchar(getSequences(seqtab)))
dim(seqtab.nochim)
saveRDS(seqtab.nochim, paste0(pathQC, "seqtab.nochim-", new_day,".rds"))
sum(seqtab.nochim)/sum(seqtab)

##Quality controls
cat("Dimensions of original sequence table",file=paste0(pathQC, "QC_run_pipeline.txt"),sep="\n")
cat(dim(seqtab),file=paste0(pathQC, "QC_run_pipeline.txt"),append=TRUE)
## Inspect distribution of sequence lengths
cat("\nDistribution of sequence lengths\n",file=paste0(pathQC, "QC_run_pipeline.txt"),append=TRUE)
write.table(table(nchar(getSequences(seqtab))),file=paste0(pathQC, "QC_run_pipeline.txt"),append=TRUE, row.names = FALSE, col.names = FALSE)
## Dimensions after removing chimeras
cat("\nDimensions of sequence table after removing chimeras\n",file=paste0(pathQC, "QC_run_pipeline.txt"),append=TRUE)
cat(dim(seqtab.nochim),file=paste0(pathQC, "QC_run_pipeline.txt"),append=TRUE)
# Percent of chimeras in the original dataset
cat("\nPercentage of chimeras in the original sequence table\n",file=paste0(pathQC, "QC_run_pipeline.txt"),append=TRUE)
cat(sum(seqtab.nochim)/sum(seqtab),file=paste0(pathQC, "QC_run_pipeline.txt"),append=TRUE)
# Write it into a csv
write.table(track, paste0(pathQC,"QC_trackReads.csv"), sep = ";", quote = F)


taxa <- assignTaxonomy(seqtab.nochim, paste0(path, "metadata/databases/silva_nr99_v138.1_wSpecies_train_set.fa.gz"), multithread = TRUE)
saveRDS(taxa, paste0(pathQC, "taxa-", new_day, ".rds"))

# And add species
taxa.plus <- addSpecies(taxa, paste0(path, "metadata/databases/silva_species_assignment_v138.1.fa.gz"))

saveRDS(taxa.plus, paste0(pathQC, "taxa-plus-", new_day, ".rds"))


st.biom <- make_biom(t(seqtab.nochim))
write_biom(st.biom, paste0(pathQC, "ASV_table.biom"))
# Also save what we have until now
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste0(pathQC, "ASVs.fa"))
# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, paste0(pathQC, "ASVs_counts.tsv"), sep="\t", quote=F, col.names=NA)
st.biom = make_biom(asv_tab)
write_biom(st.biom, paste0(pathQC, "ASV_table.biom"))
# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, paste0(pathQC, "ASVs_taxonomy.tsv"), sep="\t", quote=F, col.names=NA)

## Quality control
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
cat("\nFirst rows of assigned taxa\n",file=paste0(pathQC, "QC_run_pipeline.txt"),append=TRUE)
cat(head(taxa.print),file=paste0(pathQC, "QC_run_pipeline.txt"),append=TRUE)
cat("\nFirst rows of assigned species\n",file=paste0(pathQC, "QC_run_pipeline.txt"),append=TRUE)
cat(head(unname(taxa.plus)),file=paste0(pathQC, "QC_run_pipeline.txt"),append=TRUE)


# Filter out samples with total number of reads < 1500
# This also removes samples with less than 10 ASVs identified
clean_seqtab.nochim <- seqtab.nochim[which(rowSums(seqtab.nochim) > 1000), ]
# Filter out OTUs found in less than 5 samples
clean_seqtab.nochim <- clean_seqtab.nochim[, colSums(clean_seqtab.nochim) > 4]
# Filter out OTUs with total sum below 5 (depth of sequencing)
clean_seqtab.nochim <- clean_seqtab.nochim[,which(colSums(clean_seqtab.nochim) > 4)]

# Save clean seqtab files for analysis
saveRDS(clean_seqtab.nochim, paste0(path, "analysis/CLEAN_merged_seqtabNochim_", new_day, ".rds"))

# Asign TAXA again since our dataframe is now definitely clean
taxa <- assignTaxonomy(clean_seqtab.nochim, paste0(path, "metadata/databases/silva_nr99_v138.1_wSpecies_train_set.fa.gz"), multithread = TRUE)
taxa.plus <- addSpecies(taxa, paste0(path, "metadata/databases/silva_species_assignment_v138.1.fa.gz"))
saveRDS(taxa, paste0(path, "analysis/CLEAN_taxonomy_merged_", new_day, ".rds"))
saveRDS(taxa.plus, paste0(path, "analysis/CLEAN_taxa-plus_merged_", new_day, ".rds"))

# Getting intermediate files
st.biom <- make_biom(t(clean_seqtab.nochim))
write_biom(st.biom, paste0(path, "analysis/CLEAN_ASV_table_merged_", new_day, ".biom"))
# Also save what we have until now
asv_seqs <- colnames(clean_seqtab.nochim)
asv_headers <- vector(dim(clean_seqtab.nochim)[2], mode="character")
for (i in 1:dim(clean_seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste0(path, "analysis/CLEAN_ASVs_merged_", new_day, ".fa"))
# count table:
asv_tab <- t(clean_seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, paste0(path, "analysis/CLEAN_ASVs_counts_merged_", new_day, ".tsv"), sep="\t", quote=F, col.names=NA)
st.biom = make_biom(asv_tab)
write_biom(st.biom, paste0(path, "analysis/CLEAN_ASV_table_merged2_", new_day, ".biom"))
# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, paste0(path, "analysis/CLEAN_ASVs_taxonomy_merged_", new_day, ".tsv"), sep="\t", quote=F, col.names=NA)





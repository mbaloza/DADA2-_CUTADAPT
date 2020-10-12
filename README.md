# Getting ready

if (!requireNamespace("BiocManager", quietly=TRUE))
   install.packages("BiocManager")

if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")

BiocManager::install("ShortRead") 

BiocManager::install("additional_missing_package_1")# Or you can use earlier versions of Biocondcutor, e.g. if you have a pre-3.6 version of R

library("devtools")
devtools::install_github("benjjneb/dada2")

library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

path <- "/home/bios1/mbaloza/PS118_Seq/MyData_test/"  ## CHANGE ME to the directory containing the fastq files.
list.files(path)

fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))

fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Identify primers

FWD <- "CCTACGGGNGGCWGCAG"  ## CHANGE ME to your forward primer sequence

REV <- "GACTACHVGGGTATCTAATCC"  ## CHANGE ME...


allOrients <- function(primer) {
  #Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory

fnRs.filtN <- file.path(path, "filtN", basename(fnRs))


fn <- filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)


primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0)) # Counts number of reads in which the primer is found
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# Remove Primers

cutadapt <- "/home/bios1/mbaloza/.local/bin/cutadapt"  # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R


path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)

REV.RC <- dada2:::rc(REV)

#Trim FWD and the reverse-complement of REV off of R1 (forward reads)

R1.flags <- paste("-g", FWD, "-a", REV.RC) 

#Trim REV and the reverse-complement of FWD off of R2 (reverse reads)

R2.flags <- paste("-G", REV, "-A", FWD.RC) 

#Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))


#Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))

cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

#Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# Inspect read quality profiles

plotQualityProfile(cutFs[1:2])

plotQualityProfile(cutRs[1:2])

# Filter and trim

filtFs <- file.path(path.cut, "filtered", basename(cutFs))

filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(240,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # on windows, set multithread = FALSE
head(out)

# Learn the Error Rates

errF <- learnErrors(filtFs, multithread = TRUE)

errR <- learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)

# Sample Inference

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)


# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

#Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct Sequence Table

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

table(nchar(getSequences(seqtab.nochim)))
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


# Track reads through the pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
#If processing a single sample, remove the sapply calls: e.g. replace
#sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names
head(track)

write.table (track, file = "/home/bios1/mbaloza/PS118_Seq/Track_reads.csv", row.names = F, sep = ",")

# Assign taxonomy
unite.ref <- "/home/bios1/mbaloza/PS118_Seq/MyData_test/silva_nr99_v138_train_set.fa.gz"  # CHANGE ME to location on your machine
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

write.table (taxa.print, file = "/home/bios1/mbaloza/PS118_Seq/Taxa.csv", row.names = F, sep = ",")

# Assign Species
taxa_Sp <- addSpecies(taxa, "/home/bios1/mbaloza/PS118_Seq/MyData_test/silva_species_assignment_v138.fa.gz")

taxa_Sp <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)

taxa_Sp.print <- taxa_Sp  # Removing sequence rownames for display only
rownames(taxa_Sp.print) <- NULL
head(taxa_Sp.print)
write.table (taxa_Sp.print, file = "/home/bios1/mbaloza/PS118_Seq/Taxa_SP.csv", row.names = F, sep = ",")

# Evaluate accuracy

rownames(seqtab.nochim) 

unqs.mock <- seqtab.nochim["PS118-2-MUC1-C3-D4",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")




#download the file
download.file(
  url = "https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip",
  destfile = "miseqsopdata.zip")

#unzip the to get individual fastq files

unzip("miseqsopdata.zip")

#record path to data

data_dir <-"MiSeq_SOP"

list.files(data_dir)
read.delim(base::paste(data_dir,"mouse.dpw.metadata", sep = "/" ))[
  order(read.delim(base::paste(data_dir,"mouse.dpw.metadata", sep = "/" ))[[2]]),]
#Taxonomy file
download.file(
  url = "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz?download=1",
  destfile = "silva_nr99_v138.1_wSpecies_train_set.fa.gz")
BiocManager::install("dada2", force=TRUE)

BiocManager::install("phyloseq",force=TRUE)
lapply(c("dada2","phyloseq","vegan","tidyverse",
         "Biostrings"), library, character.only = T)
# All forward and reverse fastq file names have the exact format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
#Note the pattern

fnFs <- sort(list.files(data_dir, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(data_dir, pattern="_R2_001.fastq", full.names = TRUE))
#basename removes the folder names
sample.names <- stringr::str_split_fixed(basename(fnFs),"_",2)[,1]

sample.names
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
plotQualityProfile(fnFs[1:38])
plotQualityProfile(fnRs[1:38])

#R1
filtFs <- file.path("filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names
#R2
filtRs <- file.path("filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ = 2 ,rm.phix=TRUE,
                     compress=TRUE, multithread = TRUE)

#Some reshuffling of the data, get percent of reads left per sample as well as absolute numbers
out2 <- as.data.frame(out) %>% mutate(prop = (100*(reads.out/reads.in)))

# So far looks good
head(out2)
#Mean percent of seqs left
mean(out2$prop)
#Min percent of seqs left
min(out2$prop)
#But what's the sample with the lowest number of reads?
min(out2$reads.out)
errF <- learnErrors(filtFs, multithread=T)
errR <- learnErrors(filtRs, multithread=T)
plotErrors(errF, nominalQ=TRUE)
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
plotErrors(errR, nominalQ=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread = TRUE)
dadaRs <- dada(filtRs, err=errF, multithread = TRUE)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
DT::datatable(mergers$F3D0,options = list(scrollX = TRUE, autoWidth = TRUE))
seqtab <- makeSequenceTable(mergers)

dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=F, verbose=TRUE)
#Proportion ASVs were removed due to chimeras?
100-(dim(seqtab.nochim)[2]/dim(seqtab)[2]*100)
sum(seqtab.nochim)/sum(seqtab)*100
# function that sums the number of unique occurrences in this case seqs
getN <- function(x){ 
  
  sum(getUniques(x))
  
}

#Using apply to run the function on all list elements.
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

#adding colnames
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

#Adding row names
rownames(track) <- sample.names

track <- as.data.frame(track) %>% mutate(prop_kept = round(nonchim/input*100,2))

DT::datatable(track)
ASV_tbl <- as.data.frame(seqtab.nochim) %>% 
  `colnames<-`(base::paste("ASV", seq(1:ncol(seqtab.nochim)), sep = ""))

DT::datatable(ASV_tbl, options = list(scrollX = TRUE, autoWidth = TRUE))
ASV_tbl <- as.data.frame(t(ASV_tbl[-nrow(ASV_tbl),]))

#lets remove the ASVs that were only present in the mock communities

ASV_tbl <- ASV_tbl %>% mutate(rs = rowSums(ASV_tbl)) %>% 
  filter(rs > 0) %>% dplyr::select(-rs)
#Using DNAStringSet from Biostrings to read the seqs 
ASV_seqs <- DNAStringSet(base::colnames(seqtab.nochim) %>% 
                           `names<-`(base::paste("ASV",seq(1:ncol(seqtab.nochim)),sep = "")))
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)
taxa2 <- as.data.frame(taxa) %>% 
  `row.names<-`(base::paste("ASV", seq(1:nrow(taxa)), sep = ""))

DT::datatable(taxa2,options = list(scrollX = TRUE, autoWidth = TRUE))
#Remove non_prokaryotic ASVs
taxa2 <- taxa2 %>% filter(Family != "Mitochondria") %>% 
  filter(Order != "Chloroplast")

#Prepare for joining
taxa3 <- taxa2 %>% mutate(ASVs = row.names(taxa2)) %>% 
  dplyr::select(ASVs, everything())

ASV_tbl2 <- ASV_tbl %>% mutate(ASVs=row.names(ASV_tbl)) %>%
  dplyr::select(ASVs, everything())


#Join taxa table this will be used for analysis outside of phyloseq
taxa3 <- taxa3 %>% inner_join(ASV_tbl2, by = "ASVs")

#taxa2 will be used for Phyloseq
taxa2 <- taxa2 %>% filter(rownames(taxa2)%in%taxa3$ASVs)

#Filter the original ASV table 
ASV_tbl <- ASV_tbl %>% filter(row.names(ASV_tbl) %in% taxa3$ASVs)

# get a table of relative abundance of taxa per sample
#Since the samples are columns gotta use the sweep function rather than just divide

taxa3[,(ncol(taxa2)+2):ncol(taxa3)] <- sweep(
  taxa3[,(ncol(taxa2)+2):ncol(taxa3)],2,
  colSums(taxa3[,(ncol(taxa2)+2):ncol(taxa3)]),'/')*100

#Sanity check, do they sum to 100%
colSums(taxa3[,(ncol(taxa2)+2):ncol(taxa3)])
ASV_seqs[names(ASV_seqs)%in%ASV_tbl2$ASVs]
writeXStringSet(ASV_seqs[names(ASV_seqs)%in%ASV_tbl2$ASVs], "mouse_ASVs.fa")

#lets check that the file is ok

readLines("mouse_ASVs.fa")[1:10]

# Firmicutes
firmicutes_ids <- rownames(taxa2)[taxa2$Phylum == "Firmicutes"]
firmicutes_seqs <- ASV_seqs[names(ASV_seqs) %in% firmicutes_ids]
genus <- as.character(taxa2[firmicutes_ids, "Genus"])
names(firmicutes_seqs) <- paste(firmicutes_ids, genus)
writeXStringSet(firmicutes_seqs, "Firmicutes_ASVs_identified.fa")
head(firmicutes_seqs)

# Actinobacteria
actino_ids <- rownames(taxa2)[taxa2$Phylum %in% c("Actinobacteria", "Actinobacteriota")]
actino_seqs <- ASV_seqs[names(ASV_seqs) %in% actino_ids]
genus_actino <- as.character(taxa2[actino_ids, "Genus"])
names(actino_seqs) <- paste(actino_ids, genus_actino)
writeXStringSet(actino_seqs, "Actinobacteria_ASVs_identified.fa")
actino_seqs

unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
mock.ref <- getSequences(file.path(data_dir, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
#lets see what it looks like

mdta <- read.delim(base::paste(data_dir,"mouse.dpw.metadata", sep = "/" ))
#Are the sample names similar to the sample names we have in the ASV table?
length(mdta$group %in% colnames(ASV_tbl))
#set the row names and massage the table

mdta <- mdta %>% mutate(category = if_else(dpw > 100 , "Late","Early")) %>% 
  mutate(gender = if_else(str_detect(group, "F"),"F","M")) %>%
  mutate(subject = str_extract(str_split_fixed(group, "D", 2)[,1],"(\\d)+")) %>%
  `row.names<-`(mdta$group)
#Generating the phyloseq object
ps <- phyloseq(otu_table(as.matrix(ASV_tbl), taxa_are_rows=T), 
               sample_data(mdta), 
               tax_table(as.matrix(taxa2)))

ps <- merge_phyloseq(ps, ASV_seqs)

ps
#Lets calculate richness and diversity
#Shannon and Simpson are diversity estimators, ACE is a richness estimator, observed is how many ASVs we have (198)
rnd <- estimate_richness(
  ps, measures=c("Observed", "ACE","Shannon", "Simpson")) %>% 
  dplyr::select(-se.ACE)


DT::datatable(rnd[1:2])
#QUESTION4
data.table::fwrite(rnd[, 1:2], "datatable_export.csv")
rnd <- rnd %>% mutate(group = row.names(rnd))%>% inner_join(mdta, "group")
#make a long table
rndl <- rnd %>% gather("Estimator", "Value", 1:4) %>% 
  mutate(rd = if_else(Estimator == "Observed" | Estimator == "ACE", 
                      "Richness","Diversity"))

ggplot(rndl, aes(x=dpw, y=Value, color = Estimator, fill= Estimator))+
  geom_bar(stat = "identity", position = position_dodge())+
  facet_grid(rd~category, scales = "free")
#QUSTION 5

library(ggplot2)
library(RColorBrewer)

cb_palette <- brewer.pal(n = length(unique(rndl$Estimator)), name = "Set2")

ggplot(rndl, aes(x = dpw, y = Value, color = Estimator, fill = Estimator)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_grid(rd ~ category, scales = "free") +
  scale_color_manual(values = cb_palette) +
  scale_fill_manual(values = cb_palette) +
  theme_bw()
# Log transforming the data prior to using bray curtis so that overly abundant ASVs will not skew the results
ps.prop <- transform_sample_counts(ps, function(otu) log1p(otu))
#Perform the ordination
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="category", title="Bray NMDS")
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU)*100)
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="dpw", fill="Family") + facet_wrap(~category, scales="free_x") +
  ylab("Relative abundance (% of seqs/sample)")
fam <- as.data.frame(taxa3 %>% dplyr::select(Family, where(is.numeric)) %>% 
                       group_by(Family) %>% summarise_all(sum))

head(fam)
#Does it sum to 100%
colSums(fam[,2:ncol(fam)])
fam_o <- fam

fam[2:ncol(fam)][fam[2:ncol(fam)]<1] <- 0
fam_o[2:ncol(fam_o)][fam[2:ncol(fam_o)]>1] <- 0

Others <-colSums(fam_o[2:ncol(fam_o)])

fam[nrow(fam)+1,2:ncol(fam)] <- Others
fam[nrow(fam),1] <- "Others"

tail(fam)
#Get a long table needed for ggplot2
fam_l <- fam %>% gather("group","relab",2:ncol(fam)) %>%
  inner_join(mdta,"group")
txplot <- ggplot(fam_l, aes(x = dpw, y=Family, size = relab))+ 
  geom_point()+
  facet_wrap(~category, scales = "free_x")
txplot
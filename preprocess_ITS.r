library("tidyverse")
library("qiime2R")
library ("phyloseq")

### Load data to R
## Begin

dir = as.character(read.table("dir.txt"))
setwd(dir)



asv_meta = read_qza("./ITS/dada2/FeatureTable[Frequency]_ITS.qza") # Table of metagenome data reads for each ASV by samples

asv_meta = cbind(row.names(asv_meta$data), asv_meta$data)
asv_meta = as.data.frame(asv_meta)
colnames(asv_meta)[1] = "Species"
row.names(asv_meta) = asv_meta[,1]
asv_meta = asv_meta[,-1]
asv_meta[] = apply(asv_meta[], 2, as.numeric)



tax_meta = read_qza("./ITS/feature-classifier_classify-sklearn/FeatureData[Taxonomy]_ITS.qza") # Taxonomy table for metagenome data
tax_meta = parse_taxonomy(tax_meta$data, trim_extra=FALSE)
tax_meta[is.na(tax_meta)] <- "uncultured"

# Swap values ​​from Swap_unidentified in tax_meta to previous taxon level
Swap_unidentified = c("uncultured", "unidentified", "metagenome", "bacteriap25", "Unknown")
for(Swap in Swap_unidentified) {
    for(j in 1:7) {
        for (i in 1:nrow(tax_meta)) {
            if(grepl(Swap, tax_meta[i,j])) {
                tax_meta[i,j] = tax_meta[i,j-1] }
                else { tax_meta[i,j] = tax_meta[i,j]
            }       
        }
    }
}


#

# Filtering metagenome data from non-significant taxa
tax_meta = filter(tax_meta, 
Phylum != "d__Bacteria" & 
Phylum != "Unassigned" &
Kingdom != "d__Eukaryota" & 
Phylum != "k__Fungi" & 
Kingdom != "d__Archaea" &
Kingdom != "k__Viridiplantae" &
Genus != "g__Mitochondria" &
Genus != "g__Chloroplast")
#

tax_meta["Other",] = c("Other","Other","Other","Other","Other","Other","Other")

# Name ASVs to Genus level
tax_meta[grep("f__", tax_meta[,"Genus"]),"Genus"] = make.unique(tax_meta[grep("f__", tax_meta[,"Genus"]), "Genus"])
tax_meta[grep("o__", tax_meta[,"Genus"]),"Genus"] = make.unique(tax_meta[grep("o__", tax_meta[,"Genus"]), "Genus"])
tax_meta[grep("c__", tax_meta[,"Genus"]),"Genus"] = make.unique(tax_meta[grep("c__", tax_meta[,"Genus"]), "Genus"])
tax_meta[grep("p__", tax_meta[,"Genus"]),"Genus"] = make.unique(tax_meta[grep("p__", tax_meta[,"Genus"]), "Genus"])
tax_meta[,"Genus"] = str_replace_all(tax_meta[,"Genus"], "(\\d+)$", function(x) as.numeric(x) + 1)
#

sampledata_ITS = read.table(file="./Metagenome metadata/sampledata_ITS_meta.txt", header = TRUE, colClasses = "character", sep="\t") # Sample metadata for the metagenome
row.names(sampledata_ITS) = sampledata_ITS$SRA


# Creating a phyloseq object for metagenome data
sampledata_meta = sample_data(sampledata_ITS)

tax_ITS = tax_table(as.matrix(tax_meta))
asv_ITS = otu_table(asv_meta, taxa_are_rows = TRUE)
physeq_ITS = phyloseq(asv_ITS, tax_ITS, sampledata_meta)


# phyloseq object for sow (ITS)
sampledata_ITS_sow = read.table(file="./ITS/Culture-dependent data/sampledata_ITS_sow.txt", header = TRUE, row.names = 1, colClasses = "character", sep="\t")
asv_ITS_sow = read.table(file="./ITS/Culture-dependent data/asv_ITS_sow.txt", header = TRUE, row.names = 1, sep="\t")
tax_ITS_sow = read.table(file="./ITS/Culture-dependent data/tax_ITS_sow.txt", header = TRUE, row.names = 1, sep="\t")
sampledata_ITS_sow = sample_data(sampledata_ITS_sow)
tax_ITS_sow = tax_table(as.matrix(tax_ITS_sow))
asv_ITS_sow = otu_table(asv_ITS_sow, taxa_are_rows = TRUE)
physeq_ITS_sow = phyloseq(asv_ITS_sow, tax_ITS_sow, sampledata_ITS_sow)


filtered_data = as.data.frame(colSums(otu_table(physeq_ITS)))
filtered_data = as.data.frame(t(t(filtered_data)))
colnames(filtered_data) = "Sequences after filtration used in analysis"
filtered_data$SRA = row.names(filtered_data)
raw_data = read_qza("./ITS/dada2/ITS-stats-dada2.qza")$data
raw_data$SRA = row.names(raw_data)


data_ITS = left_join(raw_data,filtered_data, by="SRA")

data_ITS = data_ITS[,c("SRA", "input", "Sequences after filtration used in analysis")]
colnames(data_ITS) = c("SRA", "Raw paired-end reads", "Sequences after filtration used in analysis")
data_ITS$Sum_raw = sum(data_ITS[,"Raw paired-end reads"])
data_ITS$Mean_raw = mean(data_ITS[,"Raw paired-end reads"])
data_ITS$Median_raw = median(data_ITS[,"Raw paired-end reads"])
data_ITS$Sum_filtered = sum(data_ITS[,"Sequences after filtration used in analysis"])
data_ITS$Mean_filtered = mean(data_ITS[,"Sequences after filtration used in analysis"])
data_ITS$Median_filtered = median(data_ITS[,"Sequences after filtration used in analysis"])

sampledata_ITS = left_join(sampledata_ITS, data_ITS, by="SRA")

write.table(sampledata_ITS, file="Sampledata_ITS_with_summary.txt", quote=FALSE, sep="\t")

save.image(file='./ITS/physeq_ITS.RData')
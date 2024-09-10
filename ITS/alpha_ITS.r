library("tidyverse")
library ("phyloseq")
library("vegan")

dir = as.character(read.table("dir.txt"))
setwd(dir)
load("./ITS/physeq_ITS.RData")



alpha_diversity = estimate_richness(physeq_ITS, split = TRUE, measures = NULL)
alpha_diversity$SRA = rownames(alpha_diversity)
alpha_diversity = left_join(sampledata_ITS, alpha_diversity, by = "SRA")


# Figure 5b. Shannon′s alpha diversity among organs (ITS)
alpha_diversity$Organ_material = factor(alpha_diversity$Organ_material, level = c("Leaf","Stem","Root","Flower","Seed"))
ggplot(data = alpha_diversity, aes(x = Organ_material, y = Shannon)) + 
	geom_boxplot() +
	geom_jitter(width = 0.2) +
	theme_bw() +
	scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7), limits = c(0,7)) +
	theme(text = element_text(size=25)) +
	labs(x = "", y = "Shannon diversity measurement") +
	theme(axis.text.x = element_text(angle=75, vjust = 0.5, face = "italic")) +
	ggtitle("Shannon diversity index")
ggsave("./Figure 5b. Shannon′s alpha diversity among organs (ITS).png", width = 7, height = 10)
#


# Pairwise Wilcoxon rank sum test (organs) (ITS)
Shannon = pairwise.wilcox.test(alpha_diversity$Shannon, alpha_diversity$Organ_material, p.adjust.method = 'fdr')
print("Pairwise Wilcoxon rank sum test (organs) (ITS)")
Shannon$p.value
write.table(Shannon$p.value, sep = "\t", quote = FALSE, file="Table S7 (ITS).txt")



# Figure 8b. Shannon′s alpha diversity among seasons (ITS)
ggplot(data = alpha_diversity, aes(x = Season, y = Shannon)) + 
	geom_boxplot() +
	geom_jitter(width = 0.2) +
	theme_bw() +
	scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7), limits = c(0,7)) +
	theme(text = element_text(size=25)) +
	labs(x = "", y = "Shannon diversity measurement") +
	theme(axis.text.x = element_text(angle=75, vjust = 0.5, face = "italic")) +
	ggtitle("Shannon diversity index")
ggsave("./Figure 8b. Shannon′s alpha diversity among seasons (ITS).png", width = 7, height = 10)
#
dev.off()

# Pairwise Wilcoxon rank sum test (seasons) (ITS)
Shannon = pairwise.wilcox.test(alpha_diversity$Shannon, alpha_diversity$Season, p.adjust.method = 'fdr')
print("Pairwise Wilcoxon rank sum test (seasons) (ITS)")
Shannon


# Figure 10b. Shannon′s alpha diversity among locations (ITS)
ggplot(data = alpha_diversity, aes(x = location, y = Shannon)) + 
	geom_boxplot() +
	geom_jitter(width = 0.2) +
	theme_bw() +
	scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7), limits = c(0,7)) +
	theme(text = element_text(size=25)) +
	labs(x = "", y = "Shannon diversity measurement") +
	theme(axis.text.x = element_text(angle=75, vjust = 0.5, face = "italic")) +
	ggtitle("Shannon diversity index")
ggsave("./Figure 10b. Shannon′s alpha diversity among locations (ITS).png", width = 7, height = 10)
#
dev.off()

# Pairwise Wilcoxon rank sum test (locations) (ITS)
Shannon = pairwise.wilcox.test(alpha_diversity$Shannon, alpha_diversity$location, p.adjust.method = 'fdr')
print("Pairwise Wilcoxon rank sum test (locations) (ITS)")
Shannon




library("tidyverse")
library ("phyloseq")
library ("microViz")




mycolors = c("#018BFF", "#C02727", "#774f38", "#fa6900", "#FFD101", "#C08227",
 "#65216B", "#C177C7", "#EBBAE3", "#BAEBC4", "#397DD7", "#6FE8D9", 
"#69d2e7", "#a7dbd8", "#e0e4cc", "#f38630", "#fa6900", "#fe4365", "#fc9d9a", "#f9cdad", 
"#c8c8a9", "#83af9b", "#ecd078", "#d95b43", "#c02942", "#542437", "#53777a", 
"#4ecdc4", "#c7f464", "#ff6b6b", "#c44d58", "#774f38", "#e08e79", "#f1d4af", 
"#c5e0dc")



dir = as.character(read.table("dir.txt"))
setwd(dir)
load("./16s/physeq_16s.RData")


ibd = physeq_16s

sample_data(ibd)$Organ_material = factor(sample_data(ibd)$Organ_material, level = c("Leaf","Stem","Root","Flower","Seed"))

png("Figure 5c. PCoA plot organs (16S).png", width = 10, height = 10, res=300, units = "in")
ibd %>%
  tax_filter(tax_level = "Genus", min_prevalence = 0.1, verbose = FALSE) %>%
  tax_transform(rank = "Genus", trans = "identity", zero_replace = 0) %>%
  dist_calc(dist = "bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.6, size = 2, color = "Organ_material", shape = "Plant", auto_caption = NA) +
  scale_colour_manual(values = mycolors)+
  scale_shape_manual(values = c(15,16,17,18,0:2,5:7,9,10,12:14,3,4,8))+
  stat_ellipse(aes(colour = Organ_material))+
  theme_classic(12) +
  labs(color = "Organ material") +
  coord_fixed(0.7)
dev.off()


png("Figure 8c. PCoA plot seasons (16S).png", width = 10, height = 10, res=300, units = "in")
ibd %>%
  tax_filter(tax_level = "Genus", min_prevalence = 0.1, verbose = FALSE) %>%
  tax_transform(rank = "Genus", trans = "identity", zero_replace = 0) %>%
  dist_calc(dist = "bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.6, size = 2, color = "Season", shape = "Plant", auto_caption = NA) +
  scale_colour_manual(values = mycolors)+
  scale_shape_manual(values = c(15,16,17,18,0:2,5:7,9,10,12:14,3,4,8))+
  stat_ellipse(aes(colour = Season))+
  theme_classic(12) +
  coord_fixed(0.7)
dev.off()


png("Figure 10c. PCoA plot locations (16S).png", width = 10, height = 10, res=300, units = "in")
ibd %>%
  tax_filter(tax_level = "Genus", min_prevalence = 0.1, verbose = FALSE) %>%
  tax_transform(rank = "Genus", trans = "identity", zero_replace = 0) %>%
  dist_calc(dist = "bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.6, size = 2, color = "location", shape = "Plant", auto_caption = NA) +
  scale_colour_manual(values = mycolors)+
  scale_shape_manual(values = c(15,16,17,18,0:2,5:7,9,10,12:14,3,4,8))+
  stat_ellipse(aes(colour = location))+
  theme_classic(12) +
  labs(color = "Location") +
  coord_fixed(0.7)
dev.off()



Perm <- physeq_16s %>%
  tax_transform("identity", rank = "Genus") %>%
  dist_calc("bray") %>%
  dist_permanova(
  seed = 1234, 
  n_processes = 1, n_perms = 999, 
  variables = c("Season", "Organ_material", "location")
  )


write.table(Perm@permanova, sep = "\t", quote = FALSE, file="./Table S8 (16S).txt")
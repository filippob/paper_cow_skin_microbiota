## Rscript
## script to study the taxonomy

#############################################################################
## This script is mainly meant to be run locally
## where R packages can more easily be installed/updated
#############################################################################

###############
## SET UP
###############

library("knitr")
library("tidyr")
library("broom")
library("ggplot2")
library("egg")
library("data.table")
library("dplyr")

###############
## PARAMETERS
###############

HOME <- Sys.getenv("HOME")
project_folder = file.path(HOME, "paper_cow_skin_microbiota")
data_folder = "Data"
res_folder = "Results"
conf_file = "metadata_cow_skin_microbiota.csv"
otu_table = "otu_norm_cow_skin_microbiota.csv"

fname1 = file.path(project_folder, data_folder, conf_file)
metadata <- fread(fname1)
metadata = filter(metadata, sample != "sample-29")
metadata <- metadata [order(metadata$sample),]
meta <- metadata [,c(1,3)]
meta_cols = names(meta)

fname2 = file.path(project_folder, data_folder, otu_table)
otu <- fread(fname2)
otu <- filter(otu, Family !="Mitochondria")
otu <- filter(otu, Class !="Chloroplast")
otu <- filter(otu, Order !="Chloroplast")

###############
## Taxonomic analysis at phylum level
###############

otu_phylum=select(otu, 2:32)
otu_phylum$Phyla <- paste(otu$Phylum)
otu_phylum <- otu_phylum %>% group_by(Phyla) %>% summarise(across(everything(),sum))

writeLines(" - generating table of core microbiota at phylum level")
core_phyla <- otu_phylum%>%filter_all(all_vars(.!=0))
rownames(core_phyla) <- core_phyla$Phyla
core_phyla$Phyla <- NULL
core_phyla$avg_count <- rowMeans(core_phyla)

## making results folder
if(!file.exists(file.path(project_folder, res_folder))) dir.create(file.path(project_folder, res_folder), showWarnings = FALSE)
fwrite(x = core_phyla, file = file.path(project_folder, res_folder, "phyla_core.csv"))

## Taxonomic analysis at phylum level - relative abundances

otu_phylum=select(otu, 2:32,34)
otu_phylum <- otu_phylum %>% group_by(Phylum) %>% summarise(across(everything(),sum))
otu_norm <- select(otu_phylum, -1)
otu_norm <- otu_norm/colSums(otu_norm)
rownames(otu_norm) <- otu_phylum$Phylum
sample_names <- colnames(otu_norm)
taxa_names <- otu_phylum$Phylum
taxa_names <- as.data.frame(taxa_names)
otu_norm$`#OTU ID` <- cbind(taxa_names$taxa_names)
otu_norm <- otu_norm %>% dplyr::select(`#OTU ID`, everything())
otu_norm <- gather(otu_norm, key = "sample", value ="counts", -`#OTU ID`) %>% spread(key = `#OTU ID`, value = counts)
otu_norm <- otu_norm[order(otu_norm$sample),]
otu_norm <- cbind(otu_norm, timepoint=metadata$timepoint, sample=metadata$sample)
mO <- reshape2::melt(otu_norm,id.vars = meta_cols, value.name = "counts", variable.name = "taxa")
mO$counts <- as.numeric(mO$counts)
mO <- mO %>%
  group_by(sample) %>%
  mutate(tot = sum(counts), rel_abundance = counts/tot)
D <- mO %>%
  group_by(taxa,timepoint) %>%
  summarise(avg_abund = round(mean(rel_abundance),4), std = round(sd(rel_abundance),3))
D <- na.omit(D)
D <- D %>%
  group_by(taxa,timepoint) %>%
  # filter(avg_abund > 0.01) %>% #ONE FILTER HERE
  arrange(desc(avg_abund))
oldc <- D$taxa[D$avg_abund < 0.01]
newc <- rep("Lower than 1%", length(oldc))
vec <- newc[match(D$taxa,oldc)]
D$taxa <- ifelse(D$taxa %in% oldc, "Lower than 1%", as.character(D$taxa))
D$taxa <- gsub("group","",D$taxa)
D <- D %>% group_by(taxa, timepoint) %>% summarise(across(everything(),sum))
D <- D[order(D$avg_abund, decreasing = T),]
kable(D)

writeLines(" - printing a piechart of relative abundances through time point (filter at 1%)")
D <- D %>% arrange(desc(avg_abund))
p <- ggplot(D, aes(x=factor(1), y=avg_abund, fill=taxa), palette(mycolors)) + geom_bar(width=1,stat="identity", position = "fill") + 
  facet_grid(~timepoint) +
  guides(fill = guide_legend(title = "Phyla")) +
  coord_polar(theta='y') +
  xlab("Relative abundances") + ylab("Percentages")+
  scale_fill_brewer(palette = "Set3")+
  theme(text = element_text(size=20),
               axis.text = element_blank(),
               axis.ticks = element_blank(),
               legend.text=element_text(size=20),
               legend.title=element_text(size=20)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
p
ggsave(filename = file.path(project_folder, res_folder, "phyla_rel.abund.png"), plot = p, device = "png", width = 15, height = 9)

## Taxonomic analysis at phylum level - Significance of time point
mO$level <- "phylum"
m1 <- mO %>%
  filter(!is.na(level)) %>%
  arrange(level,timepoint)
save(m1, file = file.path(project_folder, res_folder, "phyla_significance.RData"))

dd_counts <- mO %>%
  group_by(level, taxa,timepoint) %>%
  summarise(avg = mean(rel_abundance)) %>%
  spread(key = "timepoint", value = "avg")

group_by(dd_counts, level) 
D <- m1 %>%
  group_by(level, taxa) %>%
  do(tidy(anova(lm(counts ~ timepoint, data = .)))) %>%
  filter(term == "timepoint")

D$level  <- factor(D$level,levels = c("phylum"))
D <- D %>%
  arrange(level,taxa)
DX <- D %>%
  filter(`p.value` <= 0.05) %>% 
  dplyr::select(c(level,taxa, `p.value`)) %>%
  arrange(level,`p.value`)
D0 <- mO %>%
  dplyr::group_by(level,taxa, timepoint) %>%
  dplyr::summarise(avg_counts = mean(counts))
to_save = list(D,DX,D0)

save(to_save, file = file.path(project_folder, res_folder, "phyla_taxonomy.RData"))

load (file.path(project_folder, res_folder, "phyla_taxonomy.RData"))
D <- to_save[[1]]
DX <- to_save[[2]]
D0 <- to_save[[3]]

dd <- spread(D0, key = timepoint, value = avg_counts)
temp <- inner_join(DX,dd, by = c("level" = "level", "taxa" = "taxa"))
fwrite(temp, file = file.path(project_folder, res_folder, "phyla_significant_otus.csv"), col.names = TRUE, sep = ",")
print (dd)

load(file.path(project_folder, res_folder, "phyla_taxonomy.RData"))
D <- to_save[[1]]
DX <- to_save[[2]]
D0 <- to_save[[3]]

D0 <- mutate(D0, avg_counts = avg_counts+1) %>% spread(key = timepoint, value = avg_counts)
D1 <- DX %>%
  inner_join(D0, by = c("level" = "level", "taxa" = "taxa")) %>%
  mutate(p.value = -log10(p.value)) %>%
  gather(key = "treatment", value = "counts", -c(level,taxa, p.value))
D1$level <- factor(D1$level, levels = c("phylum")) #, "species"
D1 <- D1 %>% group_by(level) %>% mutate(tot = sum(counts), relab = counts/tot)

## Extracting estimated coefficients for "significant" OTUs

res = data.frame("level"=NULL, "OTU"=NULL, "p_value"=NULL, "T0"=NULL,"T1"=NULL, "T2"=NULL)
DX <- DX %>% filter(level == "phylum")
m1 <- m1 %>% filter(level=="phylum")

for (name in DX$taxa) {
  
  print(paste("analysing OTU ", name))
  pval = as.numeric(DX[DX$taxa==name,"p.value"])
  level = DX[DX$taxa==name,"level"]
 
  ## estimating coefficients
  temp2 = filter(m1, taxa == name)
  temp2$timepoint <- factor(temp2$timepoint, levels = c("T0","T1", "T2"))
  g = lm(counts ~ timepoint, data = temp2)
  
  ## extracting coefficients
  coefs = g$coefficients
  coefs = coefs[!grepl("(Intercept)", names(coefs))]
  names(coefs) = gsub("treatment","",names(coefs))
  coefs = as.data.frame(t(coefs))
  
  ## adding metadata
  coefs["level"] = level
  coefs["OTU"] = name
  coefs["p_value"] = pval
 
  # saving results
  res = rbind.data.frame(res, coefs)
}

comp <- temp[,c(2,4:6)]
comp <- comp %>% pivot_longer(cols = T0:T2, names_to = "timepoint", names_prefix = "", values_to = "value")
q <- ggplot(temp, aes(x = factor(1), y = taxa, height=0.95)) + 
  geom_tile(aes(fill = p.value), colour = "white") +
  scale_fill_gradient(low = "#1aff1a", high = "#4b0092", limits=c(0,0.05)) + 
  xlab("pvalue") + ylab("Phyla") + 
  theme(axis.text.x = element_text(size = 19), axis.text.y = element_text(size = 13), axis.ticks = element_blank(), legend.position = "left", axis.title.x = element_blank(), axis.title.y = element_text(size=14))
q

comp$log <- -log10(comp$value)
px <- ggplot(comp, aes(x=taxa, y=log, fill=timepoint)) + 
  geom_bar(width=0.95,stat="identity", position=position_dodge()) + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  coord_flip() + 
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC", "#e66100", "#40b0A6" )) + 
  ylab("-log10(count)")
px

join <- ggarrange(q, px, widths=c(0.5,1), heights=c(1,1), labels = c("A","B"))
print(join)
ggsave(filename = file.path(project_folder, res_folder, "phyla_pvalues.png"), plot = join, device = "png", width = 15, height = 9)

###############
## Taxonomic analysis at genus level
###############

otu_genus=select(otu, 2:32)
otu_genus$Genera <- paste(otu$Phylum, otu$Class, otu$Order, otu$Family, otu$Genus, sep = "; ")
otu_genus <- otu_genus %>% group_by(Genera) %>% summarise(across(everything(),sum))

writeLines(" - generating table of core microbiota at genus level")
core_genera <- otu_genus%>%filter_all(all_vars(.!=0))
core_genera$avg_count <- rowMeans(subset(core_genera, select = c(2:32)))
fwrite(x = core_genera, file = file.path(project_folder, res_folder, "genera_core.csv"))

## Taxonomic analysis at genera level - relative abundances
otu_genus=select(otu, 2:32, 38)
otu_genus <- otu_genus %>% group_by(Genus) %>% summarise(across(everything(),sum))
uncult <- slice(otu_genus, 1, 481:498)
uncult$Genus <- "Uncultured or unknown"
uncult <- uncult %>%
  group_by(Genus) %>%
  summarise(across(everything(), sum))
otu_genus <- otu_genus[-c(1,  1, 481:498), ]
otu_genus <- rbind(otu_genus, uncult)
otu_norm <- select(otu_genus, -1)
otu_norm <- otu_norm/colSums(otu_norm)
rownames(otu_norm) <- otu_genus$Genus

writeLines(" - generating table of core microbiota at phylum level")
core_genera <- otu_norm%>%filter_all(all_vars(.!=0))
sample_names <- colnames(otu_norm)
taxa_names <- otu_genus$Genus
taxa_names <- as.data.frame(taxa_names)
otu_norm$`#OTU ID` <- cbind(taxa_names$taxa_names)
otu_norm <- otu_norm %>% dplyr::select(`#OTU ID`, everything())
otu_norm <- gather(otu_norm, key = "sample", value ="counts", -`#OTU ID`) %>% spread(key = `#OTU ID`, value = counts)
otu_norm <- otu_norm[order(otu_norm$sample),]
otu_norm <- cbind(otu_norm, timepoint=metadata$timepoint, sample=metadata$sample)
mO <- reshape2::melt(otu_norm,id.vars = meta_cols, value.name = "counts", variable.name = "taxa")

## Taxonomic analysis at genus level - relative abundances

mO$counts <- as.numeric(mO$counts)
mO <- mO %>%
  group_by(sample) %>%
  mutate(tot = sum(counts), rel_abundance = counts/tot)
D <- mO %>%
  group_by(taxa,timepoint) %>%
  summarise(avg_abund = round(mean(rel_abundance),4), std = round(sd(rel_abundance),3))
D <- na.omit(D)
D <- D %>%
  group_by(taxa,timepoint) %>%
  # filter(avg_abund > 0.01) %>% #ONE FILTER HERE
  arrange(desc(avg_abund))

oldc <- D$taxa[D$avg_abund < 0.02]
newc <- rep("Lower than 2%", length(oldc))
vec <- newc[match(D$taxa,oldc)]
D$taxa <- ifelse(D$taxa %in% oldc, "Lower than 2%", as.character(D$taxa))
D$taxa <- gsub("group","",D$taxa)

D <- D %>% group_by(taxa, timepoint) %>% summarise(across(everything(),sum))
kable(D)

writeLines(" - printing a piechart of relative abundances through time point (filter at 2%)")
p <- ggplot(D, aes(x=factor(1), y=avg_abund, fill=taxa)) + geom_bar(width=1,stat="identity", position = "fill") +
  facet_grid(~timepoint) +
  guides(fill = guide_legend(title = "Genera")) + coord_polar(theta='y') +
  xlab("Relative abundances") + ylab("Percentages") +
  scale_fill_brewer(palette = "Paired") +
  theme(text = element_text(size=20),
               axis.text = element_blank(),
               axis.ticks = element_blank(),
               legend.text=element_text(size=20),
               legend.title=element_text(size=20)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
ggsave(filename = file.path(project_folder, res_folder, "genera_rel.abund.png"), plot = p, device = "png", width = 15, height = 9)

## Taxonomic analysis at genus level - Significance of timepoint

mO$level <- "genus"
m1 <- mO %>%
  filter(!is.na(level)) %>%
  arrange(level,timepoint)
save(m1, file = file.path(project_folder, res_folder, "genera_significance.RData"))
dd_counts <- mO %>%
  group_by(level, taxa,timepoint) %>%
  summarise(avg = mean(rel_abundance)) %>%
  spread(key = "timepoint", value = "avg")
group_by(dd_counts, level)
D <- m1 %>%
  group_by(level, taxa) %>%
  do(tidy(anova(lm(counts ~ timepoint, data = .)))) %>%
  filter(term == "timepoint")
D$level  <- factor(D$level,levels = c("genus")) 
D <- D %>%
  arrange(level,taxa)
DX <- D %>%
  filter(`p.value` <= 0.05) %>% 
  dplyr::select(c(level,taxa, `p.value`)) %>%
  arrange(level,`p.value`)
D0 <- mO %>%
  dplyr::group_by(level,taxa, timepoint) %>%
  dplyr::summarise(avg_counts = mean(counts))
to_save = list(D,DX,D0)
save(to_save, file = file.path(project_folder, res_folder, "genera_taxonomy.RData"))

load(file.path(project_folder, res_folder, "genera_taxonomy.RData"))
D <- to_save[[1]]
DX <- to_save[[2]]
D0 <- to_save[[3]]

dd <- spread(D0, key = timepoint, value = avg_counts)
temp <- inner_join(DX,dd, by = c("level" = "level", "taxa" = "taxa"))
fwrite(temp, file = file.path(project_folder, res_folder, "genera_significant_otus.csv"), col.names = TRUE, sep = ",")
print (dd)

load(file.path(project_folder, res_folder, "genera_taxonomy.RData"))
D <- to_save[[1]]
DX <- to_save[[2]]
D0 <- to_save[[3]]
D0 <- mutate(D0, avg_counts = avg_counts+1) %>% spread(key = timepoint, value = avg_counts)
D1 <- DX %>%
  inner_join(D0, by = c("level" = "level", "taxa" = "taxa")) %>%
  mutate(p.value = -log10(p.value)) %>%
  gather(key = "treatment", value = "counts", -c(level,taxa, p.value))
D1$level <- factor(D1$level, levels = c("phylum")) #, "species"
D1 <- D1 %>% group_by(level) %>% mutate(tot = sum(counts), relab = counts/tot)

## Extracting estimated coefficients for "significant" OTUs
res = data.frame("level"=NULL, "OTU"=NULL, "p_value"=NULL, "T0"=NULL,"T1"=NULL, "T2"=NULL)
DX <- DX %>% filter(level == "genus")
m1 <- m1 %>% filter(level=="genus")

for (name in DX$taxa) {
  
  print(paste("analysing OTU ", name))
  pval = as.numeric(DX[DX$taxa==name,"p.value"])
  level = DX[DX$taxa==name,"level"]
 
  ## estimating coefficients
  temp2 = filter(m1, taxa == name)
  temp2$timepoint <- factor(temp2$timepoint, levels = c("T0","T1", "T2"))
  g = lm(counts ~ timepoint, data = temp2)
  
  ## extracting coefficients
  coefs = g$coefficients
  coefs = coefs[!grepl("(Intercept)", names(coefs))]
  names(coefs) = gsub("treatment","",names(coefs))
  coefs = as.data.frame(t(coefs))
  
  ## adding metadata
  coefs["level"] = level
  coefs["OTU"] = name
  coefs["p_value"] = pval
 
  # saving results
  res = rbind.data.frame(res, coefs)
}

comp <- temp[,c(2,4:6)]
comp <- comp %>% pivot_longer(cols = T0:T2, names_to = "timepoint", names_prefix = "", values_to = "value")
q <- ggplot(temp, aes(x = factor(1), y = taxa, height=0.95)) +
  geom_tile(aes(fill = p.value), colour = "white") +
  scale_fill_gradient(low = "#1aff1a", high = "#4b0092", limits=c(0,0.05)) +
  xlab("pvalue") + ylab("Genera") +
  theme(axis.text.x = element_text(size = 19), axis.text.y = element_text(size = 13), axis.ticks = element_blank(), legend.position = "left", axis.title.x = element_blank(), axis.title.y = element_text(size=14))
q

comp$log <- -log10(comp$value)
px <- ggplot(comp, aes(x=taxa, y=log, fill=timepoint)) + 
  geom_bar(width=0.95,stat="identity", position=position_dodge()) + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) + 
  coord_flip() + scale_fill_manual(values = c("#FFC20A", "#0C7BDC", "#e66100", "#40b0A6" )) + ylab("-log10(count)")
px

join <- ggarrange(q, px, widths=c(0.5,1), heights=c(1,1), labels = c("A","B"))
print(join)
ggsave(filename = file.path(project_folder, res_folder, "genera_pvalues.png"), plot = join, device = "png", width = 15, height = 9)
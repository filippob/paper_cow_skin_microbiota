## Rscript
## script to normalize the filtered OTU table

#############################################################################
## This script is mainly meant to be run locally
## where R packages can more easily be installed/updated
#############################################################################

###############
## SET UP
###############

library("phyloseq")
library("tidyverse")
library("data.table")
library("metagenomeSeq")

###############
## PARAMETERS
###############

HOME <- Sys.getenv("HOME")
project_folder = file.path(HOME, "paper_cow_skin_microbiota")
data_folder = "Data"
outdir = file.path(data_folder)
support = "support_functions"
conf_file = "metadata_cow_skin_microbiota.csv"
otu_table = "otu_table_cow_skin_microbiota.biom"
min_tot_counts = 500 ## minimum number of total counts per sample to be included in the analysis

source(file.path(project_folder, support, "dist2list.R")) ## from: https://github.com/vmikk/metagMisc/
source(file.path(project_folder, support, "phyloseq_transform.R")) ## from: https://github.com/vmikk/metagMisc/

###############
## Reading OTU table in biom format and metadata
###############

writeLines(" - reading the filtered (OTU-wise) biom file into phyloseq")
biom_otu_tax <- phyloseq::import_biom(BIOMfilename = file.path(project_folder,data_folder,otu_table))

writeLines(" - removing samples with too few total counts")
biom_otu_tax = prune_samples(sample_sums(biom_otu_tax)>=min_tot_counts, biom_otu_tax)

otu = otu_table(biom_otu_tax, taxa_are_rows = TRUE)
taxa = tax_table(biom_otu_tax)

print(paste("N. of OTUs read from biom file is:", nrow(otu)))
print(paste("N .of samples retained after filtering is:", ncol(otu)))

colnames(otu) <- paste("sample-",colnames(otu),sep="")
print(head(otu))

writeLines(" - change the names of taxonomic levels to Kngdom, Class etc.")
colnames(taxa) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species") #if number does not fit, add "" as blank spaces to solve the problem
print(head(taxa))

## metadata
writeLines(" - reading the metadata")
metadata = fread(file.path(project_folder, data_folder,conf_file))
metadata <- as.data.frame(metadata)
row.names(metadata) <- metadata$sample
metadata$sample <- NULL

## read into phyloseq
writeLines(" - add metadata to the phyloseq object")
samples = sample_data(metadata)
otu_tax_sample = phyloseq(otu,taxa,samples)
sample_data(otu_tax_sample)

###############
## Alpha diversity
###############

## alpha diversity is calculated on the original count data, not normalised 
## (see https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis)
## (see https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531)

writeLines(" - calculate alpha diversity indices")
alpha = estimate_richness(otu_tax_sample, split = TRUE)
alpha$sample = row.names(alpha)
alpha = relocate(alpha, sample)
alpha$sample <-  gsub("sample.","sample-", as.character(alpha$sample))
fwrite(x = alpha, file = file.path(project_folder, data_folder, "alpha_div_cow_skin_microbiota.csv"))

writeLines(" - CSS normalization")
otu_tax_sample_norm = phyloseq_transform_css(otu_tax_sample, norm = TRUE, log = FALSE)

## add taxonomy to normalised counts
otu_css_norm = base::as.data.frame(otu_table(otu_tax_sample_norm))
otu_css_norm$tax_id = row.names(otu_css_norm)
otu_css_norm <- relocate(otu_css_norm, tax_id)
taxonomy = as.data.frame(tax_table(otu_tax_sample_norm))
taxonomy$tax_id = row.names(taxonomy)
taxonomy <- relocate(taxonomy, tax_id)
otu_css_norm = otu_css_norm %>% inner_join(taxonomy, by = "tax_id")

writeLines(" - writing out the CSS normalized OTU table")
fwrite(x = otu_css_norm, file = file.path(project_folder, outdir, "otu_norm_cow_skin_microbiota.csv"))

## relative abundances
otu_relative = transform_sample_counts(otu_tax_sample, function(x) x/sum(x) )
otu_rel_filtered = filter_taxa(otu_relative, function(x) mean(x) > 5e-3, TRUE)
nrow(otu_table(otu_rel_filtered))

###############
## distances 
###############

writeLines(" - beta diversity: distance matrices")
writeLines(" - available distance metrics")
dist_methods <- unlist(distanceMethodList)
print(dist_methods)

## bray-curtis
writeLines(" - calculate Bray-Curtis distances")
distances = distance(otu_tax_sample_norm, method="bray", type = "samples")
iMDS  <- ordinate(otu_tax_sample_norm, "MDS", distance=distances)
writeLines(" - write out distance matrix")
dd = dist2list(distances, tri = FALSE)
dx = spread(dd, key = "col", value = "value")
fwrite(x = dx, file = file.path(project_folder, data_folder, "bray_curtis_distances_cow_skin_microbiota.csv"))
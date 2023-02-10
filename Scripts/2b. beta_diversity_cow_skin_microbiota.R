## Rscript
## Script to run beta diversity analysis

#############################################################################
## This script is mainly meant to be run locally
## where R packages can more easily be installed/updated
#############################################################################

###############
## SET UP
###############

library("dplyr")
library("ggplot2")
library("data.table")

###############
## Parameters
###############

HOME <- Sys.getenv("HOME")
project_folder = file.path(HOME, "paper_cow_skin_microbiota")
data_folder = "Data"
analysis_folder = "results"
conf_file = "metadata_cow_skin_microbiota.csv"
bray_curtis = "bray_curtis_distances_cow_skin_microbiota.csv"
fname1 = file.path(project_folder, data_folder, conf_file)

metadata <- fread(fname1)
metadata = filter(metadata, sample != "sample-29")

fname2 = file.path(project_folder, data_folder, bray_curtis)
matrice <- read.table(fname2, row.names = 1, sep = ",") 
colnames(matrice) <- matrice[1,] 
matrice <- matrice[-1,]

vec <- select(metadata, "sample") %>% pull()
vex <- names(matrice) %in% vec
mat_milk = matrice[vex,vex]
mat_milk$timepoint <- as.character(metadata$timepoint[match(row.names(mat_milk),metadata$sample)])
matx= data.matrix(select(mat_milk, -c(timepoint)))

###############
## MDS
###############

mds <- cmdscale(as.dist(matx))
mds <- as.data.frame(mds)

mds$timepoint <- metadata$timepoint[match(rownames(mds), metadata$sample)]
p <- ggplot(mds, aes(V1,V2)) + geom_point(aes(colour = timepoint), size = 3) + stat_ellipse(aes(x=V1, y=V2,color=timepoint), type="norm")
p <- p + xlab("dim1") + ylab("dim2")
p

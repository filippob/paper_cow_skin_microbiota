## Rscript
## Script to run alpha diversity analysis

#############################################################################
## This script is mainly meant to be run locally
## where R packages can more easily be installed/updated
#############################################################################


###############
## SET UP
###############

library("DT")
library("dplyr")
library("tidyr")
library("broom")
library("scales")
library("data.table")

###############
## Parameters
###############

HOME <- Sys.getenv("HOME")
project_folder = file.path(HOME, "paper_cow_skin_microbiota")
data_folder = "Data"
analysis_folder = "results"
conf_file = "metadata_cow_skin_microbiota.csv"
alpha = "alpha_div_cow_skin_microbiota.csv"
fname1 = file.path(project_folder, data_folder, conf_file)

metadata <- fread(fname1)
metadata = filter(metadata, sample != "sample-29")

fname2 = file.path(project_folder, data_folder, alpha)
alpha <- fread(fname2)
names(alpha)[1] <- "sample"
alpha$`sample` <- gsub("\\.", "-", alpha$`sample`)
alpha = select(alpha, -c(se.chao1, se.ACE))

###############
## Alpha diversity
###############

mAlpha <- reshape2::melt(alpha, id.vars = "sample", variable.name = "metric", value.name = "value")
mAlpha$timepoint <- metadata$timepoint[match(mAlpha$sample,metadata$sample)]

C <- mAlpha %>%
  group_by(metric, timepoint) %>%
  summarize(N=n(),avg=round(mean(value),3)) %>%
  spread(key = metric, value = avg)

D <- mAlpha %>%
  group_by(metric) %>%
  do(tidy(lm(value ~ timepoint, data = .))) %>%
  filter(term != "(Intercept)")

#indices <- colnames(alpha)[-1]
#for (k in indices)

fit <- lm(value~timepoint, filter(mAlpha, metric == "InvSimpson")) #replace with index name
summary(fit) #retrieve overall pvalue and single index per timepoint pvalue

datatable(D, options = list(pageLength=100)) %>% 
  formatStyle('p.value', backgroundColor = styleInterval(0.05, c('yellow', 'white')))

mAlpha <- mAlpha %>% group_by(metric) %>%
  mutate(scaled_value = rescale(value, to = c(0,100)))

p <- ggplot(data = D, mapping= aes(x=term, y=p.value)) +
  geom_point(aes(color = metric, stroke = 1), position=position_jitter(h=0, w=0.27)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", linewidth=0.5) +
  geom_hline(yintercept=0.10, linetype="dashed", color = "darkorange", linewidth=0.5) +
  coord_trans(y="log2") +
  scale_y_continuous(breaks=pretty_breaks(n=20)) + 
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 1)) + 
  theme(axis.text.x = element_text(angle=90))
p

###############
## Bootstrapping
###############

boot_sample = function(data,index) {
  n = nrow(data)
  vec = sample(1:n, n, replace = TRUE)
  temp = data[vec,]
  return(temp) }

indices <- colnames(alpha)[-1]

res = data.frame("index"=NULL, "stat"=NULL, "pval"=NULL, "coef"=NULL, "timepoint"=NULL)
for (k in indices) {
  mm = mAlpha[mAlpha$metric == k,]
  for (i in 1:10) {
    print(paste("bootstrap replicate n.", i))
    temp = boot_sample(mm, k)
    tbl = temp %>%
      group_by(timepoint) %>%
      summarise(N=n()) %>%
      spread(key = timepoint, value = N)
    tmp <- temp %>%
      group_by(metric) %>%
      do(tidy(lm(scaled_value ~ timepoint, data = .))) %>%
      filter(term != "(Intercept)") 
    res = bind_rows(res, tmp) }}

s <- ggplot(data = res, mapping= aes(x=metric, y=p.value)) +
  geom_boxplot(aes(fill=metric)) + 
  facet_wrap(~term) +
  theme(axis.text.x = element_text(angle=90))
s

D = res %>%
  group_by(term, metric) %>%
  summarise(avg = median(p.value), std = sd(p.value), wt = 1/std)

px <- ggplot(data = D, mapping= aes(x=term, y=avg)) +
  geom_point(aes(color = metric, size = wt), position=position_jitter(h=0, w=0.27)) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=0.5) +
  geom_hline(yintercept=0.10, linetype="dashed", color = "darkorange", size=0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_trans(y="log2") +
  scale_y_continuous(breaks=pretty_breaks(n=20)) +
  scale_color_manual(values = c("#f8766d", "#d39200", "#93aa00", "#00ba38", "#00c19f", "#00b9e3", "#619cff", "#db72fb", "#ff61c3"))  +
  scale_shape_manual(values = c(0, 1, 2, 23, 3,4,7, 8, 10))  +
  theme(legend.key.size = unit(0.05, 'cm'), axis.text.x = element_text(angle=90)) +
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 1))
px
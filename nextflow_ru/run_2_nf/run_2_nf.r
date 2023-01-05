# This is the processing of Ru's RNA-seq data
# See run_2.nf for the preprocessing of the read data. I simply ran it through fastp and then kallisto.
# The transcriptome was downloaded form here: https://www.ncbi.nlm.nih.gov/genome/?term=txid556484[Organism:noexp]
# I have created a tx2gene file using: gffread -E GCF_000150955.2_ASM15095v2_genomic.gff --table @id,@geneid > tx2gene.raw
# That I got from here: https://www.biostars.org/p/9513656/#9513659
# I then further processed this file with: awk -F "\t" 'BEGIN{print "TXNAME\tGENEID"} {gsub("rna-", "", $1); gsub("gene-", "", $2); print $1 "\t" $2}' tx2gene.raw > tx2gene.txt
# To produce the final tx2gene file

# I also created a samples meta file: /home/humebc/projects/ru/nextflow_ru/run_1_nf/samples.csv

library(dplyr)
library(stringr)
library(tximport)
library(DESeq2)
library("pheatmap")
library(ggplot2)
library(ggvenn)

# If the script has already been run through then you can use this load.
# load(file = "nextflow_ru/run_2_nf/run_1_data.RData")

tx2gene = read.table("/home/humebc/projects/ru/reference/tx2gene.txt", header=TRUE)
samples = read.csv("/home/humebc/projects/ru/nextflow_ru/run_2_nf/samples_run_2.csv", header=TRUE)
samples = samples %>% mutate(axenic = as.factor(axenic), time_hr = as.factor(time_hr), shaking = as.factor(shaking), cell_line = as.factor(cell_line)) %>% mutate(axenic = relevel(axenic, "TRUE"))
# Create a column that is 'mutant' so that we can easily do the comparisons that ru made for the shaking samples
samples = samples %>% mutate(mutant = as.factor(case_when(cell_line == "WT" ~ FALSE, TRUE ~ TRUE)))


# Make a vector that contains the full paths to the abundance.h5 files
kallisto.base.dir = "/home/humebc/projects/ru/nextflow_ru/run_2_nf/results/kallisto"


# There are two initial results that we want to recapitulate

### SHAKING mutant vs WT for each of 0, 3 and 24h ###
# The first is the number of differential genes for each of the time points
# for mutant vs axenic

# We will want to run a DESEQ2 with a subset of the samples to get the DE for each of the time points
shaking_times = c(0, 3, 24)

# To collect the results of how many genes are up and down regulated for
# each of the time points we will make some empty vectors
contrast_time = c()
num_genes = c()
up_down = c()

# The loop for doing each of the subsets
for (time in shaking_times){
    samples_sub = samples %>% dplyr::filter(time_hr==time & shaking!="FALSE")

    files <- file.path(kallisto.base.dir, samples_sub$dir_name, "abundance.h5")

    # Finally we can use tximport to read in the abundance tables
    # and perform the normalizations
    txi = tximport(files, type = "kallisto", tx2gene = tx2gene)

    # Create the DESEQ2 object
    dds = DESeqDataSetFromTximport(txi, colData = samples_sub, design = ~ mutant)

    # Filter out those genes with <10 counts in more than 1/4 of the samples
    keep <- rowSums(counts(dds) >= 10) >= ceiling(dim(samples_sub)[[1]]/4)
    dds <- dds[keep,]

    # Fit the model and run the DE analysis
    dds = DESeq(dds)
    
    res = results(dds)
    
    up = as.data.frame(res) %>% dplyr::filter(log2FoldChange > 1 & padj < 0.05)
    contrast_time = append(contrast_time, time); num_genes = append(num_genes, dim(up)[[1]]); up_down = append(up_down, "up");
    
    down = as.data.frame(res) %>% dplyr::filter(log2FoldChange < -1 & padj < 0.05)
    contrast_time = append(contrast_time, time); num_genes = append(num_genes, dim(down)[[1]]); up_down = append(up_down, "down");
    

}


# Create the df that we will use for plotting
plotting_df_shaking = data.frame(contrast_time=as.factor(contrast_time), num_genes=num_genes, up_down=up_down)

# Plot up the results in the same way as Ru.
# I.e. stacked bar plot, one stack for each time with up and down regulated
# gene counts for each bar.
ggplot(plotting_df_shaking, aes(fill=up_down, y=num_genes, x=contrast_time, label=num_genes)) + 
    geom_bar(position="stack", stat="identity") + geom_text(size = 10, position = position_stack(vjust = 0.5)) + theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20), plot.title = element_text(size = 40), legend.text = element_text(size=25)) + ggtitle("WT vs Mutant shaking") +
  scale_fill_manual(values=c("up" = "#404040", "down" = "#AFABAB"))
ggsave("nextflow_ru/run_2_nf/rna_2_shaking_stacked_bars.filtered.png")



### SHAKING PCA ###
# Filter down the samples
# Next we want to produce a PCA for all of the samples. This means we'll
# want to make a DESeq2 object that contains all of the samples
samples_shaking = samples %>% dplyr::filter(shaking!="FALSE")
files_shaking <- file.path(kallisto.base.dir, samples_shaking$dir_name, "abundance.h5")

# Finally we can use tximport to read in the abundance tables
# and perform the normalizations
txi_shaking = tximport(files_shaking, type = "kallisto", tx2gene = tx2gene)

# Create the DESEQ2 object
dds_shaking = DESeqDataSetFromTximport(txi_shaking, colData = samples_shaking, design = ~ mutant)

# Fit the model and run the DE analysis
dds_shaking = DESeq(dds_shaking)

# Filter out those genes with <10 counts in more than 1/4 of the samples
keep_shaking <- rowSums(counts(dds_shaking) >= 10) >= ceiling(dim(samples_shaking)[[1]]/4)
dds_shaking <- dds_shaking[keep_shaking,]

vsd_shaking <- vst(dds_shaking, blind=FALSE)
rld_shaking <- rlog(dds_shaking, blind=FALSE)
head(assay(vsd_shaking), 3)

pcaData_shaking = plotPCA(vsd_shaking, intgroup=c("time_hr", "cell_line"), returnData=TRUE)

percentVar_shaking <- round(100 * attr(pcaData_shaking, "percentVar"))

ggplot(pcaData_shaking, aes(PC1, PC2, color=cell_line, shape=time_hr)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar_shaking[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_shaking[2],"% variance")) + 
  coord_fixed() + ggtitle("2nd RNA-seq shaking PCA") + scale_color_manual(values=c("A7" = "#082BD2", "A8" = "#65D527", "WT" = "#000000"))

ggsave("nextflow_ru/run_2_nf/rna_2_shaking_pca.filtered.png")
# This PCA looks quite different to the one that Ru produced.


### NON-SHAKING ###

# Next we want to do the same but for the non-shaking samples
# This includes making a stacked bar plot for the WT co-cultured vs WT axenic for each of the time points
# And for the Mutant vs wild type co-cultured.
# To collect the results of how many genes are up and down regulated for
# each of the time points we will make some empty vectors
contrast_non_shaking = c()
num_genes_non_shaking = c()
up_down_non_shaking = c()

# The loop for doing each of the subsets
for (time in c(24, 48)){
  # First we are interested in the WT co-cultured vs axenic
  samples_sub = samples %>% dplyr::filter(time_hr==time & shaking!="TRUE" & cell_line == "WT")
  files <- file.path(kallisto.base.dir, samples_sub$dir_name, "abundance.h5")
  txi = tximport(files, type = "kallisto", tx2gene = tx2gene)
  dds = DESeqDataSetFromTximport(txi, colData = samples_sub, design = ~ axenic)
  # Filter out those genes with <10 counts in more than 1/4 of the samples
  keep <- rowSums(counts(dds) >= 10) >= ceiling(dim(samples_sub)[[1]]/4)
  dds <- dds[keep,]
  dds = DESeq(dds)
  res = results(dds)

  up = as.data.frame(res) %>% dplyr::filter(log2FoldChange > 2 & padj < 0.05)
  contrast_non_shaking = append(contrast_non_shaking, paste0(time, "h_WT_CO_vs_WT_AX")); num_genes_non_shaking = append(num_genes_non_shaking, dim(up)[[1]]); up_down_non_shaking = append(up_down_non_shaking, "up");

  down = as.data.frame(res) %>% dplyr::filter(log2FoldChange < -2 & padj < 0.05)
  contrast_non_shaking = append(contrast_non_shaking, paste0(time, "h_WT_CO_vs_WT_AX")); num_genes_non_shaking = append(num_genes_non_shaking, dim(down)[[1]]); up_down_non_shaking = append(up_down_non_shaking, "down");



  # Then we will be interested in the WT co-cultured vs the mutants co-cultured
  samples_sub = samples %>% dplyr::filter(time_hr==time & shaking!="TRUE")
  files <- file.path(kallisto.base.dir, samples_sub$dir_name, "abundance.h5")
  txi = tximport(files, type = "kallisto", tx2gene = tx2gene)
  dds = DESeqDataSetFromTximport(txi, colData = samples_sub, design = ~ mutant)
  # Filter out those genes with <10 counts in more than 1/4 of the samples
  keep <- rowSums(counts(dds) >= 10) >= ceiling(dim(samples_sub)[[1]]/4)
  dds <- dds[keep,]
  dds = DESeq(dds)
  res = results(dds)

  up = as.data.frame(res) %>% dplyr::filter(log2FoldChange > 1 & padj < 0.05)
  contrast_non_shaking = append(contrast_non_shaking, paste0(time, "h_mutants_vs_WT")); num_genes_non_shaking = append(num_genes_non_shaking, dim(up)[[1]]); up_down_non_shaking = append(up_down_non_shaking, "up");

  down = as.data.frame(res) %>% dplyr::filter(log2FoldChange < -1 & padj < 0.05)
  contrast_non_shaking = append(contrast_non_shaking, paste0(time, "h_mutants_vs_WT")); num_genes_non_shaking = append(num_genes_non_shaking, dim(down)[[1]]); up_down_non_shaking = append(up_down_non_shaking, "down");
}

# Create the df that we will use for plotting
plotting_df_non_shaking = data.frame(contrast_non_shaking=as.factor(contrast_non_shaking), num_genes_non_shaking=num_genes_non_shaking, up_down_non_shaking=up_down_non_shaking)
# Reorder the columns so that they are in the same order as Ru's output
plotting_df_non_shaking = plotting_df_non_shaking %>% mutate(contrast_non_shaking = factor(contrast_non_shaking, levels=c("24h_WT_CO_vs_WT_AX", "24h_mutants_vs_WT", "48h_WT_CO_vs_WT_AX", "48h_mutants_vs_WT")))


# Plot up the results in the same way as Ru.
# I.e. stacked bar plot, one stack for each time with up and down regulated
# gene counts for each bar.
ggplot(plotting_df_non_shaking, aes(fill=up_down_non_shaking, y=num_genes_non_shaking, x=contrast_non_shaking, label=num_genes_non_shaking)) + 
    geom_bar(position="stack", stat="identity") + geom_text(size = 10, position = position_stack(vjust = 0.5)) + theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20), plot.title = element_text(size = 40), legend.text = element_text(size=25), axis.text.x = element_text(angle = -45, vjust=-0/5, hjust=0)) + ggtitle("non shaking") +
        scale_fill_manual(values=c("up" = "#404040", "down" = "#AFABAB"))

ggsave("nextflow_ru/run_2_nf/rna_2_non_shaking_stacked_bars.filtered.png")

### NON-SHAKING PCA 
# Filter down the samples 
# Next we want to produce a PCA for all of the samples. This means we'll
# want to make a DESeq2 object that contains all of the samples
samples_non_shaking = samples %>% dplyr::filter(shaking!="TRUE")

# We need to create a new variable that combines the time and the axenic state
# Instead of using the TRUE or FALSE of the of the axenic we will
# work with 'AX' and 'CO' in compliance with what RU has.
samples_non_shaking = samples_non_shaking %>% mutate(treatment = as.factor(str_c(case_when(axenic == "TRUE" ~ "AX", TRUE ~ "CO"), as.character(time_hr))))

files_non_shaking <- file.path(kallisto.base.dir, samples_non_shaking$dir_name, "abundance.h5")

# Finally we can use tximport to read in the abundance tables
# and perform the normalizations
txi_non_shaking = tximport(files_non_shaking, type = "kallisto", tx2gene = tx2gene)

# Create the DESEQ2 object
dds_non_shaking = DESeqDataSetFromTximport(txi_non_shaking, colData = samples_non_shaking, design = ~ mutant)

# Filter out those genes with <10 counts in more than 1/4 of the samples
keep_non_shaking <- rowSums(counts(dds_non_shaking) >= 10) >= ceiling(dim(samples_non_shaking)[[1]]/4)
dds_non_shaking <- dds_non_shaking[keep_non_shaking,]

# Fit the model and run the DE analysis
dds_non_shaking = DESeq(dds_non_shaking)

vsd_non_shaking <- vst(dds_non_shaking, blind=FALSE)
rld_non_shaking <- rlog(dds_non_shaking, blind=FALSE)
head(assay(vsd_non_shaking), 3)

pcaData_non_shaking = plotPCA(vsd_non_shaking, intgroup=c("cell_line", "treatment"), returnData=TRUE)

percentVar_non_shaking <- round(100 * attr(pcaData_non_shaking, "percentVar"))
ggplot(pcaData_non_shaking, aes(PC1, PC2, color=cell_line, shape=treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar_non_shaking[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_non_shaking[2],"% variance")) + 
  coord_fixed() + scale_color_manual(values=c("A7" = "#082BD2", "A8" = "#65D527", "WT" = "#000000")) + ggtitle("2nd RNA-seq non-shaking PCA")

ggsave("nextflow_ru/run_2_nf/rna_2_non_shaking_pca.filtered.png")

save.image(file="/home/humebc/projects/ru/nextflow_ru/run_2_nf/run_2_nd.RData")



########## ADDITIONAL ANALYSES #########
# A strong result would be to compare the differential expressed genes between axenic WT and co-cultured WT
# (i.e. look for the effect of axenic state) to the DE for mutant vs wild type (both co-cultured).
# In the first comparison, the network would would be active in the WT co but not in the WT ax and thus,
# DE. In the second comparison, the network should be active in the WT co-culture, but nocked out in the
# mutant samples, so again, the network genes should be DE. Therefore I would expect an overlap of all 3 (first_de, WT_ax_vs_co_de, mut_vs_wt_co_de)
# of those groups because they should have the network genes as DE genes.
samples_ax_vs_co = samples %>% dplyr::filter(mutant=="FALSE" & shaking!="TRUE")

files_ax_vs_co <- file.path(kallisto.base.dir, samples_ax_vs_co$dir_name, "abundance.h5")

# Finally we can use tximport to read in the abundance tables
# and perform the normalizations
txi_ax_vs_co = tximport(files_ax_vs_co, type = "kallisto", tx2gene = tx2gene)

# Create the DESEQ2 object
dds_ax_vs_co = DESeqDataSetFromTximport(txi_ax_vs_co, colData = samples_ax_vs_co, design = ~ axenic)

# Filter out those genes with <10 counts in more than 1/4 of the samples
keep_ax_vs_co <- rowSums(counts(dds_ax_vs_co) >= 10) >= ceiling(dim(samples_ax_vs_co)[[1]]/4)
dds_ax_vs_co <- dds_ax_vs_co[keep_ax_vs_co,]

# Fit the model and run the DE analysis
# We will set the minRep argument to see if we can replace the outlier that is causing
# the NA value for the 43365 sample.
dds_ax_vs_co = DESeq(dds_ax_vs_co, minReplicatesForReplace=3)

# There is a single sample with a very high count that is causing the NA
# for the padj and pval. We will force replacement using the minReplicatesForReplace
# attribute of the DESeq method above.
# > assay(dds_ax_vs_co)["PHATRDRAFT_43365",]
#     WT_0h1     WT_0h2     WT_0h3   WT_24hB1   WT_24hB2   WT_24hB3 WT_24hBAX1 
#         97         94         46         53         52         47         45 
# WT_24hBAX2 WT_24hBAX3   WT_48hB1   WT_48hB2   WT_48hB3 WT_48hBAX1 WT_48hBAX2 
#         17         29       8509        163         71         81         69 
# WT_48hBAX3 
#         97

res_ax_vs_co = results(dds_ax_vs_co)

res_ax_vs_co = as.data.frame(res_ax_vs_co)

res_ax_vs_co.filtered = res_ax_vs_co %>% dplyr::filter(padj <= 0.01) %>% dplyr::arrange(padj)
# > dim(res_ax_vs_co.filtered)
# [1] 538   6
# There are many more DE expressed genes than for the 1st RNA seq comparison.

# The pval and padj have been set to NA for 43365.
# According to the docs this is something to do with an extreme count outlier
# See above. Now we have fixed this.
res_ax_vs_co["PHATRDRAFT_43365",]
# > res_ax_vs_co["PHATRDRAFT_43365",]
#                  baseMean log2FoldChange     lfcSE      stat   pvalue      padj
# PHATRDRAFT_43365 67.16041      0.2730666 0.3610843 0.7562406 0.449505 0.6615257
# The gene is no longer DE in the non-shaking samples!

# TODO have a look at the normalized gene counts and compare them to the RNA-seq 1 experiment.
# We should be able to see the same sorts of changes in the numbers

# we want to see if the other genes in the DE network that we created from the 
# 1st RNA-seq run are found in the DE genes of this study. This will be achieved
# by the Venn approach hopefully. See below.

# save this df for use later
save(res_ax_vs_co.filtered, file="/home/humebc/projects/ru/nextflow_ru/run_2_nf/run_2_res_ax_vs_co.filtered.RData")

########### Mutant DE ###############
# Get the mutant DE genes and make a comparison of the overlap to
# the set of genes above.
samples_mutant_vs_WT_co = samples %>% dplyr::filter(shaking!="TRUE" & axenic=="FALSE")

files_mutant_vs_WT_co <- file.path(kallisto.base.dir, samples_mutant_vs_WT_co$dir_name, "abundance.h5")

# Finally we can use tximport to read in the abundance tables
# and perform the normalizations
txi_mutant_vs_WT_co = tximport(files_mutant_vs_WT_co, type = "kallisto", tx2gene = tx2gene)

# Create the DESEQ2 object
dds_mutant_vs_WT_co = DESeqDataSetFromTximport(txi_mutant_vs_WT_co, colData = samples_mutant_vs_WT_co, design = ~ mutant)

# Filter out those genes with <10 counts in more than 1/4 of the samples
keep_mutant_vs_WT_co <- rowSums(counts(dds_mutant_vs_WT_co) >= 10) >= ceiling(dim(samples_mutant_vs_WT_co)[[1]]/4)
dds_mutant_vs_WT_co <- dds_mutant_vs_WT_co[keep_mutant_vs_WT_co,]

# Fit the model and run the DE analysis
dds_mutant_vs_WT_co = DESeq(dds_mutant_vs_WT_co, minReplicatesForReplace=3)

res_mutant_vs_WT_co = results(dds_mutant_vs_WT_co)

res_mutant_vs_WT_co = as.data.frame(res_mutant_vs_WT_co)

res_mutant_vs_WT_co.filtered = res_mutant_vs_WT_co %>% dplyr::filter(padj <= 0.01) %>% dplyr::arrange(padj)
# > dim(res_mutant_vs_WT_co.filtered)
# [1] 694   6
res_mutant_vs_WT_co.filtered["PHATRDRAFT_43365",]
# > res_mutant_vs_WT_co.filtered["PHATRDRAFT_43365",]
#    baseMean log2FoldChange lfcSE stat pvalue padj
# NA       NA             NA    NA   NA     NA   NA
res_mutant_vs_WT_co["PHATRDRAFT_43365",]
# > res_mutant_vs_WT_co["PHATRDRAFT_43365",]
#                  baseMean log2FoldChange     lfcSE     stat    pvalue      padj
# PHATRDRAFT_43365 159.3906       1.073634 0.6981413 1.537846 0.1240863 0.4490131

# The PHATRDRAFT_43365 is not differentially expressed between the mutants and the WT
# This may be believable, because both the mutants and the WT are co-cultured.
# However, we would expect the network that we identified in the 1st RNA-seq run
# to be knocked out. I.e. although the PHATRDRAFT_43365 transcript may be expressed
# its protein should be non effective and therefore the network should be knocked out.

assay(dds_mutant_vs_WT_co)["PHATRDRAFT_43365",]
# > assay(dds_mutant_vs_WT_co)["PHATRDRAFT_43365",]
# A7_24hB1 A7_24hB2 A7_24hB3 A7_48hB1 A7_48hB2 A7_48hB3 A8_24hB1 A8_24hB2 
#      166      165      227      737     1012      773       12       16 
# A8_24hB3 A8_48hB1 A8_48hB2 A8_48hB3 WT_24hB1 WT_24hB2 WT_24hB3 WT_48hB1 
#       20      126       60       28       53       52       47     8509 
# WT_48hB2 WT_48hB3 
#      163       71

# save this df for use later
save(res_mutant_vs_WT_co.filtered, file="/home/humebc/projects/ru/nextflow_ru/run_2_nf/run_2_res_mutant_vs_WT_co.filtered.RData")


# There is one more combination of samples to investigate and that is the
# DE between the mutant co-cultured and the axenic WT samples. For this comparison we would
# expect there to be a very minimal number of DE expressed genes.

samples_mutant_co_vs_WT_ax = samples %>% dplyr::filter(shaking!="TRUE" & !(cell_line=="WT" & axenic=="FALSE"))

files_mutant_co_vs_WT_ax <- file.path(kallisto.base.dir, samples_mutant_co_vs_WT_ax$dir_name, "abundance.h5")

# Finally we can use tximport to read in the abundance tables
# and perform the normalizations
txi_mutant_co_vs_WT_ax = tximport(files_mutant_co_vs_WT_ax, type = "kallisto", tx2gene = tx2gene)

# Create the DESEQ2 object
dds_mutant_co_vs_WT_ax = DESeqDataSetFromTximport(txi_mutant_co_vs_WT_ax, colData = samples_mutant_co_vs_WT_ax, design = ~ mutant)

# Filter out those genes with <10 counts in more than 1/4 of the samples
keep_mutant_co_vs_WT_ax <- rowSums(counts(dds_mutant_co_vs_WT_ax) >= 10) >= ceiling(dim(samples_mutant_co_vs_WT_ax)[[1]]/4)
dds_mutant_co_vs_WT_ax <- dds_mutant_co_vs_WT_ax[keep_mutant_co_vs_WT_ax,]

# Fit the model and run the DE analysis
dds_mutant_co_vs_WT_ax = DESeq(dds_mutant_co_vs_WT_ax, minReplicatesForReplace=3)

res_mutant_co_vs_WT_ax = results(dds_mutant_co_vs_WT_ax)

res_mutant_co_vs_WT_ax = as.data.frame(res_mutant_co_vs_WT_ax)

res_mutant_co_vs_WT_ax.filtered = res_mutant_co_vs_WT_ax %>% dplyr::filter(padj <= 0.01) %>% dplyr::arrange(padj)
# > dim(res_mutant_co_vs_WT_ax.filtered)
# [1] 1091    6
# This shows that there is a huge number of DE genes between the mutant co-cultured
# compared to the WT axenic which is not what we expected to see.

res_mutant_co_vs_WT_ax.filtered["PHATRDRAFT_43365",]
# > res_mutant_co_vs_WT_ax.filtered["PHATRDRAFT_43365",]
#    baseMean log2FoldChange lfcSE stat pvalue padj
# NA       NA             NA    NA   NA     NA   NA
res_mutant_co_vs_WT_ax["PHATRDRAFT_43365",]
# > res_mutant_co_vs_WT_ax["PHATRDRAFT_43365",]
#                  baseMean log2FoldChange     lfcSE     stat      pvalue       padj
# PHATRDRAFT_43365 157.4127       1.689483 0.6136731 2.753066 0.005903995 0.03726966

head(res_mutant_co_vs_WT_ax.filtered)
# > head(res_mutant_co_vs_WT_ax.filtered)
#                   baseMean log2FoldChange     lfcSE      stat        pvalue           padj
# PHATRDRAFT_51129 7811.1126     -11.614298 0.2200020 -52.79178  0.000000e+00   0.000000e+00
# PHATRDRAFT_37231 1535.1849       5.457324 0.1997412  27.32197 2.325922e-164  1.181220e-160
# PHATRDRAFT_38705  383.0974     -12.620949 0.4842504 -26.06286 9.619450e-150  3.256825e-146
# PHATRDRAFT_48286  270.0884     -12.116128 0.4953144 -24.46149 3.797690e-132  9.643285e-129
# PHATRDRAFT_46458  406.5408      -5.287903 0.2316197 -22.83010 2.304168e-115  4.680687e-112
# PHATRDRAFT_22332  229.1821     -11.075503 0.5082598 -21.79103 2.822602e-105  4.778195e-102

# There are HUGE differences between the mutants and the wild types

assay(dds_mutant_co_vs_WT_ax)["PHATRDRAFT_43365",]
# > assay(dds_mutant_co_vs_WT_ax)["PHATRDRAFT_43365",]
#     A7_0h1     A7_0h2     A7_0h3   A7_24hB1   A7_24hB2   A7_24hB3   A7_48hB1 
#        156        262        157        166        165        227        737 
#   A7_48hB2   A7_48hB3     A8_0h1     A8_0h2     A8_0h3   A8_24hB1   A8_24hB2 
#       1012        773         15         20          8         12         16 
#   A8_24hB3   A8_48hB1   A8_48hB2   A8_48hB3     WT_0h1     WT_0h2     WT_0h3 
#         20        126         60         28         97         94         46 
# WT_24hBAX1 WT_24hBAX2 WT_24hBAX3 WT_48hBAX1 WT_48hBAX2 WT_48hBAX3 
#         45         17         29         81         69         97

# save this df for use later
save(res_mutant_co_vs_WT_ax.filtered, file="/home/humebc/projects/ru/nextflow_ru/run_2_nf/run_2_res_mutant_co_vs_WT_ax.filtered.RData")

# I then want to do and overlap comparison with the 3 comparisons above and of the DE genes from the 1st RNA-seq.
# Load in the de_genes from the 1st RNA-seq run
load("/home/humebc/projects/ru/nextflow_ru/run_1_nf/de_genes.subset.RData")
head(de_genes)
de_l = list(
  first_de = rownames(de_genes),
  WT_ax_vs_co_de = rownames(res_ax_vs_co.filtered),
  mut_vs_wt_co_de = rownames(res_mutant_vs_WT_co.filtered),
  mut_co_vs_WT_ax = rownames(res_mutant_co_vs_WT_ax.filtered)
  )


ggvenn(
  de_l, 
  fill_color = c("#D721BB", "#144BD5", "#3CCA23", "#21CBCA"),
  stroke_size = 0.5, set_name_size = 4
  ) + theme_bw() + xlim(c(-2.5,3.25))
ggsave("/home/humebc/projects/ru/nextflow_ru/run_2_nf/de_genes.Venn.png")
# This is not a good result.
# I would expect that the DE genes of the WT_ax_vs_co_de_list and mut_vs_vt_co_de_list comparisons
# to have a large overlap
# That is, the axenic state difference (i.e. WT co vs WT ax) would be the same as the
# mutant vs WT (co-cultured) difference. But there is only a very small over lap.
# I would also expect there to be an almost complete overlap of the DE genes of 
# the first analysis with those of the second ax vs co de list. But there is 0 overlap! 
# This second comparison is essentially a control for the experiment.


# I don't see any point in pursuing the WGCNA analysis as there is no way that the network
# is maintained as the fist_de_list vs WT_ax_vs_co_de_list has 0 overlap.
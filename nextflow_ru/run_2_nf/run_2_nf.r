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

    # Fit the model and run the DE analysis
    dds = DESeq(dds)
    
    res = results(dds)
    
    up = as.data.frame(res) %>% dplyr::filter(log2FoldChange > 2 & padj < 0.05)
    contrast_time = append(contrast_time, time); num_genes = append(num_genes, dim(up)[[1]]); up_down = append(up_down, "up");
    
    down = as.data.frame(res) %>% dplyr::filter(log2FoldChange < -2 & padj < 0.05)
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
ggsave("nextflow_ru/run_2_nf/rna_2_shaking_stacked_bars.png")



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

ggsave("nextflow_ru/run_2_nf/rna_2_shaking_pca.png")
# This PCA looks quite different to that of the one that Ru produced.


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
  dds = DESeq(dds)
  res = results(dds)

  up = as.data.frame(res) %>% dplyr::filter(log2FoldChange > 2 & padj < 0.05)
  contrast_non_shaking = append(contrast_non_shaking, paste0(time, "h_mutants_vs_WT")); num_genes_non_shaking = append(num_genes_non_shaking, dim(up)[[1]]); up_down_non_shaking = append(up_down_non_shaking, "up");

  down = as.data.frame(res) %>% dplyr::filter(log2FoldChange < -2 & padj < 0.05)
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

ggsave("nextflow_ru/run_2_nf/rna_2_non_shaking_stacked_bars.png")

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

ggsave("nextflow_ru/run_2_nf/rna_2_non_shaking_pca.png")

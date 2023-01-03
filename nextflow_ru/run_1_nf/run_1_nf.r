# This is the processing of Ru's RNA-seq data
# See run_1.nf for the preprocessing of the read data. I simply ran it through fastp and then kallisto.
# The transcriptome was downloaded form here: https://www.ncbi.nlm.nih.gov/genome/?term=txid556484[Organism:noexp]
# I have created a tx2gene file using: gffread -E GCF_000150955.2_ASM15095v2_genomic.gff --table @id,@geneid > tx2gene.raw
# That I got from here: https://www.biostars.org/p/9513656/#9513659
# I then further processed this file with: awk -F "\t" 'BEGIN{print "TXNAME\tGENEID"} {gsub("rna-", "", $1); gsub("gene-", "", $2); print $1 "\t" $2}' tx2gene.raw > tx2gene.txt
# To produce the final tx2gene file

# I also created a samples meta file: /home/humebc/projects/ru/nextflow_ru/run_1_nf/samples_run_1.csv

library(dplyr)
library(stringr)
library(tximport)
library(DESeq2)
library("pheatmap")
library(ggplot2)
library(ggvenn)

# If the script has already been run through then you can use this load.
load(file = "nextflow_ru/run_1_nf/run_1_data.RData")

tx2gene = read.table("/home/humebc/projects/ru/reference/tx2gene.txt", header=TRUE)
samples = read.csv("/home/humebc/projects/ru/nextflow_ru/run_1_nf/samples_run_1.csv", header=TRUE)
samples = samples %>% mutate(axenic = as.factor(axenic), time_hr = as.factor(time_hr)) %>% mutate(axenic = relevel(axenic, "TRUE"))


# Make a vector that contains the full paths to the abundance.h5 files
kallisto.base.dir = "/home/humebc/projects/ru/nextflow_ru/run_1_nf/results/kallisto"


# There are three initial results that we want to recapitulate

### AXENIC VS CO-CULTURED STACKED BAR ###
# The first is the number of differential genes realised for each of the time points
# comparing control to axenic

# We will want to run a DESEQ2 with a subset of the samples to get the DE for each of the time points
times = c(0.5, 3, 24, 48)

# To collect the results of how many genes are up and down regulated for
# each of the time points we will make some empty vectors
contrast_time = c()
num_genes = c()
up_down = c()
# We will also want to collect the list of genes that are
# DE as we want to find those genes that are DE in common
# across each of the time comparisons
de_genes = c()

# The loop for doing each of the subsets
for (time in times){
    samples_sub = samples %>% dplyr::filter(time_hr==time)

    files <- file.path(kallisto.base.dir, samples_sub$dir_name, "abundance.h5")

    # Finally we can use tximport to read in the abundance tables
    # and perform the normalizations
    txi = tximport(files, type = "kallisto", tx2gene = tx2gene)

    # Create the DESEQ2 object
    dds = DESeqDataSetFromTximport(txi, colData = samples_sub, design = ~ axenic)

    # Filter out those genes with <10 counts in more than 1/4 of the samples
    keep <- rowSums(counts(dds) >= 10) >= ceiling(dim(samples_sub)[[1]]/4)
    dds <- dds[keep,]

    # Fit the model and run the DE analysis
    dds = DESeq(dds)
    
    res = results(dds)
    
    up = as.data.frame(res) %>% dplyr::filter(log2FoldChange > 2 & padj < 0.05)
    contrast_time = append(contrast_time, time); num_genes = append(num_genes, dim(up)[[1]]); up_down = append(up_down, "up");
    
    down = as.data.frame(res) %>% dplyr::filter(log2FoldChange < -2 & padj < 0.05)
    contrast_time = append(contrast_time, time); num_genes = append(num_genes, dim(down)[[1]]); up_down = append(up_down, "down");
    
    # Collect the differentially expressed genes
    de_genes[[as.character(time)]] = c(rownames(up), rownames(down))
}


# Create the df that we will use for plotting
plotting_df = data.frame(contrast_time=as.factor(contrast_time), num_genes=num_genes, up_down=up_down)

# Plot up the results in the same way as Ru.
# I.e. stacked bar plot, one stack for each time with up and down regulated
# gene counts for each bar.
ggplot(plotting_df, aes(fill=up_down, y=num_genes, x=contrast_time, label=num_genes)) + 
    geom_bar(position="stack", stat="identity") + geom_text(size = 10, position = position_stack(vjust = 0.5)) + 
    scale_fill_manual(values=c("up" = "#404040", "down" = "#AFABAB")) + ggtitle("1st RNA-seq run DEGs")

ggsave("nextflow_ru/run_1_nf/rna_1_stacked_bars.filtered.png")



### PCA ###
# Next we want to produce a PCA for all of the samples. This means we'll
# want to make a DESeq2 object that contains all of the samples
files_all_samples <- file.path(kallisto.base.dir, samples$dir_name, "abundance.h5")

# Finally we can use tximport to read in the abundance tables
# and perform the normalizations
txi_all_samples = tximport(files_all_samples, type = "kallisto", tx2gene = tx2gene)

# Create the DESEQ2 object
dds_all_samples = DESeqDataSetFromTximport(txi_all_samples, colData = samples, design = ~ axenic)

# Filter out those genes with <10 counts in more than 1/4 of the samples
keep_all_samples <- rowSums(counts(dds_all_samples) >= 10) >= ceiling(dim(samples)[[1]]/4)
dds_all_samples <- dds_all_samples[keep_all_samples,]

# Fit the model and run the DE analysis
dds_all_samples = DESeq(dds_all_samples)

vsd_all_samples <- vst(dds_all_samples, blind=FALSE)
rld_all_samples <- rlog(dds_all_samples, blind=FALSE)
head(assay(vsd_all_samples), 3)

pcaData = plotPCA(vsd_all_samples, intgroup=c("time_hr", "axenic"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=time_hr, shape=axenic)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + scale_color_manual(values=c("0" = "#000000", "0.5" = "#D721BB", "3" = "#144BD5", "24" = "#3CCA23", "48" = "#21CBCA")) +
  ggtitle("1st RNA-seq PCA all samples")

ggsave("nextflow_ru/run_1_nf/rna_1_all_sample_pca.filtered.png")

# The PCA seems to be in good agreement with Ru's PCA

# Next we want to find the genes that are DE for each
# of the time comparisons
# Turns out that there are none. In particular we will look up the results for the PHATRDRAFT_43365 gene.
common_genes = Reduce(intersect, de_genes) # Empty.
# PHATRDRAFT_43365 is not a DEG in the 48 time point according to our analysis
#                   baseMean log2FoldChange     lfcSE      stat    pvalue      padj
#                  <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
# PHATRDRAFT_43365   582.7921       1.417574 0.8136896 1.742155 0.08148128   0.2282742

# Make a Venn of the gene overlap
ggvenn(
  de_genes, 
  fill_color = c("#D721BB", "#144BD5", "#3CCA23", "#21CBCA"),
  stroke_size = 0.5, set_name_size = 8
  )
ggsave("nextflow_ru/run_1_nf/rna_1_gene_venn.filtered.png")

# Let's make some heat maps to look at the gene expressions across the samples
annotation_df_all_samples <- as.data.frame(colData(dds_all_samples)[,c("time_hr","axenic")])
rownames(annotation_df_all_samples) = colData(dds_all_samples)$dir_name
# Assign the assay object to a variable so that we can name the rows of the heat map
# (columns of the vsd_assay df)
vsd_assay_all_samples = assay(vsd_all_samples)
colnames(vsd_assay_all_samples) = rownames(annotation_df_all_samples)

# Get a unique list of the DE genes across all comparisons
de_genes_unique = unique(Reduce(c, de_genes))

# First the heat map of the 43365 gene across all samples
# We make the heat map with one other random gene as the heat map doesn't work with 1 row.
heat = pheatmap(vsd_assay_all_samples[c("PHATRDRAFT_43365", "PHATRDRAFT_55137"),], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=TRUE, scale="row", annotation_col = annotation_df_all_samples)
ggsave("nextflow_ru/run_1_nf/rna_1_heatmap_PHATRDRAFT_43365.png", heat)

# Then a heat map of all DE genes across all samples
heat = pheatmap(vsd_assay_all_samples[de_genes_unique,], cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, scale="row", annotation_col = annotation_df_all_samples)
ggsave("nextflow_ru/run_1_nf/rna_1_heatmap_all_DE.png", heat)

# And finally do one with only the 48hr samples
time_48 = 48
samples_sub_48 = samples %>% dplyr::filter(time_hr==time_48)

files_only_48 <- file.path(kallisto.base.dir, samples_sub_48$dir_name, "abundance.h5")

# Finally we can use tximport to read in the abundance tables
# and perform the normalizations
txi_only_48 = tximport(files_only_48, type = "kallisto", tx2gene = tx2gene)

# Create the DESEQ2 object
dds_only_48 = DESeqDataSetFromTximport(txi_only_48, colData = samples_sub_48, design = ~ axenic)

# Fit the model and run the DE analysis
dds_only_48 = DESeq(dds_only_48)
vsd_only_48 <- vst(dds_only_48, blind=FALSE)

annotation_df_only_48 <- as.data.frame(colData(dds_only_48)[,c("time_hr","axenic")])
rownames(df) = colData(dds_only_48)$dir_name
vsd_assay_only_48 = assay(vsd_only_48)
colnames(vsd_assay_only_48) = rownames(annotation_df_only_48)

# We make the heat map with one other random gene as the heat map doesn't work with 1 row.
heat = pheatmap(vsd_assay_only_48[c("PHATRDRAFT_43365", "PHATRDRAFT_55137"),], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=TRUE, scale="row", annotation_col = annotation_df_only_48)
ggsave("nextflow_ru/run_1_nf/rna_1_heatmap_PHATRDRAFT_43365_only_48.png", heat)

# Finally as an additional approach to this analysis I want to run a DE expression
# using some subsets of the data
# We have seen from the PCA that the difference between the axenic states is
# largest for the 24 hour samples but I would also include the 3 and 48 hour samples
# and see if we can find significant genes.
# Let's start by running the DE with all of the time points except 0.
time_points = c(3, 24, 48)
tx2gene = read.table("/home/humebc/projects/ru/reference/tx2gene.txt", header=TRUE)
samples = read.csv("/home/humebc/projects/ru/nextflow_ru/run_1_nf/samples_run_1.csv", header=TRUE)
samples = samples %>% mutate(axenic = as.factor(axenic), time_hr = as.factor(time_hr)) %>% mutate(axenic = relevel(axenic, "TRUE"))
samples_sub = samples %>% dplyr::filter(time_hr %in% time_points)

# Make a vector that contains the full paths to the abundance.h5 files
kallisto.base.dir = "/home/humebc/projects/ru/nextflow_ru/run_1_nf/results/kallisto"

files <- file.path(kallisto.base.dir, samples_sub$dir_name, "abundance.h5")

# Finally we can use tximport to read in the abundance tables
# and perform the normalizations
txi = tximport(files, type = "kallisto", tx2gene = tx2gene)

# Create the DESEQ2 object
dds = DESeqDataSetFromTximport(txi, colData = samples_sub, design = ~ axenic)

# Filter out those genes with <10 counts in more than 1/4 of the samples
keep <- rowSums(counts(dds) >= 10) >= ceiling(dim(samples_sub)[[1]]/4)
dds <- dds[keep,]

# Fit the model and run the DE analysis
dds = DESeq(dds)

res = results(dds)

res.df = as.data.frame(res) %>% dplyr::arrange(padj)
head(res.df, 20)
# > head(res.df)
#                    baseMean log2FoldChange      lfcSE      stat       pvalue
# PHATRDRAFT_43365  1548.7911      3.8539681 0.37605531 10.248408 1.203091e-24
# PHATRDRAFT_48554  4218.4944      4.5133006 0.51335526  8.791768 1.472240e-18
# PHATRDRAFT_50288  7614.6849      1.1720971 0.16612762  7.055402 1.721023e-12
# PHATRDRAFT_44925   187.1997      4.3949089 0.72402836  6.070078 1.278478e-09
# PHATRDRAFT_43020 13960.1351      0.4619692 0.08115779  5.692235 1.253868e-08
# PHATRDRAFT_48069  3156.8690      1.1386302 0.19977389  5.699595 1.200926e-08
#                          padj
# PHATRDRAFT_43365 1.218370e-20
# PHATRDRAFT_48554 7.454685e-15
# PHATRDRAFT_50288 5.809599e-09
# PHATRDRAFT_44925 3.236787e-06
# PHATRDRAFT_43020 2.116320e-05
# PHATRDRAFT_48069 2.116320e-05

# Write out the list of DE genes P<0.01
de_genes = res.df %>% dplyr::filter(padj<0.01) %>% arrange(padj)
save(de_genes, file="/home/humebc/projects/ru/nextflow_ru/run_1_nf/de_genes.subset.RData")

# Now compare this to the adjacency score of 43365
load("/home/humebc/projects/ru/nextflow_ru/run_1_nf/adj_plot_df.RData")
head(adj_plot_df)
# I want to plot the significance of the DE expressed genes and highlight those with high adjacency in red
high_adj_genes = row.names(adj_plot_df %>% dplyr::filter(adj > 0.1))
res.df.plot = res.df %>% mutate(x_val=1, rn=row.names(res.df), col=dplyr::case_when(rn %in% high_adj_genes ~ "red", TRUE ~ "black")) %>% dplyr::filter(padj<0.01)
ggplot(res.df.plot, aes(x=x_val, y=-log10(padj), lab=rn, color=col)) + geom_point() + geom_label_repel(aes(label=rn), label.size=0.01) + scale_color_manual(values = c("red" = "red", "black" = "black")) + guides(color = "none") +
ggtitle("DE genes (P<0.01); red == adjacency to PHATRDRAFT_43365 > 0.1")
ggsave("/home/humebc/projects/ru/nextflow_ru/run_1_nf/DEgenes.adjacency.subset.png", height=40, width=20, units="cm")
# This confirms the findings that the 43365 gene is the right gene to knock out.

# SUMMARY: The overall conclusions of the R work are that the gene they knocked out is
# was a good choice as it shows up in the WGCNA and DESEQ2 analyses. However,
# we were not able to identify gene modules that relate to the axenic state using
# the three time points.

# Save the data
save.image(file = "nextflow_ru/run_1_nf/run_1_data.RData")

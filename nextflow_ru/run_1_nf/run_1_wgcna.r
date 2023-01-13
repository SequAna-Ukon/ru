# This is the WGCNA-based analysis of Ru's 1st RNA-seq dataset.
# We have already run a prelimiary analysis of the 1st RNA-seq dataset
# in the run_1_nf.r R script.


# In this script we will:
#   Identify modules
#   Compute significance with the axenic/co-cultured trait
#   Inspect the genes of anysignificant modules.

# In particular we are interested in seeing where the PHATRDRAFT_43365 gene
# that Ru identified as DE between the different time points lies in the networks.

# We can either create two separate networks, 1 for anexic 1 for co-cultured
# or we can run trait correlations for each of the modules.


# We will work with tximport and DESEQ2 to generate the normalized abundances
# for use with the WGCNA analysis.
# Here is a resource about filtering the counts:
# https://support.bioconductor.org/p/115583/
# It essentially suggests filtering out very low level transcripts
# but we will have to check to see if the gene of interest has sufficient
# abundance to be retained.

#==============================================
#
# Read in the abundance tables for all samples
# using tximport and normalize with DESEQ2
# Filter according to abundance
#
#==============================================

library(dplyr)
library(stringr)
library(tximport)
library(DESeq2)
library("pheatmap")
library(ggplot2)
library(ggvenn)
library(WGCNA)
library(gridExtra)
library(ggdendro)
library(ggrepel)

# variable to choose which transcriptome to work with
# values can be either "NCBI" or "ensembl"
transcriptome = "ensembl"

# Read in the tx2gene mapping file
if (transcriptome == "ensembl") {
  tx2gene = read.table("/home/humebc/projects/ru/reference/alt_transcriptome/tx2gene.txt", header=TRUE)
}
if (transcriptome == "NCBI") {
  tx2gene = read.table("/home/humebc/projects/ru/reference/tx2gene.txt", header=TRUE)
}
# Read in the samples meta info
samples = read.csv("/home/humebc/projects/ru/nextflow_ru/run_1_nf/samples_run_1.csv", header=TRUE)
samples = samples %>% mutate(axenic = as.factor(axenic), time_hr = as.factor(time_hr)) %>% mutate(axenic = relevel(axenic, "TRUE"))

rownames(samples) = lapply(lapply(str_split(samples$dir_name, "_"), FUN=head, -4), str_c, collapse="_")

# Make a vector that contains the full paths to the abundance.h5 files
if (transcriptome == "ensembl") {
  kallisto.base.dir = "/home/humebc/projects/ru/nextflow_ru/run_1_nf/results_alt_ref/kallisto"
}
if (transcriptome == "NCBI") {
  kallisto.base.dir = "/home/humebc/projects/ru/nextflow_ru/run_1_nf/results/kallisto"
}

files <- file.path(kallisto.base.dir, samples$dir_name, "abundance.h5")
txi = tximport(files, type = "kallisto", tx2gene = tx2gene)

# Create the DESEQ2 object
dds = DESeqDataSetFromTximport(txi, colData = samples, design = ~ axenic)

# Filter out those genes with <10 counts in more than 1/4 of the samples
keep <- rowSums(counts(dds) >= 10) >= ceiling(dim(samples)[[1]]/4)
dds <- dds[keep,]

# Fit the model and run the DE analysis
dds = DESeq(dds)

# Now transform and normalize the data
# This will return log2 transformed data normalized by library size
# and with the experiment-wide mean-variance trend removed.
vsd <- vst(dds, blind=FALSE)

summary(as.data.frame(assay(vsd)))

wgcna_in = t(as.data.frame(assay(vsd)))

#==============================================
#
# Begin the WGCNA analysis
# Visualize the sample similarity
#
#==============================================
# Look for outliers
hc = hclust(dist(wgcna_in), method = "average")
ggdendrogram(hc, rotate = FALSE, size = 2)
# All looks good but we can once again see that the 30_min_1 sample looks pretty weird (NG-23864_C_30min_1)
# Strikes me that this might be the same sample that was strange at 48 hours?

# Check that all of the samples are good with regards to missing data
gsg = goodSamplesGenes(wgcna_in, verbose = 3)
gsg$allOK # All OK

# Create a traits file from the samples file that contains the factors we are interested in (axenic, time_hr)
# converted to numeric
traits = samples %>% select(axenic, time_hr) %>% mutate(axenic = as.numeric(axenic), time_hr = as.numeric(time_hr))

traitColors = numbers2colors(traits, signed = FALSE)

plotDendroAndColors(hc, traitColors,
                    groupLabels = names(traits), 
                    main = "Sample dendrogram and trait heatmap")

#==============================================
#
# Construct WGCNA network
#
#==============================================

###############################################
# CHOOSING SOFT THRESHOLD POWER
###############################################
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(wgcna_in, powerVector = powers, verbose = 5)

# Plot scale-free toplogy model fit R^2 value against power
text_size = 5
sft_plotting_df = sft$fitIndices
r_sqr_plot = ggplot(sft_plotting_df, aes(x=Power, y=-sign(slope)*SFT.R.sq, label=Power)) + 
geom_point(size=0.5) + geom_text(hjust=0, vjust=0, size=3, color=ifelse(sft_plotting_df$Power==16, 'red', 'black')) + ggtitle("Scale independence") + 
    theme(plot.title = element_text(size = text_size, face = "bold"), 
    axis.title=element_text(size=text_size,face="bold"), 
    axis.text=element_text(size=text_size)) + xlab("Soft Threshold (Power)") + geom_hline(yintercept=0.8, color="red")


mean_conn_plot = ggplot(sft_plotting_df, aes(x=Power, y=mean.k., label=Power)) + 
geom_point(size=0.5) + geom_text(hjust=0, vjust=0, size=3, color=ifelse(sft_plotting_df$Power==16, 'red', 'black')) + ggtitle("Mean connectivity") + 
    theme(plot.title = element_text(size = text_size, face = "bold"),
    axis.title=element_text(size=text_size,face="bold"), axis.text=element_text(size=text_size)) + 
    xlab("Soft Threshold (Power)")


threshold_plot = grid.arrange(r_sqr_plot, mean_conn_plot, ncol=2)
# Make a vector that contains the full paths to the abundance.h5 files
if (transcriptome == "ensembl") {
  ggsave("nextflow_ru/run_1_nf/run_1_wgcna_threshold_selection.ensembl.png", plot=threshold_plot, width=20, height=10, units="cm")
}
if (transcriptome == "NCBI") {
  ggsave("nextflow_ru/run_1_nf/run_1_wgcna_threshold_selection.png", plot=threshold_plot, width=20, height=10, units="cm")
}

# We will work with the soft threshold power value of 16

###############################################
# MAKING THE NETWORK
###############################################

# Let's try enabling the  multithreading
enableWGCNAThreads(nThreads = 40) # As far as I can tell this doesn't do much

if (transcriptome == "ensembl") {
    if(file.exists("/home/humebc/projects/ru/nextflow_ru/run_1_nf/adjacency_TOM.ensembl.RData")){
        load("/home/humebc/projects/ru/nextflow_ru/run_1_nf/adjacency_TOM.ensembl.RData")
    }else{
        # Make adjacency matrix
        adjacency = adjacency(wgcna_in, power = 16)
        # Turn adjacency into topological overlap
        TOM = TOMsimilarity(adjacency)
        save(adjacency, TOM, file="/home/humebc/projects/ru/nextflow_ru/run_1_nf/adjacency_TOM.ensembl.RData")
    }
}

if (transcriptome == "NCBI") {
    if(file.exists("/home/humebc/projects/ru/nextflow_ru/run_1_nf/adjacency_TOM.RData")){
        load("/home/humebc/projects/ru/nextflow_ru/run_1_nf/adjacency_TOM.RData")
    }else{
        # Make adjacency matrix
        adjacency = adjacency(wgcna_in, power = 16)
        # Turn adjacency into topological overlap
        TOM = TOMsimilarity(adjacency)
        save(adjacency, TOM, file="/home/humebc/projects/ru/nextflow_ru/run_1_nf/adjacency_TOM.RData")
    }
}


dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
table(dynamicMods)
# 0 is automatically assigned to grey in the labels2colors method.
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(wgcna_in, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result

plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")

MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(wgcna_in, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(200));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts

if (transcriptome == "ensembl") {
    save(MEs, moduleLabels, moduleColors, geneTree, file = "/home/humebc/projects/ru/nextflow_ru/run_1_nf/network_objs.ensembl.RData")
}
if (transcriptome == "NCBI") {
    save(MEs, moduleLabels, moduleColors, geneTree, file = "/home/humebc/projects/ru/nextflow_ru/run_1_nf/network_objs.RData")
}


#==============================================
#
# Module-trait analysis
#
#==============================================

# At this point I'm thinking it's probably best if we exclude the time
# point 0 samples as the first comparison we're interested in is the
# comparison of axenic against co-cultured. The time 0 samples
# will not be informative to this comparison.

# For the time being we will continue with all of the samples
# Define numbers of genes and samples
nGenes = ncol(wgcna_in);
nSamples = nrow(wgcna_in);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(wgcna_in, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# This is quite exciting there are some strong hits in both the time and the axenic
# when working with all of the samples:
# > moduleTraitCor
#                        axenic     time_hr
# MEbrown           -0.13485185 -0.82363556
# MEdarkorange2     -0.23719593  0.09242001
# MEwhite            0.01461977  0.23831976
# MEblue2           -0.12012400  0.45958048
# MEdarkolivegreen  -0.07279308  0.82099905
# MEantiquewhite4   -0.32174810 -0.14973524
# MElavenderblush3  -0.12060523  0.24447477
# MEmaroon          -0.03588692  0.10905490
# MEskyblue3        -0.18347382 -0.47190625
# MEfirebrick4       0.29725367  0.42234363
# MEhoneydew1        0.29794686  0.81547830
# MEcoral1           0.15595787  0.33994439
# MEdarkolivegreen4 -0.12416275 -0.03507466
# MEgrey            -0.78242792  0.26835019
# > moduleTraitPvalue
#                         axenic      time_hr
# MEbrown           5.024578e-01 1.310439e-07
# MEdarkorange2     2.335439e-01 6.466018e-01
# MEwhite           9.423030e-01 2.312707e-01
# MEblue2           5.506311e-01 1.587711e-02
# MEdarkolivegreen  7.182301e-01 1.553702e-07
# MEantiquewhite4   1.017149e-01 4.560003e-01
# MElavenderblush3  5.490240e-01 2.190824e-01
# MEmaroon          8.589530e-01 5.881827e-01
# MEskyblue3        3.596374e-01 1.294918e-02
# MEfirebrick4      1.321283e-01 2.819464e-02
# MEhoneydew1       1.311859e-01 2.199968e-07
# MEcoral1          4.372732e-01 8.275592e-02
# MEdarkolivegreen4 5.372115e-01 8.621151e-01
# MEgrey            1.422931e-06 1.759325e-01


textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)

labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(traits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))

#==============================================
#
# Gene Module Membership
#
#==============================================

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(wgcna_in, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

#==============================================
#
# Gene Trait Significance
#
#==============================================
# Define variable axenic containing the axenic column of the trait df
axenic = as.data.frame(traits$axenic)
names(axenic) = "axenic"
geneTraitSignificance_axenic = as.data.frame(cor(wgcna_in, axenic, use = "p"))
GSPvalue_axenic = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_axenic), nSamples))
names(geneTraitSignificance_axenic) = paste("GS.", names(axenic), sep="")
names(GSPvalue_axenic) = paste("p.GS.", names(axenic), sep="")

time_hr = as.data.frame(traits$time_hr)
names(time_hr) = "time_hr"
geneTraitSignificance_time_hr = as.data.frame(cor(wgcna_in, time_hr, use = "p"))
GSPvalue_time_hr = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_time_hr), nSamples))
names(geneTraitSignificance_time_hr) = paste("GS.", names(time_hr), sep="")
names(GSPvalue_time_hr) = paste("p.GS.", names(time_hr), sep="")

# Put all four of the dataframes together into a single df
l = list(geneTraitSignificance_axenic, GSPvalue_axenic, geneTraitSignificance_time_hr, GSPvalue_time_hr)
geneTraitSignificance = Reduce(merge, lapply(l, function(x) data.frame(x, rn = row.names(x))))
row.names(geneTraitSignificance) = geneTraitSignificance$rn
geneTraitSignificance = geneTraitSignificance %>% select(-rn)

# Search for the top 10 genes that relate to axenic state
head(geneTraitSignificance %>% arrange(desc(abs(GS.axenic))),10)
# Specifically look at the PHATRDRAFT_43365 gene
if (transcriptome == "ensembl") {
    geneTraitSignificance["Phatr3_J43365",]
}
if (transcriptome == "NCBI") {
    geneTraitSignificance["PHATRDRAFT_43365",]
}

# > geneTraitSignificance["PHATRDRAFT_43365",]
#                  GS.axenic p.GS.axenic GS.time_hr p.GS.time_hr
# PHATRDRAFT_43365 0.5717391 0.001836068 -0.4218402   0.02840287

# > geneTraitSignificance["Phatr3_J43365",]
#               GS.axenic p.GS.axenic GS.time_hr p.GS.time_hr
# Phatr3_J43365 0.5699869 0.001910303 -0.4253445   0.02697886


# Find the genes with the highest connectivity (i.e. hub gene of the module)
# we will move on at this point as the grey module is the 'outlier' module

# I think it would be a good idea to rerun the above code but for a limited set of samples.
# Looking at the rna_1_all_sample_pca.filtered.png PCA it seems as though there is little
# difference at 0 and 0.5 hours between the axenic states. By contrast the biggest difference
# appears to come at 24 hours. So I thought it would be good to work with the 3, 24 and 48 samples

# I will run this code in a separate R script run_1_wgcna.follow.up.1.r
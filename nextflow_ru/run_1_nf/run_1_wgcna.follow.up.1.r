# This is the WGCNA-based analysis of Ru's 1st RNA-seq dataset.
# We have already run prelimiary analyses of the 1st RNA-seq dataset
# in the run_1_nf.r and run_1_wgcna.r R scripts.

# In this script we will redo much of what we did in run_1_wgcna.r but for a subset of the samples
# We will refactor much of the code used in run_1_wgcna.r for this purpose.



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

sample_subsets = list(c(3, 24, 48), c(3, 24), c(24))
sample_subset = c(3, 24, 48)
# Read in the tx2gene mapping file
tx2gene = read.table("/home/humebc/projects/ru/reference/tx2gene.txt", header=TRUE)
# Read in the samples meta info
samples = read.csv("/home/humebc/projects/ru/nextflow_ru/run_1_nf/samples_run_1.csv", header=TRUE)
samples = samples %>% mutate(axenic = as.factor(axenic), time_hr = as.factor(time_hr)) %>% mutate(axenic = relevel(axenic, "TRUE"))

# subset the samples to the relevant subset
samples = samples %>% dplyr::filter(time_hr %in% sample_subset)

rownames(samples) = lapply(lapply(str_split(samples$dir_name, "_"), FUN=head, -4), str_c, collapse="_")

# Make a vector that contains the full paths to the abundance.h5 files
kallisto.base.dir = "/home/humebc/projects/ru/nextflow_ru/run_1_nf/results/kallisto"

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

# Worryingly the R^2 value plateaus a little over 0.75.

# Call the network topology analysis function
sft = pickSoftThreshold(wgcna_in, powerVector = powers, verbose = 5, blockSize=20000)

# Plot scale-free toplogy model fit R^2 value against power
text_size = 5
sft_plotting_df = sft$fitIndices
r_sqr_plot = ggplot(sft_plotting_df, aes(x=Power, y=-sign(slope)*SFT.R.sq, label=Power)) + 
geom_point(size=0.5) + geom_text(hjust=0, vjust=0, size=3, color=ifelse(sft_plotting_df$Power==12, 'red', 'black')) + ggtitle("Scale independence") + 
    theme(plot.title = element_text(size = text_size, face = "bold"), 
    axis.title=element_text(size=text_size,face="bold"), 
    axis.text=element_text(size=text_size)) + xlab("Soft Threshold (Power)") + geom_hline(yintercept=0.8, color="red")


mean_conn_plot = ggplot(sft_plotting_df, aes(x=Power, y=mean.k., label=Power)) + 
geom_point(size=0.5) + geom_text(hjust=0, vjust=0, size=3, color=ifelse(sft_plotting_df$Power==12, 'red', 'black')) + ggtitle("Mean connectivity") + 
    theme(plot.title = element_text(size = text_size, face = "bold"),
    axis.title=element_text(size=text_size,face="bold"), axis.text=element_text(size=text_size)) + 
    xlab("Soft Threshold (Power)")


threshold_plot = grid.arrange(r_sqr_plot, mean_conn_plot, ncol=2)
ggsave("nextflow_ru/run_1_nf/run_1_wgcna_threshold_selection.subset.png", plot=threshold_plot, width=20, height=10, units="cm")

# We will work with the soft threshold power value of 12

###############################################
# MAKING THE NETWORK
###############################################

# Let's try enabling the  multithreading
enableWGCNAThreads(nThreads = 40) # As far as I can tell this doesn't do much


# Make adjacency matrix
adjacency = adjacency(wgcna_in, power = 12)
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency)


dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
ggdendrogram(geneTree, rotate = FALSE, size = 2)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
table(dynamicMods)
# 0 is automatically assigned to grey by default
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")

# Run the module-trait analysis on the non-merged modules
nGenes = ncol(wgcna_in);
nSamples = nrow(wgcna_in);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(wgcna_in, colors = dynamicColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# Similar to the merged module-traits that I have calculated below
# we see that there are non modules that really correlate strongly with
# the axenic state.
# > moduleTraitCor
#                         axenic     time_hr
# MEcoral2           0.216650421 -0.72191909
# MEdarkgreen       -0.021861562 -0.90754354
# MEgreen           -0.006842228 -0.96180656
# MEorangered4      -0.025494918 -0.78528799
# MEcyan            -0.025047136 -0.93602960
# MEred             -0.024039038 -0.85591912
# MEbrown            0.143712910 -0.94205619
# MEturquoise        0.029696196 -0.98585824
# MEivory           -0.033012474 -0.86351007
# MEthistle1         0.227278889 -0.90772168
# MEblack           -0.076003588 -0.29415752
# MEtan             -0.079598590 -0.02151416
# MEgrey60          -0.148094567 -0.41615847
# MEyellow          -0.030021702 -0.58097130
# MEdarkorange      -0.168799370 -0.76525688
# MEgreenyellow     -0.180893143 -0.61657579
# MEplum1           -0.125072704 -0.75207553
# MEfloralwhite      0.154064457 -0.49907431
# MEpurple           0.168753335 -0.80676986
# MElavenderblush3  -0.190589457  0.68939051
# MEcoral1           0.016975446  0.43181526
# MEdarkseagreen4   -0.054431158  0.19874878
# MEdarkolivegreen   0.126474889  0.12910571
# MEyellowgreen      0.208951377  0.14214966
# MEviolet           0.184053081  0.11493164
# MEbisque4          0.201826451 -0.20432023
# MElightpink4       0.137960394 -0.13690662
# MElightsteelblue1  0.070472638  0.01059375
# MEpaleturquoise    0.060359394  0.65653329
# MEsteelblue        0.032552061  0.56972727
# MEthistle2         0.052555367  0.73603238
# MEdarkturquoise    0.117643098  0.34280859
# MElightcyan        0.160425266  0.34863605
# MEroyalblue        0.127760499  0.55771977
# MEsienna3         -0.322895516  0.79283885
# MElightgreen      -0.133628241  0.76847195
# MEorange          -0.046040725  0.83990182
# MEpink            -0.020953792  0.94766169
# MEmagenta         -0.080067808  0.24841271
# MEdarkred         -0.112035212  0.46230637
# MEsalmon          -0.112143517  0.64281569
# MEbrown4          -0.101515163  0.33930587
# MEdarkorange2     -0.198741665  0.17266662
# MElightcyan1      -0.217403436  0.43136774
# MEskyblue3        -0.156045769  0.60844191
# MEhoneydew1       -0.080789987  0.72774705
# MEblue             0.064403546  0.89933476
# MEmidnightblue    -0.047564465  0.76141499
# MEsalmon4         -0.223052045  0.72512874
# MEdarkslateblue   -0.255552693  0.87626475
# MEantiquewhite4   -0.258817097  0.89501037
# MElightyellow     -0.169540060  0.92615633
# MEpalevioletred3  -0.004169589  0.33379816
# MEwhite           -0.057165878  0.60623759
# MEdarkmagenta     -0.078289235  0.80985669
# MEdarkgrey        -0.090051999  0.91162831
# MEplum2           -0.043344323  0.84926063
# MEmediumpurple3    0.094644761 -0.58786906
# MEnavajowhite2     0.046660019 -0.15682166
# MEmaroon           0.089589113  0.25887072
# MEskyblue          0.232092494 -0.13048241
# MEmediumorchid     0.353694886  0.28649173
# MEsaddlebrown      0.077444880  0.65091154
# MEgrey             0.108775146 -0.01828212

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
main = paste("Module-trait relationships (unmerged)"))


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
table(mergedColors)
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



#==============================================
#
# Module-trait analysis
#
#==============================================

# For the time being we will continue with all of the samples
# Define numbers of genes and samples
nGenes = ncol(wgcna_in);
nSamples = nrow(wgcna_in);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(wgcna_in, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# This shows strong associations to time, but not axenic
# when working with the c(3, 24, 48) samples:
# > moduleTraitCor
#                       axenic     time_hr
# MEmediumpurple3   0.07720808 -0.44170599
# MEcoral2          0.03420079 -0.98501585
# MEfloralwhite     0.16846036 -0.73523210
# MEblack          -0.10807242 -0.49709456
# MEsienna3        -0.12896334  0.64912676
# MEhoneydew1      -0.06870903  0.90156372
# MEmaroon          0.16684523  0.33415061
# MEdarkolivegreen  0.13665768  0.33889093
# MElavenderblush3 -0.08557305  0.47026786
# MEgrey            0.10877515 -0.01828212



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
main = paste("Module-trait relationships (merged)"))

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
geneTraitSignificance["PHATRDRAFT_43365",]
# > head(geneTraitSignificance %>% arrange(desc(abs(GS.axenic))),10)
#                   GS.axenic  p.GS.axenic   GS.time_hr p.GS.time_hr
# PHATRDRAFT_47845  0.8925782 6.333572e-07 -0.116677426    0.6447537
# PHATRDRAFT_46195 -0.8629767 4.028727e-06  0.201185036    0.4234132
# PHATRDRAFT_43365  0.8487865 8.456006e-06 -0.104868197    0.6787848
# PHATRDRAFT_43020  0.8481265 8.736711e-06 -0.172685526    0.4932044
# PHATRDRAFT_48069  0.8283735 2.175682e-05 -0.024390803    0.9234682
# PHATR_46647       0.8075654 5.067855e-05  0.283210420    0.2547829
# PHATRDRAFT_48554  0.7989125 6.997334e-05  0.032759618    0.8973236
# PHATRDRAFT_47597  0.7917863 9.024530e-05 -0.109146687    0.6663827
# PHATR_44100       0.7909817 9.281887e-05  0.004895515    0.9846187
# PHATRDRAFT_46056  0.7874332 1.049258e-04 -0.303422802    0.2209485

# > geneTraitSignificance["PHATRDRAFT_43365",]
#                  GS.axenic  p.GS.axenic GS.time_hr p.GS.time_hr
# PHATRDRAFT_43365 0.8487865 8.456006e-06 -0.1048682    0.6787848

# Get a connectivity score for each of the genes
adjacency_sum = rowSums(adjacency)
temp_merged = merge(geneTraitSignificance, adjacency_sum, by="row.names")
rownames(temp_merged) = temp_merged$Row.names
geneTraitSignificance_adj = temp_merged %>% mutate(adj_score=y) %>% select(-Row.names, -y) %>% arrange(desc(abs(GS.axenic)))
head(geneTraitSignificance_adj)
head(geneTraitSignificance_adj,10)
# > head(geneTraitSignificance_adj %>% arrange(desc(abs(GS.axenic))),10)
#                   GS.axenic  p.GS.axenic   GS.time_hr p.GS.time_hr adj_score
# PHATRDRAFT_47845  0.8925782 6.333572e-07 -0.116677426    0.6447537  3.722438
# PHATRDRAFT_46195 -0.8629767 4.028727e-06  0.201185036    0.4234132  4.618721
# PHATRDRAFT_43365  0.8487865 8.456006e-06 -0.104868197    0.6787848  7.060198
# PHATRDRAFT_43020  0.8481265 8.736711e-06 -0.172685526    0.4932044  4.937492
# PHATRDRAFT_48069  0.8283735 2.175682e-05 -0.024390803    0.9234682  6.855097
# PHATR_46647       0.8075654 5.067855e-05  0.283210420    0.2547829 10.687355
# PHATRDRAFT_48554  0.7989125 6.997334e-05  0.032759618    0.8973236  7.205187
# PHATRDRAFT_47597  0.7917863 9.024530e-05 -0.109146687    0.6663827  2.753997
# PHATR_44100       0.7909817 9.281887e-05  0.004895515    0.9846187 11.055521
# PHATRDRAFT_46056  0.7874332 1.049258e-04 -0.303422802    0.2209485  7.588502

geneTraitSignificance_adj["PHATRDRAFT_43365",]
# > geneTraitSignificance_adj["PHATRDRAFT_43365",]
#                  GS.axenic  p.GS.axenic GS.time_hr p.GS.time_hr adj_score
# PHATRDRAFT_43365 0.8487865 8.456006e-06 -0.1048682    0.6787848  7.060198

# TODO It's also useful to have a look at the module membership of PHATRDRAFT_43365
# to show that it is not a strong member of any of the modules.
# TODO we can also look to see if there are any genes to which PHATRDRAFT_43365
# has a particularly strong connection as these would be the genes that 
# it is likely involved in the regulation of.
geneModuleMembership["PHATRDRAFT_43365",]
# > geneModuleMembership["PHATRDRAFT_43365",]
#                  MMmediumpurple3  MMcoral2 MMfloralwhite     MMblack  MMsienna3
# PHATRDRAFT_43365        0.073501 0.1188302     0.2519848 -0.05693819 -0.1751843
#                  MMhoneydew1   MMmaroon MMdarkolivegreen MMlavenderblush3
# PHATRDRAFT_43365   -0.133679 0.02852935        0.1136577       -0.1022929
#                    MMgrey
# PHATRDRAFT_43365 0.512664
# THis is confusing because according to the modulemembership, PHATRDRAFT_43365 would be a member
# of the Grey group. However in the module assignment it is a member of the coral2 group.
# This is apparently due to how the module assignments are done (https://www.biostars.org/p/76611/)

# Now let's look at which genes PHATRDRAFT_43365 has a strong connection to
adj_plot_df = data.frame(adj=as.numeric(adjacency["PHATRDRAFT_43365",]))
row.names(adj_plot_df) = names(adjacency["PHATRDRAFT_43365",])
summary(adj_plot_df$adj)
head(adj_plot_df %>% dplyr::arrange(desc(adj)), 20)
adj_plot_df = adj_plot_df %>% mutate(gene = "PHATRDRAFT_43365", label=row.names(adj_plot_df))
ggplot(adj_plot_df, aes(x=gene, y=adj, lab=label)) + geom_point() + geom_label_repel(aes(label=ifelse(adj > .1, label, ""))) + ggtitle("Adjacency score to PHATRDRAFT_43365")
# > head(adj_plot_df %>% dplyr::arrange(desc(adj)), 20)
#                         adj
# PHATRDRAFT_43365 1.00000000
# PHATRDRAFT_48554 0.80303596
# PHATRDRAFT_48069 0.61900382
# PHATRDRAFT_55070 0.45620549
# PHATRDRAFT_46444 0.26431079
# PHATRDRAFT_50288 0.25122068
# PHATRDRAFT_33892 0.24997839
# PHATRDRAFT_43020 0.24476284
# PHATRDRAFT_45862 0.17656256
# PHATRDRAFT_43081 0.17503532
# PHATRDRAFT_49510 0.17010541
# PHATRDRAFT_40174 0.16466467
# PHATRDRAFT_33664 0.14630032
# PHATRDRAFT_48064 0.12164598
# PHATRDRAFT_48558 0.10943451
# PHATRDRAFT_43363 0.09803264
# PHATR_44100      0.08924753
# PHATRDRAFT_49590 0.08711789
# PHATRDRAFT_47845 0.07600548
# PHATRDRAFT_15935 0.07491525
# This is super interesting. We can compare this list to the list of 
# DE expressed genes and we see that there is a large overlap.
# I will output this list so that we can compare it to the DE genes
# I have made the figure here: /home/humebc/projects/ru/nextflow_ru/run_1_nf/DEgenes.adjacency.subset.png
# All in all this shows us that approximately half of the significant genes are related to PHATRDRAFT_43365 each other.
save(adj_plot_df, file="/home/humebc/projects/ru/nextflow_ru/run_1_nf/adj_plot_df.RData")

# next most signficacnt gene that wasn't related to PHATRDRAFT_43365 was PHATRDRAFT_44925.
# I want to do 2 things. I want to check its adjacency to PHATRDRAFT_43365.
# > adj_plot_df["PHATRDRAFT_44925",]
#                          adj             gene            label
# PHATRDRAFT_44925 0.003807441 PHATRDRAFT_43365 PHATRDRAFT_44925
# Then I want to pull up the adjacency scores of PHATRDRAFT_44925 and see if it relates to
# the other DE genes.
adj_plot_df_44925 = data.frame(adj=as.numeric(adjacency["PHATRDRAFT_44925",]))
row.names(adj_plot_df_44925) = names(adjacency["PHATRDRAFT_44925",])
adj_plot_df_44925 = adj_plot_df_44925 %>% mutate(gene = "PHATRDRAFT_44925", label=row.names(adj_plot_df_44925)) %>% arrange(desc(adj)) %>% dplyr::filter(adj>0.1)
summary(adj_plot_df_44925$adj)

ggplot(adj_plot_df_44925, aes(x=gene, y=adj, lab=label)) + geom_point() + geom_label_repel(aes(label=ifelse(adj > .1, label, ""))) + ggtitle("Adjacency score to PHATRDRAFT_44925")

# Read in the list of DE genes
load("/home/humebc/projects/ru/nextflow_ru/run_1_nf/de_genes.subset.RData")
de_genes_names = row.names(de_genes)
sum(row.names(adj_plot_df_44925) %in% de_genes_names)
# > sum(row.names(adj_plot_df_44925) %in% de_genes_names)
# [1] 1
# Of 44925's adjacent genes it is the only gene that is DE expressed.

# We can work up a quick method for this which will give us the number of genes
# that are DE expressed in those genes with an adjacency score > 0.1 for the given gene

get_de_network_size = function(gene_in_q){
    adj_plot_df_gene_in_q = data.frame(adj=as.numeric(adjacency[gene_in_q,]))
    row.names(adj_plot_df_gene_in_q) = names(adjacency[gene_in_q,])
    adj_plot_df_gene_in_q = adj_plot_df_gene_in_q %>% mutate(gene = gene_in_q, label=row.names(adj_plot_df_gene_in_q)) %>% arrange(desc(adj)) %>% dplyr::filter(adj>0.1)

    list(sum(row.names(adj_plot_df_gene_in_q) %in% de_genes_names), row.names(adj_plot_df_gene_in_q)[row.names(adj_plot_df_gene_in_q) %in% de_genes_names])
}
net_scores = lapply(as.list(de_genes_names), FUN=get_de_network_size)
names(net_scores) = de_genes_names
# This shows us that the DE genes are largely connected with each other
# and there appears to be some structure to the networks.
# From here we would like to do some graphs of the connectivity.
# I think for sake of keeping things clear
# I will do this in yet another R script (run_1_networks.r) in which I will 
# load the adjacency matrix and the de_genes df.
# We already have the de_genes saved here: "/home/humebc/projects/ru/nextflow_ru/run_1_nf/de_genes.subset.RData"
# and then we will save the adjacency matrix:
save(adjacency, file="/home/humebc/projects/ru/nextflow_ru/run_1_nf/adjacency.mat.subset.RData")

# Look up the high adjacency genes in the de_genes
row.names(adj_plot_df_gene_in_q)

# After that I want to see if it's possible to plot up this part of the network.


# Somewhat worryingly the adjacency score of the PHATRDRAFT_43365 is relatively low compared to 
# the dataset average
summary(geneTraitSignificance_adj)
save(geneTraitSignificance_adj, file="/home/humebc/projects/ru/nextflow_ru/run_1_nf/geneTraitSignificance_adj.subset.RData")

# > summary(geneTraitSignificance_adj)
#    GS.axenic          p.GS.axenic          GS.time_hr       
#  Min.   :-0.862977   Min.   :0.0000006   Min.   :-0.989673  
#  1st Qu.:-0.123098   1st Qu.:0.3847051   1st Qu.:-0.660033  
#  Median : 0.006599   Median :0.6151589   Median : 0.032198  
#  Mean   : 0.002956   Mean   :0.5842519   Mean   :-0.009227  
#  3rd Qu.: 0.130386   3rd Qu.:0.8141989   3rd Qu.: 0.615400  
#  Max.   : 0.892578   Max.   :0.9999322   Max.   : 0.995665  
#   p.GS.time_hr         adj_score       
#  Min.   :0.0000000   Min.   :   1.014  
#  1st Qu.:0.0000157   1st Qu.: 129.395  
#  Median :0.0046698   Median : 326.689  
#  Mean   :0.1366744   Mean   : 360.916  
#  3rd Qu.:0.1473624   3rd Qu.: 553.755  
#  Max.   :0.9997037   Max.   :1083.547 

# The final thing I want to look at is see which if the GS.axenic genes belong to
# a specific module or not. The module membership is held in moduleColors
names(moduleColors) = colnames(wgcna_in)
head(moduleColors)
geneTraitSignificance_adj_col = (merge(geneTraitSignificance_adj, moduleColors, by="row.names"))
rownames(geneTraitSignificance_adj_col) = geneTraitSignificance_adj_col$Row.names
geneTraitSignificance_adj_col = geneTraitSignificance_adj_col %>% arrange(desc(abs(GS.axenic))) %>% select(-Row.names)
head(geneTraitSignificance_adj_col, 20)
# > head(geneTraitSignificance_adj_col, 20)
#                   GS.axenic  p.GS.axenic   GS.time_hr p.GS.time_hr adj_score
# PHATRDRAFT_47845  0.8925782 6.333572e-07 -0.116677426    0.6447537  3.722438
# PHATRDRAFT_46195 -0.8629767 4.028727e-06  0.201185036    0.4234132  4.618721
# PHATRDRAFT_43365  0.8487865 8.456006e-06 -0.104868197    0.6787848  7.060198
# PHATRDRAFT_43020  0.8481265 8.736711e-06 -0.172685526    0.4932044  4.937492
# PHATRDRAFT_48069  0.8283735 2.175682e-05 -0.024390803    0.9234682  6.855097
# PHATR_46647       0.8075654 5.067855e-05  0.283210420    0.2547829 10.687355
# PHATRDRAFT_48554  0.7989125 6.997334e-05  0.032759618    0.8973236  7.205187
# PHATRDRAFT_47597  0.7917863 9.024530e-05 -0.109146687    0.6663827  2.753997
# PHATR_44100       0.7909817 9.281887e-05  0.004895515    0.9846187 11.055521
# PHATRDRAFT_46056  0.7874332 1.049258e-04 -0.303422802    0.2209485  7.588502
# PHATR_46881       0.7867420 1.074336e-04 -0.127726825    0.6135080  4.707569
# PHATRDRAFT_50288  0.7799556 1.348737e-04 -0.342088323    0.1646775 10.898140
# PHATRDRAFT_50116 -0.7737138 1.651339e-04 -0.286700932    0.2487204  9.625186
# PHATRDRAFT_43073  0.7732615 1.675333e-04  0.386672876    0.1129386 12.852649
# PHATRDRAFT_49309 -0.7663273 2.081752e-04  0.033894049    0.8937867  6.107471
# PHATRDRAFT_33892 -0.7548818 2.933904e-04  0.382252569    0.1174721  9.998348
# PHATRDRAFT_49510  0.7531664 3.083907e-04  0.264105384    0.2895943  9.010355
# PHATRDRAFT_42653  0.7470036 3.677204e-04 -0.124148009    0.6235623  2.637504
# PHATRDRAFT_9509   0.7389373 4.596193e-04 -0.252891997    0.3113065  9.161192
# PHATRDRAFT_48558  0.7364438 4.916511e-04 -0.360191947    0.1420306  9.297885
#                               y
# PHATRDRAFT_47845        sienna3
# PHATRDRAFT_46195        sienna3
# PHATRDRAFT_43365         coral2
# PHATRDRAFT_43020           grey
# PHATRDRAFT_48069         coral2
# PHATR_46647              coral2
# PHATRDRAFT_48554         coral2
# PHATRDRAFT_47597      honeydew1
# PHATR_44100              coral2
# PHATRDRAFT_46056        sienna3
# PHATR_46881         floralwhite
# PHATRDRAFT_50288    floralwhite
# PHATRDRAFT_50116          black
# PHATRDRAFT_43073          black
# PHATRDRAFT_49309 darkolivegreen
# PHATRDRAFT_33892    floralwhite
# PHATRDRAFT_49510         coral2
# PHATRDRAFT_42653        sienna3
# PHATRDRAFT_9509       honeydew1
# PHATRDRAFT_48558      honeydew1

# We can see that the high GS.axenic related genes are spread across several modules
# and when we look inside the coral2 module we see a negative correlation (see below) between group membership
# and gene significance for axenic. This is simply a repeat result from the module-trait heat map
# results.
# PHATRDRAFT_43365 is in module coral2
module_of_interest = "coral2"
column = match(module_of_interest, modNames)
moduleGenes = moduleColors==module_of_interest

plotting.df = data.frame(MM=abs(geneModuleMembership[moduleGenes, column]), GS=abs(geneTraitSignificance[moduleGenes, "GS.axenic"]))
row.names(plotting.df) = row.names(geneModuleMembership)[moduleGenes]
plotting.df$gene_names = row.names(geneModuleMembership)[moduleGenes]

ggplot(plotting.df, aes(x=MM, y=GS, label=gene_names)) + 
    geom_point(color=dplyr::case_when(plotting.df$GS > 0.75 ~ "red", TRUE ~ "black")) + 
    geom_label_repel(aes(label=ifelse(GS > .75, gene_names, ""))) + 
    xlab(paste("Module Membership in", module_of_interest, "module")) + ylab("Gene significance for axenic") + 
    ggtitle("Module membership vs. gene significance\n")

ggsave("/home/humebc/projects/ru/nextflow_ru/run_1_nf/coral2_MM_GS_scatter.subset.png")


# THe same negative correlation can be seen for sienna3
module_of_interest = "sienna3"
column = match(module_of_interest, modNames)
moduleGenes = moduleColors==module_of_interest

plotting.df = data.frame(MM=abs(geneModuleMembership[moduleGenes, column]), GS=abs(geneTraitSignificance[moduleGenes, "GS.axenic"]))
row.names(plotting.df) = row.names(geneModuleMembership)[moduleGenes]
plotting.df$gene_names = row.names(geneModuleMembership)[moduleGenes]

ggplot(plotting.df, aes(x=MM, y=GS, label=gene_names)) + 
    geom_point(color=dplyr::case_when(plotting.df$GS > 0.75 ~ "red", TRUE ~ "black")) + 
    geom_label_repel(aes(label=ifelse(GS > .75, gene_names, ""))) + 
    xlab(paste("Module Membership in", module_of_interest, "module")) + ylab("Gene significance for axenic") + 
    ggtitle("Module membership vs. gene significance\n")

ggsave("/home/humebc/projects/ru/nextflow_ru/run_1_nf/sienna3_MM_GS_scatter.subset.png")

#### THIS MEANS: that the WGCNA does not result in us identifying modules that relate to
# Axenic state.

save.image("/home/humebc/projects/ru/nextflow_ru/run_1_nf/wgcna.subset.RData")

# There are some really highly correllated genes with the axenic
# state, but that the modules that are being defined are not related
# strongly to axenic state but rather to time.
# We tried running the analyses on only the 24 hour samples
# but you end up with nonsense results. The recommended minimum
# number of samples is 20 and so the 6 is way below that. With the 
# 3 time points we are at least at 3*6 = 18 samples. If we were to
# include the 0.5 hour samples this would improve the sample number but I think
# the lack of effect at this time would jeaporadise our ability to detect signal.
# Moving forwards, I think we should try looking at the pre-merging modules
# to see if any of them have high significance to the axenic state.

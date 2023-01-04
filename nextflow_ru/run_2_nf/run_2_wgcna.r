# In this script, merely as a sanity check I will try to recreate
# the WGCNA module trait findings that Ru came up with.

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

# If the script has already been run through then you can use this load.
# load(file = "nextflow_ru/run_2_nf/run_1_data.RData")

tx2gene = read.table("/home/humebc/projects/ru/reference/tx2gene.txt", header=TRUE)
samples = read.csv("/home/humebc/projects/ru/nextflow_ru/run_2_nf/samples_run_2.csv", header=TRUE)
samples = samples %>% mutate(axenic = as.factor(axenic), time_hr = as.factor(time_hr), shaking = as.factor(shaking), cell_line = as.factor(cell_line)) %>% mutate(axenic = relevel(axenic, "TRUE"))
# Create a column that is 'mutant' so that we can easily do the comparisons that ru made for the shaking samples
samples = samples %>% mutate(mutant = as.factor(case_when(cell_line == "WT" ~ FALSE, TRUE ~ TRUE)))
rownames(samples) = samples$dir_name

# Make a vector that contains the full paths to the abundance.h5 files
kallisto.base.dir = "/home/humebc/projects/ru/nextflow_ru/run_2_nf/results/kallisto"

files <- file.path(kallisto.base.dir, samples$dir_name, "abundance.h5")

# Finally we can use tximport to read in the abundance tables
# and perform the normalizations
txi = tximport(files, type = "kallisto", tx2gene = tx2gene)

# Create the DESEQ2 object
dds = DESeqDataSetFromTximport(txi, colData = samples, design = ~ mutant)

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
traits = samples %>% select(axenic, time_hr, mutant) %>% mutate(axenic = as.numeric(axenic), time_hr = as.numeric(time_hr), mutant = as.numeric(mutant))

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
sft = pickSoftThreshold(wgcna_in, powerVector = powers, verbose = 5, blockSize=20000)

# Plot scale-free toplogy model fit R^2 value against power
text_size = 5
sft_plotting_df = sft$fitIndices
r_sqr_plot = ggplot(sft_plotting_df, aes(x=Power, y=-sign(slope)*SFT.R.sq, label=Power)) + 
geom_point(size=0.5) + geom_text(hjust=0, vjust=0, size=3, color=ifelse(sft_plotting_df$Power==8, 'red', 'black')) + ggtitle("Scale independence") + 
    theme(plot.title = element_text(size = text_size, face = "bold"), 
    axis.title=element_text(size=text_size,face="bold"), 
    axis.text=element_text(size=text_size)) + xlab("Soft Threshold (Power)") + geom_hline(yintercept=0.8, color="red")


mean_conn_plot = ggplot(sft_plotting_df, aes(x=Power, y=mean.k., label=Power)) + 
geom_point(size=0.5) + geom_text(hjust=0, vjust=0, size=3, color=ifelse(sft_plotting_df$Power==8, 'red', 'black')) + ggtitle("Mean connectivity") + 
    theme(plot.title = element_text(size = text_size, face = "bold"),
    axis.title=element_text(size=text_size,face="bold"), axis.text=element_text(size=text_size)) + 
    xlab("Soft Threshold (Power)")


threshold_plot = grid.arrange(r_sqr_plot, mean_conn_plot, ncol=2)
ggsave("nextflow_ru/run_2_nf/run_2_wgcna_threshold_selection.png", plot=threshold_plot, width=20, height=10, units="cm")

# We will work with the soft threshold power value of 16

###############################################
# MAKING THE NETWORK
###############################################

# Let's try enabling the  multithreading
enableWGCNAThreads(nThreads = 40) # As far as I can tell this doesn't do much


if(file.exists("/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_TOM.RData")){
    load("/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_TOM.RData")
}else{
    # Make adjacency matrix
    adjacency = adjacency(wgcna_in, power = 8)
    # Turn adjacency into topological overlap
    TOM = TOMsimilarity(adjacency)
    save(adjacency, TOM, file="/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_TOM.RData")
}

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

# There is a strong module association with mutant
# > moduleTraitCor
#                       axenic     time_hr        mutant
# MEfloralwhite     0.32810563  0.05839342  0.0552339685
# MEdarkslateblue   0.45697792  0.06358095  0.0501108387
# MEsienna3         0.15624843  0.12960521 -0.2121726348
# MEdarkolivegreen  0.18983135 -0.09689823  0.2209222595
# MEdarkorange2    -0.05431407 -0.55095508 -0.0005054042
# MEivory          -0.04592322 -0.46927334  0.1069277747
# MEdarkgreen      -0.24628798  0.11020883 -0.9940476738
# MEdarkorange     -0.53318372  0.10510987 -0.2047655855
# MEdarkred        -0.50272849 -0.40445257 -0.2037588832
# MEgreen           0.14298778  0.49418042  0.0868418987
# MElightgreen      0.12884870 -0.21118530  0.0051119573
# MEblack           0.40647406  0.28901818  0.0514247121
# MEroyalblue       0.38371646  0.39589160  0.3218840437
# MEgrey           -0.02168026  0.01144192 -0.0413383915

# > moduleTraitPvalue
#                        axenic      time_hr       mutant
# MEfloralwhite    1.874649e-02 6.839920e-01 7.002636e-01
# MEdarkslateblue  7.489863e-04 6.575841e-01 7.269312e-01
# MEsienna3        2.735535e-01 3.646968e-01 1.349853e-01
# MEdarkolivegreen 1.821294e-01 4.987671e-01 1.192554e-01
# MEdarkorange2    7.050267e-01 2.795852e-05 9.971916e-01
# MEivory          7.489715e-01 5.131903e-04 4.551655e-01
# MEdarkgreen      8.147493e-02 4.413674e-01 7.593857e-49
# MEdarkorange     5.617394e-05 4.629106e-01 1.494747e-01
# MEdarkred        1.702288e-04 3.243477e-03 1.515295e-01
# MEgreen          3.168378e-01 2.281389e-04 5.445501e-01
# MElightgreen     3.675267e-01 1.368535e-01 9.715998e-01
# MEblack          3.078477e-03 3.969018e-02 7.200599e-01
# MEroyalblue      5.442179e-03 4.031971e-03 2.125625e-02
# MEgrey           8.799699e-01 9.364850e-01 7.733318e-01


textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)

png(filename="/home/humebc/projects/ru/nextflow_ru/run_2_nf/mod_trait.png", width=10, height=10, units="cm", res=300)
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
dev.off()

# I'm going to quicly recode the traits to 0 and 1 to see if it makes a difference.
# Yeah, I did that and it makes no difference which is good.

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

mutant = as.data.frame(traits$mutant)
names(mutant) = "mutant"
geneTraitSignificance_mutant = as.data.frame(cor(wgcna_in, mutant, use = "p"))
GSPvalue_mutant = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_mutant), nSamples))
names(geneTraitSignificance_mutant) = paste("GS.", names(mutant), sep="")
names(GSPvalue_mutant) = paste("p.GS.", names(mutant), sep="")

# Put all four of the dataframes together into a single df
l = list(geneTraitSignificance_axenic, GSPvalue_axenic, geneTraitSignificance_time_hr, GSPvalue_time_hr, geneTraitSignificance_mutant, GSPvalue_mutant)
geneTraitSignificance = Reduce(merge, lapply(l, function(x) data.frame(x, rn = row.names(x))))
row.names(geneTraitSignificance) = geneTraitSignificance$rn
geneTraitSignificance = geneTraitSignificance %>% select(-rn)

# Search for the top 10 genes that relate to axenic state
head(geneTraitSignificance %>% arrange(desc(abs(GS.axenic))),10)

# > head(geneTraitSignificance %>% arrange(desc(abs(GS.axenic))),10)
#                   GS.axenic  p.GS.axenic  GS.time_hr p.GS.time_hr   GS.mutant   p.GS.mutant
# PHATRDRAFT_11154 -0.8225646 1.329514e-13 -0.28379426 4.357509e-02 -0.22278686  0.1160909878
# PHATRDRAFT_46444  0.7889601 6.140370e-12  0.53328702 5.595288e-05  0.06193450  0.6659227396
# PHATRDRAFT_36600 -0.7800941 1.505969e-11 -0.44783039 9.833776e-04 -0.34035529  0.0145328986
# PHATR_46834       0.7618552 8.426973e-11  0.38164853 5.720353e-03  0.47031357  0.0004967106
# PHATRDRAFT_45995 -0.7473818 2.974898e-10 -0.14446267 3.118196e-01 -0.41873414  0.0022274568
# PHATR_44183      -0.7424169 4.497716e-10 -0.54873408 3.057232e-05 -0.28988864  0.0390714084
# PHATRDRAFT_55004 -0.7412648 4.943911e-10  0.01991404 8.896851e-01 -0.17755899  0.2125728629
# PHATRDRAFT_43360 -0.7409723 5.063643e-10 -0.65486060 1.854777e-07 -0.23933044  0.0907472407
# PHATRDRAFT_2164  -0.7277233 1.449645e-09 -0.24723090 8.027782e-02 -0.44005997  0.0012320423
# PHATRDRAFT_54257 -0.7257136 1.691396e-09 -0.67453111 5.725440e-08 -0.09419078  0.5108949527

# There are genes that have a high correlation to the axenic state (and this is before we
# subset down to work with only the non-shaking)

# Let's quickly look up the DE genes from the 1st run in this
# geneTraitSig table
load("/home/humebc/projects/ru/nextflow_ru/run_1_nf/de_genes.subset.RData")
head(de_genes)
geneTraitSignificance[rownames(de_genes),]

# > geneTraitSignificance[rownames(de_genes),]
#                    GS.axenic  p.GS.axenic  GS.time_hr p.GS.time_hr    GS.mutant   p.GS.mutant
# PHATRDRAFT_43365  0.26646704 5.874320e-02  0.33469507 1.636682e-02 -0.075469244  0.5986461379
# PHATRDRAFT_48554  0.22050198 1.199777e-01  0.15714087 2.707886e-01 -0.174884047  0.2196552502
# PHATRDRAFT_50288  0.33593586 1.594878e-02  0.02787275 8.460551e-01 -0.025093635  0.8612453080
# NA                        NA           NA          NA           NA           NA            NA
# PHATRDRAFT_43020  0.63069864 7.026026e-07  0.38360985 5.456226e-03  0.008598289  0.9522489350
# PHATRDRAFT_48069  0.16899716 2.358154e-01 -0.08178488 5.683125e-01 -0.114163443  0.4250497899
# PHATR_44100       0.04876472 7.339932e-01 -0.16134722 2.580070e-01 -0.023327843  0.8709230789
# PHATRDRAFT_49510  0.35819214 9.855882e-03  0.34195377 1.404794e-02 -0.168446202  0.2373685047
# PHATRDRAFT_55070  0.06074226 6.719861e-01 -0.14814283 2.995208e-01 -0.152652608  0.2848829926
# PHATRDRAFT_47845  0.02018313 8.882039e-01 -0.27464736 5.112401e-02 -0.194669211  0.1710343221
# PHATRDRAFT_40174  0.11550288 4.196014e-01 -0.11438792 4.241339e-01 -0.296036072  0.0349238026
# PHATRDRAFT_48558  0.47907853 3.757201e-04  0.61590855 1.503053e-06  0.087788888  0.5401540404
# PHATRDRAFT_45862  0.37701988 6.388403e-03  0.53231414 5.806694e-05 -0.208739795  0.1415639183
# PHATRDRAFT_43363  0.28580034 4.204795e-02  0.27433609 5.139848e-02 -0.195825038  0.1684577510
# PHATRDRAFT_40731 -0.06303997 6.603195e-01  0.10123922 4.796360e-01 -0.156978978  0.2712887843
# PHATRDRAFT_49511  0.28759325 4.072048e-02  0.24560487 8.235092e-02 -0.378698394  0.0061386765
# PHATRDRAFT_50470  0.25657667 6.914724e-02  0.38685010 5.043192e-03 -0.357594897  0.0099883411
# PHATRDRAFT_48179  0.13609994 3.409483e-01 -0.12753153 3.724853e-01 -0.102407301  0.4745552280
# PHATRDRAFT_49001  0.21930345 1.220555e-01  0.51311259 1.180228e-04  0.073662129  0.6074546449
# PHATRDRAFT_50347  0.19971496 1.599938e-01  0.20752804 1.439422e-01 -0.160471412  0.2606342799
# PHATR_46647       0.44887199 9.537192e-04  0.17196913 2.275580e-01 -0.005411216  0.9699380116
# PHATRDRAFT_52325 -0.16349999 2.516250e-01  0.17492025 2.195583e-01 -0.478921346  0.0003776307
# PHATRDRAFT_48063  0.56389717 1.639627e-05  0.18334339 1.978079e-01 -0.034124715  0.8120929924
# Yeah this shows us that the network genes are really not doing the same function in this network
# I.e. being significantly linked to GS.axenic. We could go on to plot the network of the same genes
# but I'm not sure what that would achieve.

# Now let's have a look at that mutant trait
# head(geneTraitSignificance %>% arrange(desc(abs(GS.mutant))),10)
#                    GS.axenic p.GS.axenic  GS.time_hr p.GS.time_hr  GS.mutant  p.GS.mutant
# PHATRDRAFT_51129  -0.2423440  0.08663491  0.11601639    0.4175232 -0.9990094 6.683373e-68
# PHATRDRAFT_38705  -0.2443900  0.08392715  0.11969278    0.4028174 -0.9983870 1.020477e-62
# PHATRDRAFT_48287  -0.2344821  0.09768058  0.12771468    0.3717934 -0.9981560 2.704121e-61
# PHATRDRAFT_38702  -0.2200665  0.12072947  0.12559684    0.3798410 -0.9956438 3.686009e-52
# PHATRDRAFT_22332  -0.2331425  0.09966687  0.11303323    0.4296781 -0.9940404 7.822802e-49
# PHATRDRAFT_bd1122 -0.2590200  0.06644992  0.11411959    0.4252288 -0.9896579 5.456393e-43
# PHATRDRAFT_bd1296 -0.2590200  0.06644992  0.11411959    0.4252288 -0.9896579 5.456393e-43
# PHATRDRAFT_48286  -0.2313868  0.10231715  0.11306941    0.4295295 -0.9895854 6.470457e-43
# PHATRDRAFT_37231   0.2768382  0.04922571 -0.05942966    0.6786856  0.9862629 5.504699e-40
# PHATRDRAFT_46458  -0.2713304  0.05411067  0.07714220    0.5905420 -0.9822344 2.865060e-37
# WOW there are some incredibly stronly related genes.
# These will be very interesting to look at.

#==============================================
#
# Intra module investigation
#
#==============================================
# Let's look at the genes of darkgreen
# that had the strongest relation to the mutant factor
module_of_interest = "darkgreen"
column = match(module_of_interest, modNames)
moduleGenes = moduleColors==module_of_interest

plotting.df = data.frame(MM=abs(geneModuleMembership[moduleGenes, column]), GS=abs(geneTraitSignificance[moduleGenes, "GS.mutant"]))
row.names(plotting.df) = row.names(geneModuleMembership)[moduleGenes]
plotting.df$gene_names = row.names(geneModuleMembership)[moduleGenes]

ggplot(plotting.df, aes(x=MM, y=GS, label=gene_names)) + 
    geom_point(color=dplyr::case_when(plotting.df$GS > 0.9 ~ "red", TRUE ~ "black")) + 
    geom_label_repel(aes(label=ifelse(GS > 0.9, gene_names, ""))) + 
    xlab(paste("Module Membership in", module_of_interest, "module")) + ylab("Gene significance for axenic") + 
    ggtitle("Module membership vs. gene significance\n")


ggsave("/home/humebc/projects/ru/nextflow_ru/run_2_nf/dark_green_MM_GS_scatter.png")

# TODO we will want to get a list of the mutant DE genes to build a network of the
# mutant DE genes.

save.image("/home/humebc/projects/ru/nextflow_ru/run_2_nf/run_2_wgcna.RData")

# From here we should now subset to the non-shaking a rerun the wgcna analysis.
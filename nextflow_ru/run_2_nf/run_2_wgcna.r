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
library(pheatmap)

# If the script has already been run through then you can use this load.
# load(file = "nextflow_ru/run_2_nf/run_1_data.RData")

# variable to choose which transcriptome to work with
# values can be either "NCBI" or "ensembl"
transcriptome = "ensembl"

if (transcriptome == "ensembl") {
  tx2gene = read.table("/home/humebc/projects/ru/reference/alt_transcriptome/tx2gene.txt", header=TRUE)
}
if (transcriptome == "NCBI") {
  tx2gene = read.table("/home/humebc/projects/ru/reference/tx2gene.txt", header=TRUE)
}
samples = read.csv("/home/humebc/projects/ru/nextflow_ru/run_2_nf/samples_run_2.csv", header=TRUE)
samples = samples %>% mutate(axenic = as.factor(axenic), time_hr = as.factor(time_hr), shaking = as.factor(shaking), cell_line = as.factor(cell_line)) %>% mutate(axenic = relevel(axenic, "TRUE"))
# Create a column that is 'mutant' so that we can easily do the comparisons that ru made for the shaking samples
samples = samples %>% mutate(mutant = as.factor(case_when(cell_line == "WT" ~ FALSE, TRUE ~ TRUE)))
rownames(samples) = samples$dir_name

# Make a vector that contains the full paths to the abundance.h5 files
# Make a vector that contains the full paths to the abundance.h5 files
if (transcriptome == "ensembl") {
  kallisto.base.dir = "/home/humebc/projects/ru/nextflow_ru/run_2_nf/results_alt_ref/kallisto"
}
if (transcriptome == "NCBI") {
  kallisto.base.dir = "/home/humebc/projects/ru/nextflow_ru/run_2_nf/results/kallisto"
}

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

if (transcriptome == "ensembl") {
  ggsave("nextflow_ru/run_2_nf/run_2_wgcna_threshold_selection.ensembl.png", plot=threshold_plot, width=20, height=10, units="cm")
}
if (transcriptome == "NCBI") {
  ggsave("nextflow_ru/run_2_nf/run_2_wgcna_threshold_selection.png", plot=threshold_plot, width=20, height=10, units="cm")
}


# We will work with the soft threshold power value of 16

###############################################
# MAKING THE NETWORK
###############################################

# Let's try enabling the  multithreading
enableWGCNAThreads(nThreads = 40) # As far as I can tell this doesn't do much

if (transcriptome == "ensembl") {
    if(file.exists("/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_TOM.ensembl.RData")){
        load("/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_TOM.ensembl.RData")
    }else{
        # Make adjacency matrix
        adjacency = adjacency(wgcna_in, power = 8)
        # Turn adjacency into topological overlap
        TOM = TOMsimilarity(adjacency)
        save(adjacency, TOM, file="/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_TOM.ensembl.RData")
    }
}
if (transcriptome == "NCBI") {
  if(file.exists("/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_TOM.RData")){
    load("/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_TOM.RData")
    }else{
        # Make adjacency matrix
        adjacency = adjacency(wgcna_in, power = 8)
        # Turn adjacency into topological overlap
        TOM = TOMsimilarity(adjacency)
        save(adjacency, TOM, file="/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_TOM.RData")
    }
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
if (transcriptome == "ensembl") {
  png(filename="/home/humebc/projects/ru/nextflow_ru/run_2_nf/mod_trait.ensembl.png", width=10, height=10, units="cm", res=300)
}
if (transcriptome == "NCBI") {
  png(filename="/home/humebc/projects/ru/nextflow_ru/run_2_nf/mod_trait.png", width=10, height=10, units="cm", res=300)
}

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

# > head(geneTraitSignificance %>% arrange(desc(abs(GS.axenic))),10)
#                 GS.axenic  p.GS.axenic  GS.time_hr p.GS.time_hr   GS.mutant  p.GS.mutant
# Phatr3_J11154  -0.8137122 3.929623e-13 -0.30278300 3.079793e-02 -0.24561111 0.0823428921
# Phatr3_J46444   0.7883787 6.520909e-12  0.52964027 6.426031e-05  0.06522593 0.6492939772
# Phatr3_J14898  -0.7807943 1.405067e-11 -0.18627365 1.906123e-01 -0.27779505 0.0484148471
# Phatr3_J36600  -0.7758469 2.281214e-11 -0.44592451 1.039794e-03 -0.34077698 0.0144035948
# Phatr3_J46834   0.7666027 5.464081e-11  0.38876048 4.812733e-03  0.46214257 0.0006401033
# Phatr3_J14442  -0.7470702 3.053927e-10 -0.47074576 4.900053e-04 -0.43440627 0.0014468694
# Phatr3_J55004  -0.7414034 4.888096e-10  0.01556316 9.136829e-01 -0.17615601 0.2162673168
# Phatr3_J44183  -0.7374532 6.736710e-10 -0.56245525 1.742004e-05 -0.28487797 0.0427445703
# Phatr3_J45995  -0.7311611 1.109957e-09 -0.11757618 4.112469e-01 -0.41424314 0.0025112539
# Phatr3_EG00253  0.7306803 1.152472e-09  0.56510758 1.558007e-05  0.30558872 0.0292058591


# There are genes that have a high correlation to the axenic state (and this is before we
# subset down to work with only the non-shaking)

# Let's quickly look up the DE genes from the 1st run in this
# geneTraitSig table
if (transcriptome == "ensembl") {
  load("/home/humebc/projects/ru/nextflow_ru/run_1_nf/de_genes.subset.ensembl.RData")
}
if (transcriptome == "NCBI") {
  load("/home/humebc/projects/ru/nextflow_ru/run_1_nf/de_genes.subset.RData")
}

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

# > geneTraitSignificance[rownames(de_genes),]
#                  GS.axenic  p.GS.axenic  GS.time_hr p.GS.time_hr     GS.mutant  p.GS.mutant
# Phatr3_J43365   0.26373660 0.0614807247  0.32873364 1.850773e-02 -0.0797566413 5.779761e-01
# Phatr3_J48554   0.20635086 0.1462808675  0.14786472 3.004391e-01 -0.1838007889 1.966722e-01
# Phatr3_J46195  -0.01668859 0.9074672557 -0.22944451 1.053120e-01 -0.0259336083 8.566487e-01
# Phatr3_EG00672  0.41812891 0.0022639627  0.12089952 3.980568e-01  0.5267662508 7.158803e-05
# Phatr3_J48069   0.16934014 0.2348521289 -0.08704310 5.436146e-01 -0.1132814983 4.286590e-01
# Phatr3_J47845   0.02095407 0.8839623472 -0.28388440 4.350552e-02 -0.1861103882 1.910082e-01
# Phatr3_J50288   0.33604991 0.0159108118  0.02359677 8.694480e-01 -0.0238685401 8.679576e-01
# NA                      NA           NA          NA           NA            NA           NA
# Phatr3_J46647   0.44782866 0.0009834276  0.16294742 2.532529e-01 -0.0009230026 9.948711e-01
# Phatr3_J55070   0.06047727 0.6733365645 -0.14941434 2.953454e-01 -0.1578748956 2.685284e-01
# Phatr3_J43073  -0.16212184 0.2556982046 -0.24149967 8.777211e-02 -0.1526939723 2.847509e-01
# Phatr3_J44100   0.04729894 0.7417078861 -0.16979494 2.335790e-01 -0.0222953533 8.765906e-01
# Phatr3_J48558   0.48287400 0.0003321621  0.61870072 1.305943e-06  0.0806156975 5.738740e-01
 

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

# > head(geneTraitSignificance %>% arrange(desc(abs(GS.mutant))),10)
#                 GS.axenic p.GS.axenic GS.time_hr p.GS.time_hr  GS.mutant  p.GS.mutant
# Phatr3_J48287  -0.2332579  0.09949459 0.12524051    0.3812051 -0.9983112 3.140393e-62
# Phatr3_J38705  -0.2441490  0.08424250 0.11866120    0.4069131 -0.9983083 3.275919e-62
# Phatr3_J29456  -0.2497614  0.07713357 0.10905029    0.4462128 -0.9981838 1.864672e-61
# Phatr3_J38702  -0.2192781  0.12209978 0.12520374    0.3813460 -0.9955431 6.447726e-52
# Phatr3_J48286  -0.2319181  0.10150945 0.11070932    0.4392833 -0.9876997 3.735390e-41
# Phatr3_J22332  -0.2294303  0.10533417 0.08292720    0.5629031 -0.9871727 1.037743e-40
# Phatr3_J46458  -0.2740757  0.05162895 0.07575044    0.5972805 -0.9815074 7.590357e-37
# Phatr3_EG01401 -0.3070454  0.02840672 0.02096898    0.8838803 -0.9810159 1.435336e-36
# Phatr3_J34902  -0.3049622  0.02955528 0.05617104    0.6954232 -0.9753249 8.289199e-34
# Phatr3_EG01364 -0.2221807  0.11711255 0.06554889    0.6476713 -0.9687093 2.588301e-31

# WOW there are some incredibly strongly related genes.
# These will be very interesting to look at.

#==============================================
#
# Intra module investigation
#
#==============================================


if (transcriptome == "ensembl") {
    # Let's look at the genes of white
    # that had the strongest relation to the mutant factor
    module_of_interest = "white"
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
  ggsave("/home/humebc/projects/ru/nextflow_ru/run_2_nf/dark_green_MM_GS_scatter.ensembl.png")
}
if (transcriptome == "NCBI") {
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
}


# From here we should now subset to the non-shaking a rerun the wgcna analysis.

# =================================

# Non-Shaking

# =================================

samples = read.csv("/home/humebc/projects/ru/nextflow_ru/run_2_nf/samples_run_2.csv", header=TRUE)
samples = samples %>% mutate(axenic = as.factor(axenic), time_hr = as.factor(time_hr), shaking = as.factor(shaking), cell_line = as.factor(cell_line)) %>% mutate(axenic = relevel(axenic, "TRUE"))
# Create a column that is 'mutant' so that we can easily do the comparisons that ru made for the shaking samples
samples = samples %>% mutate(mutant = as.factor(case_when(cell_line == "WT" ~ FALSE, TRUE ~ TRUE)))
rownames(samples) = samples$dir_name
samples = samples %>% dplyr::filter(shaking!="TRUE")

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
geom_point(size=0.5) + geom_text(hjust=0, vjust=0, size=3, color=ifelse(sft_plotting_df$Power==7, 'red', 'black')) + ggtitle("Scale independence") + 
    theme(plot.title = element_text(size = text_size, face = "bold"), 
    axis.title=element_text(size=text_size,face="bold"), 
    axis.text=element_text(size=text_size)) + xlab("Soft Threshold (Power)") + geom_hline(yintercept=0.8, color="red")


mean_conn_plot = ggplot(sft_plotting_df, aes(x=Power, y=mean.k., label=Power)) + 
geom_point(size=0.5) + geom_text(hjust=0, vjust=0, size=3, color=ifelse(sft_plotting_df$Power==7, 'red', 'black')) + ggtitle("Mean connectivity") + 
    theme(plot.title = element_text(size = text_size, face = "bold"),
    axis.title=element_text(size=text_size,face="bold"), axis.text=element_text(size=text_size)) + 
    xlab("Soft Threshold (Power)")


threshold_plot = grid.arrange(r_sqr_plot, mean_conn_plot, ncol=2)

if (transcriptome == "ensembl") {
  ggsave("nextflow_ru/run_2_nf/run_2_wgcna_threshold_selection.non_shaking.ensembl.png", plot=threshold_plot, width=20, height=10, units="cm")
}
if (transcriptome == "NCBI") {
  ggsave("nextflow_ru/run_2_nf/run_2_wgcna_threshold_selection.non_shaking.png", plot=threshold_plot, width=20, height=10, units="cm")
}


# We will work with the soft threshold power value of 16

###############################################
# MAKING THE NETWORK
###############################################

# Let's try enabling the  multithreading
enableWGCNAThreads(nThreads = 40) # As far as I can tell this doesn't do much

if (transcriptome == "ensembl") {
  if(file.exists("/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_TOM.non_shaking.ensembl.RData")){
        load("/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_TOM.non_shaking.ensembl.RData")
    }else{
        # Make adjacency matrix
        adjacency = adjacency(wgcna_in, power = 7)
        # Turn adjacency into topological overlap
        TOM = TOMsimilarity(adjacency)
        save(adjacency, TOM, file="/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_TOM.non_shaking.ensembl.RData")
    }
}
if (transcriptome == "NCBI") {
  if(file.exists("/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_TOM.non_shaking.RData")){
    load("/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_TOM.non_shaking.RData")
    }else{
        # Make adjacency matrix
        adjacency = adjacency(wgcna_in, power = 7)
        # Turn adjacency into topological overlap
        TOM = TOMsimilarity(adjacency)
        save(adjacency, TOM, file="/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_TOM.non_shaking.RData")
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


# Run the module-trait analysis on the non-merged modules
nGenes = ncol(wgcna_in);
nSamples = nrow(wgcna_in);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(wgcna_in, colors = dynamicColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# There are several modules that correspond strongly with axenic and time_hr,
# but the most exciting is the mutant lightcycan module -.99 association.
# > moduleTraitCor
#                     axenic     time_hr      mutant
# MEblue         -0.66138779 -0.62792553 -0.35080554
# MEbrown        -0.60673176 -0.78950119 -0.08686924
# MEblack        -0.81820681 -0.82410245 -0.20165701
# MEpurple       -0.83645583 -0.73839180 -0.21074968
# MEcyan         -0.70597470 -0.91738453  0.09642685
# MElightgreen   -0.85065861 -0.89626055 -0.14490431
# MEmagenta       0.84758972  0.59167789  0.55155124
# MEtan           0.68947556  0.36366179  0.39741255
# MEgreen         0.10132240 -0.09774425 -0.06175185
# MEgreenyellow   0.77674430  0.58434678  0.16769710
# MEpink          0.40507827  0.20704196 -0.03028781
# MEred           0.61706788  0.45794612  0.03983328
# MElightyellow  -0.17298699  0.13577048 -0.11478966
# MEturquoise     0.42110276  0.45826915  0.24191219
# MEyellow        0.19022260  0.29581417  0.14579474
# MEmidnightblue -0.18213012  0.32553835 -0.43125425
# MEsalmon        0.60276485  0.97482665 -0.09936760
# MEdarkred       0.10060871  0.56019501 -0.25352109
# MEroyalblue     0.43546486  0.83055025 -0.22856661
# MEgrey60       -0.13542160  0.16293482 -0.26119591
# MElightcyan    -0.25095785  0.21521749 -0.99141946
# MEgrey          0.04095908 -0.01070136 -0.06473964

textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)

if (transcriptome == "ensembl") {
  png(filename="/home/humebc/projects/ru/nextflow_ru/run_2_nf/mod_trait.non_shaking.unmerged.ensembl.png", width=10, height=30, units="cm", res=300)
}
if (transcriptome == "NCBI") {
  png(filename="/home/humebc/projects/ru/nextflow_ru/run_2_nf/mod_trait.non_shaking.unmerged.png", width=10, height=30, units="cm", res=300)
}

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
dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(wgcna_in, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result

# Look at the module-traits before merger


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
merge = merge$newMEs;

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors = mergedColors
names(moduleColors) = colnames(wgcna_in)
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(200));
moduleLabels = match(moduleColors, colorOrder)-1;
names(moduleLabels) = colnames(wgcna_in)
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

# For ensembl
# > moduleTraitCor
#                    axenic     time_hr      mutant
# MEdarkgrey     0.82737793  0.51236894  0.51414400
# MEbrown        0.33647903  0.15836995 -0.04909223
# MEorange       0.77801385  0.61956557  0.13913676
# MEcyan        -0.20742177  0.13110963 -0.27575686
# MEblue         0.29390188  0.39012161  0.18417202
# MEdarkred     -0.17335848  0.10610972 -0.05230563
# MEdarkgreen    0.14466837  0.62854596 -0.29369813
# MEpurple       0.43487055  0.88102036 -0.15045026
# MEgrey60      -0.71479103 -0.97227422  0.08440146
# MEturquoise   -0.68429979 -0.70100917 -0.28041151
# MElightcyan   -0.24196793  0.18500550 -0.99499488
# MElightyellow -0.30697785  0.14024879 -0.57208582
# MEgrey        -0.07734192  0.07841088 -0.10212585

# > moduleTraitPvalue
#                     axenic      time_hr       mutant
# MEdarkgrey    2.900734e-09 2.300398e-03 2.207515e-03
# MEbrown       5.553801e-02 3.787180e-01 7.861581e-01
# MEorange      9.899174e-08 1.206583e-04 4.399791e-01
# MEcyan        2.467605e-01 4.670575e-01 1.203482e-01
# MEblue        9.688778e-02 2.480430e-02 3.048899e-01
# MEdarkred     3.346444e-01 5.567245e-01 7.725163e-01
# MEdarkgreen   4.218296e-01 8.964550e-05 9.713030e-02
# MEpurple      1.143601e-02 1.339605e-11 4.033134e-01
# MEgrey60      2.965124e-06 3.989340e-21 6.405076e-01
# MEturquoise   1.126876e-05 5.532781e-06 1.139628e-01
# MElightcyan   1.748879e-01 3.026681e-01 1.395810e-32
# MElightyellow 8.226497e-02 4.362967e-01 5.044210e-04
# MEgrey        6.687924e-01 6.644794e-01 5.717137e-01

textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)

if (transcriptome == "ensembl") {
  png(filename="/home/humebc/projects/ru/nextflow_ru/run_2_nf/mod_trait.non_shaking.merged.ensembl.png", width=10, height=20, units="cm", res=300)
}
if (transcriptome == "NCBI") {
  png(filename="/home/humebc/projects/ru/nextflow_ru/run_2_nf/mod_trait.non_shaking.merged.png", width=10, height=10, units="cm", res=300)
}

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
#                   GS.axenic  p.GS.axenic GS.time_hr p.GS.time_hr   GS.mutant p.GS.mutant
# PHATRDRAFT_47587  0.9127855 1.363778e-13  0.6897852 8.967934e-06  0.46717692 0.006124111
# PHATR_43880       0.9001445 1.015947e-12  0.6167233 1.323042e-04  0.36985778 0.034130812
# PHATRDRAFT_48834  0.8997540 1.076303e-12  0.5923359 2.816160e-04  0.43685028 0.011024884
# PHATRDRAFT_39785  0.8978971 1.411685e-12  0.6565010 3.340825e-05  0.46797829 0.006025394
# PHATRDRAFT_47701  0.8978286 1.425726e-12  0.6861899 1.042172e-05  0.40885204 0.018156245
# PHATRDRAFT_35749 -0.8929081 2.853648e-12 -0.6065235 1.828310e-04 -0.38172218 0.028376907
# PHATRDRAFT_37590 -0.8914124 3.500406e-12 -0.6766276 1.538578e-05 -0.22743057 0.203063679
# PHATRDRAFT_43627  0.8905741 3.919960e-12  0.7621803 2.554289e-07  0.22207832 0.214187212
# PHATRDRAFT_48608  0.8904519 3.984845e-12  0.7929076 3.773588e-08  0.24584361 0.167869476
# PHATRDRAFT_43504  0.8895021 4.524490e-12  0.7211947 2.191803e-06  0.06830802 0.705645199

# ensembl
# > head(geneTraitSignificance %>% arrange(desc(abs(GS.axenic))),10)
#                 GS.axenic  p.GS.axenic GS.time_hr p.GS.time_hr  GS.mutant p.GS.mutant
# Phatr3_J43880   0.9018911 7.825133e-13  0.6200676 1.186992e-04  0.3644888  0.03702910
# Phatr3_J48834   0.8990015 1.202118e-12  0.5900756 3.011233e-04  0.4371800  0.01095764
# Phatr3_J47701   0.8986175 1.271461e-12  0.6831462 1.181590e-05  0.4086865  0.01820773
# Phatr3_J39785   0.8971660 1.568563e-12  0.6542789 3.626765e-05  0.4691582  0.00588253
# Phatr3_J16606  -0.8948104 2.191072e-12 -0.6873924 9.913297e-06 -0.2657641  0.13495006
# Phatr3_J37590  -0.8946650 2.236152e-12 -0.6825401 1.211293e-05 -0.2282120  0.20147383
# Phatr3_EG01075 -0.8940693 2.430030e-12 -0.7793294 9.117992e-08 -0.2906915  0.10076290
# Phatr3_J48608   0.8924714 3.029974e-12  0.7870136 5.577636e-08  0.2528054  0.15577159
# Phatr3_J43627   0.8908994 3.751914e-12  0.7584304 3.163596e-07  0.2293316  0.19921088
# Phatr3_J35749  -0.8897723 4.364495e-12 -0.6053063 1.898869e-04 -0.3807889  0.02879863

head(geneTraitSignificance %>% arrange(desc(abs(GS.time_hr))),10)

# > head(geneTraitSignificance %>% arrange(desc(abs(GS.time_hr))),10)
#                   GS.axenic  p.GS.axenic GS.time_hr p.GS.time_hr   GS.mutant p.GS.mutant
# PHATRDRAFT_30246 -0.6281030 9.098852e-05 -0.9896561 1.036710e-27  0.12758964  0.47920268
# PHATRDRAFT_50128 -0.4880038 3.963913e-03 -0.9632017 3.015117e-19  0.15480789  0.38966973
# PHATRDRAFT_48594  0.5599683 7.024188e-04  0.9617520 5.433077e-19 -0.08012210  0.65759716
# PHATRDRAFT_45944  0.5452027 1.034127e-03  0.9585462 1.850173e-18 -0.25542637  0.15138410
# PHATRDRAFT_48827  0.5545863 8.104426e-04  0.9575077 2.695549e-18 -0.31249557  0.07663127
# PHATR_10640      -0.5966858 2.472122e-04 -0.9574761 2.726165e-18 -0.05416375  0.76465813
# PHATRDRAFT_46547  0.5514626 8.796385e-04  0.9571220 3.092412e-18 -0.09906190  0.58336435
# PHATRDRAFT_35781  0.5220846 1.830829e-03  0.9561521 4.344645e-18 -0.22810762  0.20168563
# PHATRDRAFT_39303  0.5878418 3.215734e-04  0.9542218 8.358535e-18 -0.05244487  0.77192662
# PHATRDRAFT_54465  0.4979336 3.190663e-03  0.9513060 2.132946e-17 -0.23338634  0.19116342

# ensembl
# > head(geneTraitSignificance %>% arrange(desc(abs(GS.time_hr))),10)
#                 GS.axenic  p.GS.axenic GS.time_hr p.GS.time_hr   GS.mutant  p.GS.mutant
# Phatr3_J30246  -0.6289919 8.831134e-05 -0.9896970 9.751856e-28  0.12640614   0.48332264
# Phatr3_EG00191  0.6327498 7.775671e-05  0.9753234 6.695002e-22 -0.11926233   0.50857455
# Phatr3_J50128  -0.5012845 2.961043e-03 -0.9654571 1.149017e-19  0.16344569   0.36342747
# Phatr3_J48594   0.5598958 7.037840e-04  0.9624501 4.103407e-19 -0.08774921   0.62726301
# Phatr3_EG01286 -0.6607144 2.853808e-05 -0.9608077 7.877350e-19 -0.02723252   0.88042175
# Phatr3_J10640  -0.6056423 1.879152e-04 -0.9582007 2.099150e-18 -0.06013137   0.73957948
# Phatr3_J46547   0.5519388 8.687650e-04  0.9578777 2.359831e-18 -0.09908786   0.58326518
# Phatr3_J45944   0.5442809 1.058782e-03  0.9574646 2.737395e-18 -0.25643256   0.14972374
# Phatr3_J48827   0.5520097 8.671566e-04  0.9560144 4.556623e-18 -0.31473509   0.07443189
# Phatr3_J35781   0.5210511 1.876438e-03  0.9558957 4.747057e-18 -0.23018505   0.19749778

head(geneTraitSignificance %>% arrange(desc(abs(GS.mutant))),10)

# > head(geneTraitSignificance %>% arrange(desc(abs(GS.mutant))),10)
#                    GS.axenic p.GS.axenic  GS.time_hr p.GS.time_hr  GS.mutant  p.GS.mutant
# PHATRDRAFT_38705  -0.2612929  0.14188922  0.14159578    0.4318589 -0.9993960 8.381588e-47
# PHATRDRAFT_51129  -0.2703212  0.12813763  0.13152476    0.4656359 -0.9991485 1.715393e-44
# PHATRDRAFT_48287  -0.2596035  0.14457773  0.14879111    0.4085788 -0.9983304 5.809521e-40
# PHATRDRAFT_38702  -0.2512932  0.15834432  0.14153854    0.4320470 -0.9960969 2.979177e-34
# PHATRDRAFT_22332  -0.2419560  0.17490979  0.13737381    0.4458516 -0.9948453 2.201117e-32
# PHATRDRAFT_bd1122 -0.2989506  0.09102431  0.10856106    0.5475928 -0.9934097 9.821976e-31
# PHATRDRAFT_bd1296 -0.2989506  0.09102431  0.10856106    0.5475928 -0.9934097 9.821976e-31
# PHATRDRAFT_48286  -0.2854009  0.10740336  0.11889140    0.5099034 -0.9918792 2.474410e-29
# PHATRDRAFT_37231   0.3214944  0.06809014 -0.03968021    0.8264612  0.9857092 1.512491e-25
# PHATRDRAFT_46458  -0.3381786  0.05423859  0.09390500    0.6032071 -0.9786558 7.229616e-23

# ensembl
# > head(geneTraitSignificance %>% arrange(desc(abs(GS.mutant))),10)
#                 GS.axenic p.GS.axenic GS.time_hr p.GS.time_hr  GS.mutant  p.GS.mutant
# Phatr3_J38705  -0.2612746  0.14191814 0.14066272    0.4349303 -0.9994107 5.719960e-47
# Phatr3_J48287  -0.2579595  0.14722941 0.14489570    0.4210928 -0.9984746 1.434762e-40
# Phatr3_J29456  -0.2841466  0.10902504 0.12216881    0.4982221 -0.9978986 2.048804e-38
# Phatr3_J38702  -0.2512208  0.15846832 0.14126597    0.4329432 -0.9961998 1.970178e-34
# Phatr3_J48286  -0.2887665  0.10314203 0.11637473    0.5189648 -0.9903939 3.309194e-28
# Phatr3_J22332  -0.2473811  0.16514167 0.10821047    0.5488945 -0.9872289 2.675402e-26
# Phatr3_EG01364 -0.3425891  0.05097887 0.11412219    0.5271417 -0.9784237 8.535615e-23
# Phatr3_J46458  -0.3409109  0.05220036 0.09204654    0.6104281 -0.9779988 1.151416e-22
# Phatr3_EG01401 -0.3830436  0.02778849 0.02173042    0.9044575 -0.9754413 6.221084e-22
# Phatr3_J35054  -0.1341589  0.45666915 0.23502833    0.1879702 -0.9710818 7.599529e-21 

# There are some strong traits for each of the traits but mutatnt the strongest by far

# Let's quickly look up the DE genes from the 1st run in this
# geneTraitSig table
if (transcriptome == "ensembl") {
  load("/home/humebc/projects/ru/nextflow_ru/run_1_nf/de_genes.subset.ensembl.RData")
}
if (transcriptome == "NCBI") {
  load("/home/humebc/projects/ru/nextflow_ru/run_1_nf/de_genes.subset.RData")
}

head(de_genes)
geneTraitSignificance[rownames(de_genes),]
# For Ru, create a dataframe that holds the geneTraitSignificance and module label and module color in it
geneTraitSignificance_df = geneTraitSignificance %>% dplyr::mutate(module_color=moduleColors[rownames(geneTraitSignificance)], module_label=moduleLabels[rownames(geneTraitSignificance)])

# > geneTraitSignificance[rownames(de_genes),]
#                    GS.axenic  p.GS.axenic  GS.time_hr p.GS.time_hr   GS.mutant  p.GS.mutant
# PHATRDRAFT_43365  0.29288983 9.809693e-02  0.30015071 8.967125e-02  0.02617257 8.850446e-01
# PHATRDRAFT_48554  0.06423650 7.224777e-01  0.07149108 6.925800e-01 -0.20078832 2.625332e-01
# PHATRDRAFT_50288 -0.05202808 7.736919e-01  0.20940349 2.421737e-01 -0.18303600 3.079348e-01
# NA                        NA           NA          NA           NA          NA           NA
# PHATRDRAFT_43020  0.51704376 2.062796e-03  0.76946730 1.666553e-07 -0.16925271 3.463920e-01
# PHATRDRAFT_48069 -0.21105906 2.383857e-01 -0.37245916 3.279435e-02 -0.28942978 1.023175e-01
# PHATR_44100      -0.14466152 4.218518e-01 -0.24066543 1.772927e-01  0.06633993 7.137652e-01
# PHATRDRAFT_49510  0.23236253 1.931736e-01  0.46235262 6.748113e-03 -0.19956081 2.655228e-01
# PHATRDRAFT_55070 -0.28742107 1.048300e-01 -0.37463831 3.170788e-02 -0.12359020 4.931984e-01
# PHATRDRAFT_47845 -0.34566747 4.879756e-02 -0.30187426 8.775509e-02 -0.34974908 4.602089e-02
# PHATRDRAFT_40174 -0.14301101 4.272232e-01 -0.06728727 7.098527e-01 -0.34784871 4.729752e-02
# PHATRDRAFT_48558  0.64074064 5.898598e-05  0.64164261 5.714673e-05  0.24896028 1.623732e-01
# PHATRDRAFT_45862  0.53885388 1.214687e-03  0.57228647 5.016088e-04 -0.06067233 7.373187e-01
# PHATRDRAFT_43363  0.22970774 1.984546e-01  0.25462931 1.527088e-01 -0.31985433 6.958860e-02
# PHATRDRAFT_40731  0.02525398 8.890540e-01 -0.04778123 7.917418e-01 -0.01377047 9.393725e-01
# PHATRDRAFT_49511  0.05287112 7.701224e-01  0.38196711 2.826708e-02 -0.66490722 2.433719e-05
# PHATRDRAFT_50470  0.25914223 1.453181e-01  0.49212066 3.625668e-03 -0.34164177 5.166554e-02
# PHATRDRAFT_48179 -0.37309173 3.247589e-02  0.10553681 5.588688e-01 -0.48894815 3.884012e-03
# PHATRDRAFT_49001  0.29533609 9.519369e-02  0.48024451 4.675892e-03  0.13676319 4.478955e-01
# PHATRDRAFT_50347 -0.09943047 5.819573e-01  0.28983184 1.018202e-01 -0.28175038 1.121739e-01
# PHATR_46647       0.31520274 7.397888e-02  0.52180506 1.843069e-03 -0.22099707 2.164841e-01
# PHATRDRAFT_52325 -0.21939029 2.199284e-01  0.05851453 7.463494e-01 -0.44964755 8.656955e-03
# PHATRDRAFT_48063  0.56671151 5.851326e-04  0.45176363 8.310424e-03 -0.01615733 9.288883e-01

# ensembl
# > geneTraitSignificance[rownames(de_genes),]
#                   GS.axenic  p.GS.axenic  GS.time_hr p.GS.time_hr   GS.mutant  p.GS.mutant
# Phatr3_J43365   0.286588684 1.058846e-01  0.29429533 9.642076e-02  0.02178826 0.9042043026
# Phatr3_J48554   0.004639862 9.795554e-01  0.05525707 7.600450e-01 -0.25179420 0.1574885938
# Phatr3_J46195  -0.315804463 7.339914e-02 -0.26175192 1.411650e-01 -0.15992222 0.3740023775
# Phatr3_EG00672  0.290361530 1.011677e-01  0.21071489 2.391699e-01  0.59993946 0.0002239804
# Phatr3_J48069  -0.221972493 2.144113e-01 -0.38813510 2.561388e-02 -0.29542722 0.0950868062
# Phatr3_J47845  -0.344405156 4.968284e-02 -0.31122002 7.790627e-02 -0.34136651 0.0518664368
# Phatr3_J50288  -0.056134990 7.563467e-01  0.20245458 2.585105e-01 -0.18374561 0.3060305915
# NA                       NA           NA          NA           NA          NA           NA
# Phatr3_J46647   0.307248972 8.198090e-02  0.50922555 2.473242e-03 -0.21832401 0.2222346188
# Phatr3_J55070  -0.283386541 1.100165e-01 -0.37800096 3.008903e-02 -0.12133009 0.5011985046
# Phatr3_J43073  -0.189636737 2.905096e-01 -0.25941761 1.448757e-01 -0.21297507 0.2340517159
# Phatr3_J44100  -0.152710976 3.962015e-01 -0.25503455 1.520342e-01  0.06605562 0.7149407882
# Phatr3_J48558   0.641932883 5.656588e-05  0.64487959 5.096293e-05  0.24178009 0.1752333232
 

# Yeah this shows us the same as when we did it above that the network genes are really not
# doing the same function in this network


#==============================================
#
# Intra module investigation
#
#==============================================
# Let's look at the genes of lightcyan
# that had the strongest relation to the mutant factor
module_of_interest = "lightcyan"
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

if (transcriptome == "ensembl") {
  ggsave("/home/humebc/projects/ru/nextflow_ru/run_2_nf/lightcyan_MM_GS_scatter.non_shaking.ensembl.png")
}
if (transcriptome == "NCBI") {
  ggsave("/home/humebc/projects/ru/nextflow_ru/run_2_nf/lightcyan_MM_GS_scatter.non_shaking.png")
}


# Let's have a look at the Module Membership of the 120 genes of the cyan module
geneModuleMembership[head(rownames(geneTraitSignificance %>% arrange(desc(abs(GS.mutant)))), 100),]
# Yep, all the GS for mutatnt are members of the lightcyan module.
# I want to see if there is a natural break in the GS.mutant
gs_mut = (geneTraitSignificance %>% dplyr::filter(GS.mutant>0.5))$GS.mutant
hist_df = data.frame(x=gs_mut)
ggplot(hist_df, aes(x=x)) + geom_histogram()
# There seem to be breaks around 0.5 and 0.65 but I think its probably easiest
# to plot up the network of the lightcyan genes.

module_gene_names = rownames(geneModuleMembership)[moduleGenes]
# Let's write out these names and the adjacency matrix so that we can
# load them in to make a network in a separate script
# run_2_network_viz.r

if (transcriptome == "ensembl") {
    save(geneTraitSignificance, file="/home/humebc/projects/ru/nextflow_ru/run_2_nf/geneTraitSignificance.non_shaking.ensembl.RData")
    save(module_gene_names, file="/home/humebc/projects/ru/nextflow_ru/run_2_nf/light_cyan_module_genes_names.non_shaking.ensembl.RData")
    save(adjacency, file="/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency.non_shaking.ensembl.RData")
    save.image("/home/humebc/projects/ru/nextflow_ru/run_2_nf/run_2_wgcna.ensembl.RData")
}
if (transcriptome == "NCBI") {
    save(geneTraitSignificance, file="/home/humebc/projects/ru/nextflow_ru/run_2_nf/geneTraitSignificance.non_shaking.RData")
    save(module_gene_names, file="/home/humebc/projects/ru/nextflow_ru/run_2_nf/light_cyan_module_genes_names.non_shaking.RData")
    save(adjacency, file="/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency.non_shaking.RData")
    save.image("/home/humebc/projects/ru/nextflow_ru/run_2_nf/run_2_wgcna.RData")
}


# The network is pretty cool.


# I'd like to compare the list of de genes for the mutant contrast with the
# gene list from the network (lightcyan module).
# The contrast we're intersted in has already been done
# it was called mut_vs_wt_co_de in the venn diagram.

# Get the P<0.05 DE genes for mutatant state
# NB I was first working with P<0.01 but I will now try to recreate the networks with P<0.05
# I am loading the P<0.05 object below and will set the variable to res_mutant_vs_WT_co.filtered so that I don't have to change the code
pval = 0.05
if (pval == 0.05){
  if (transcriptome == "ensembl") {
    load("/home/humebc/projects/ru/nextflow_ru/run_2_nf/run_2_res_mutant_vs_WT_co_05.filtered.ensembl.RData")
    res_mutant_vs_WT_co.filtered = res_mutant_vs_WT_co.filtered_05
}
if (transcriptome == "NCBI") {
    load("/home/humebc/projects/ru/nextflow_ru/run_2_nf/run_2_res_mutant_vs_WT_co_05.filtered.RData")
    res_mutant_vs_WT_co.filtered = res_mutant_vs_WT_co.filtered_05
}
}
if (pval == 0.01){
    if (transcriptome == "ensembl") {
        load("/home/humebc/projects/ru/nextflow_ru/run_2_nf/run_2_res_mutant_vs_WT_co.filtered.ensembl.RData")
    }
    if (transcriptome == "NCBI") {
        load("/home/humebc/projects/ru/nextflow_ru/run_2_nf/run_2_res_mutant_vs_WT_co.filtered.RData")
    }
} 




head(res_mutant_vs_WT_co.filtered)
dim(res_mutant_vs_WT_co.filtered)
sum(module_gene_names %in% rownames(res_mutant_vs_WT_co.filtered))
# > sum(module_gene_names %in% rownames(res_mutant_vs_WT_co.filtered))
# [1] 82
length(module_gene_names)
# > length(module_gene_names)
# [1] 120
# 105 out of 120 of the module genes are DE.
# Ensembl
# 77 out of 105 of the module genes are DE
res_mutant_vs_WT_co.filtered$loc = 1
res_mutant_vs_WT_co.filtered$rn = rownames(res_mutant_vs_WT_co.filtered)
ggplot(res_mutant_vs_WT_co.filtered, aes(x=loc, y=-log10(padj))) + 
geom_point(color=dplyr::case_when(res_mutant_vs_WT_co.filtered$rn %in% module_gene_names ~ "red", TRUE ~ "black"), position = "jitter") + 
ggtitle("Mutant sig genes (red== in lightcycan module membership)")

if (transcriptome == "ensembl") {
    ggsave("/home/humebc/projects/ru/nextflow_ru/run_2_nf/lightcyan.padj.scatter.non_shaking.ensembl.png")
}
if (transcriptome == "NCBI") {
    ggsave("/home/humebc/projects/ru/nextflow_ru/run_2_nf/lightcyan.padj.scatter.non_shaking.png")
}


# Then let's plot up the relationship between GS.mutant and padj of the DE genes
# There are some DE genes that are not in the network
# First we will filter these out
dim(res_mutant_vs_WT_co.filtered)
# [1] 696   8
res_mutant_vs_WT_co.filtered.nomissing = res_mutant_vs_WT_co.filtered[rownames(res_mutant_vs_WT_co.filtered) %in% rownames(geneTraitSignificance),]
# Should drop 11 genes
dim(res_mutant_vs_WT_co.filtered.nomissing)
# [1] 685   8
res_mutant_vs_WT_co.filtered.nomissing = res_mutant_vs_WT_co.filtered[rownames(res_mutant_vs_WT_co.filtered) %in% rownames(geneTraitSignificance),]

res_mutant_vs_WT_co.filtered.nomissing$GS.mutant = geneTraitSignificance[rownames(res_mutant_vs_WT_co.filtered.nomissing), "GS.mutant"]

ggplot(res_mutant_vs_WT_co.filtered.nomissing, aes(x=-log10(padj), y=abs(GS.mutant))) + geom_point() +
geom_point(color=dplyr::case_when(res_mutant_vs_WT_co.filtered.nomissing$rn %in% module_gene_names ~ "red", TRUE ~ "black"), position = "jitter") +
scale_x_continuous(trans='log10') + ggtitle("padj vs GS.mutant for the mutant DE genes (red=lightcyan module membership)")

if (transcriptome == "ensembl") {
    ggsave("/home/humebc/projects/ru/nextflow_ru/run_2_nf/padj.GSmutant.scatter.non_shaking.ensembl.png")
}
if (transcriptome == "NCBI") {
    ggsave("/home/humebc/projects/ru/nextflow_ru/run_2_nf/padj.GSmutant.scatter.non_shaking.png")
}


# We have identified one dense and highly interconnected network.
# I would be intersted to see if there are other networks in the DE genes.
# To do this I will compare the connectivity of the genes that are in the identified
# light cyan module with the connectivity scores of those that are not in the light cyan
# network to see if we can identify any additional networks.

# For some reason 11 of DE genes are missing form the geneTraitSignificance
# > rownames(res_mutant_vs_WT_co.filtered)[!rownames(res_mutant_vs_WT_co.filtered) %in% rownames(geneTraitSignificance)]
#  [1] "PHATRDRAFT_bd1354" "PHATRDRAFT_40515"  "PHATRDRAFT_41294" 
#  [4] "PHATRDRAFT_40510"  "PHATRDRAFT_31400"  "PHATRDRAFT_54954" 
#  [7] "PHATRDRAFT_40538"  "PHATRDRAFT_36287"  "PHATRDRAFT_49705" 
# [10] "PHATRDRAFT_40651"  "PHATRDRAFT_39391"
# This is OK. The difference is due to the difference between the samples used
# in making the non-shaking network vs the samples used in the Mutant vs WT co-cultured
# DE analysis. In the first the genes in quesion are filtered out in the keep step
# whereas in the second they are not.
# We will need to take this into account by removing the DE genes that are not
# found in the network (we have done this above). I.e. if genes are found in the DE genes but not in the geneTraitSignificance
# then we should remove those genes.

de_gene_names = rownames(res_mutant_vs_WT_co.filtered.nomissing)
adjacency_de = adjacency[de_gene_names, de_gene_names]

if (transcriptome == "ensembl") {
    save(adjacency_de, file="/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_de.non_shaking.ensembl.RData")
}
if (transcriptome == "NCBI") {
    save(adjacency_de, file="/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_de.non_shaking.RData")
}

adjacency_df = data.frame(adj_sum=rowSums(adjacency_de), lc_member=rownames(adjacency_de) %in% module_gene_names)
ggplot(adjacency_df, aes(x=lc_member, y=adj_sum)) + geom_violin(trim=FALSE) + geom_point() + ggtitle("adjacency_score_vs_lightcyan_membership")

if (transcriptome == "ensembl") {
    ggsave("/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_score_vs_lightcyan_membership.non_shaking.ensembl.png")
}
if (transcriptome == "NCBI") {
    ggsave("/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_score_vs_lightcyan_membership.non_shaking.png")
}


# This shows us that there are some high connectivity genes in the DE genes so there may well be some
# further networks there.

# We will use heatmaps to look for structure in the adjacency matrix

if (transcriptome == "ensembl") {
    load("/home/humebc/projects/ru/nextflow_ru/run_2_nf/light_cyan_module_genes_names.non_shaking.ensembl.RData")
}
if (transcriptome == "NCBI") {
    load("/home/humebc/projects/ru/nextflow_ru/run_2_nf/light_cyan_module_genes_names.non_shaking.RData")
}

# These are the names of the 39 genes used in the main lighcycan network.
head(module_gene_names)
# png("/home/humebc/projects/ru/nextflow_ru/run_2_nf/mutant.degenes.adjacency.heatmap.non_shaking.png", width=30, height=30, units="cm", res=600)
# heatmap(adjacency_de, ColSideColors=ifelse(rownames(adjacency_de) %in% adjacency_sub_names, "red", "black"))
# heatmap(adjacency_de, ColSideColors=ifelse(rownames(adjacency_de) %in% module_gene_names, "red", "black"), RowSideColors=ifelse(rownames(adjacency_de) %in% module_gene_names, "red", "black"))
# dev.off()
# Clearly there is some structure in there
# I have identified the light cyan module on the plot
# tree = hclust(as.dist(1-adjacency_de), method="average")
# ggdendrogram(tree, rotate = FALSE, size = 2)
# ggsave("/home/humebc/projects/ru/nextflow_ru/run_2_nf/mutant.de.non_shaking.adjacency.dendro.png", width=40, height=10, units="cm")
# plot(tree)

# To identify the other modules we should create a tree and then perform hierarchical clustering on that tree
hc <- hclust(dist(1-adjacency_de), method = "complete")
as.dendrogram(hc) %>% plot(horiz = TRUE)
abline(v=5.8, col = "red")
gene_clusters = cutree(hc, k = 20)
gene_clusters_df = as.data.frame(gene_clusters, row.names=names(gene_clusters))

# FOR RU:
# To look up the genes in a particular cluster (cluster==12) you can do e.g.:
# rownames(gene_clusters_df %>% dplyr::filter(gene_clusters==12))

# > table(gene_clusters)
# gene_clusters
#   1   2   3   4   5   6   7   8   9  10  11  12 
#  48  17  42  62  41  31  42 223  69  46  36  28

# For ensemble we get 13 clusters
# > table(gene_clusters)
# gene_clusters
#   1   2   3   4   5   6   7   8   9  10  11  12  13 
#  47  31  36  18 244  27  15  36  54  54  33  39  83

# Clustering at 5.8 should give us 12 clusters which we can then annotate on the 
# We want to add annotations according to the light cyan membership
light_cyan_membership_df = as.data.frame(ifelse(rownames(adjacency_de) %in% module_gene_names, "yes", "no"), row.names=rownames(adjacency_de), col.names=c("lc_m"))
colnames(light_cyan_membership_df) = "lc_m"
row_annotations = merge(gene_clusters_df, light_cyan_membership_df, by="row.names")
rownames(row_annotations) = row_annotations$Row.names


if (transcriptome == "ensembl") {
    row_annotations = row_annotations %>% mutate(gene_clusters=as.factor(gene_clusters)) %>% 
    # mutate(cluster_1=as.factor(gene_clusters == 1)) %>% 
    # mutate(cluster_9=as.factor(gene_clusters == 3)) %>% 
    # mutate(cluster_11_13=as.factor(gene_clusters == 9 | gene_clusters == 12 | gene_clusters == 7)) %>% 
    dplyr::select(-Row.names)
    png("/home/humebc/projects/ru/nextflow_ru/run_2_nf/mutant.de.non_shaking.adjacency.dendro.potential_networks.ensembl.png", width=40, height=30, units="cm", res=300)
}
if (transcriptome == "NCBI") {
    row_annotations = row_annotations %>% mutate(gene_clusters=as.factor(gene_clusters)) %>% 
    mutate(cluster_1=as.factor(gene_clusters == 1)) %>% 
    mutate(cluster_3=as.factor(gene_clusters == 3)) %>% 
    mutate(cluster_9_12_7=as.factor(gene_clusters == 9 | gene_clusters == 12 | gene_clusters == 7)) %>% 
    mutate(cluster_11=as.factor(gene_clusters == 11)) %>% 
    dplyr::select(-Row.names)
    png("/home/humebc/projects/ru/nextflow_ru/run_2_nf/mutant.de.non_shaking.adjacency.dendro.potential_networks.png", width=40, height=30, units="cm", res=300)
}

pheatmap(as.dist(1-adjacency_de), annotation_row=row_annotations, main="DE genes for mutant trait, clustered by adjacency matrix")
dev.off()
# We can see from the heatmap that there is considerable structure and what looks like some further discrete networks.
# I would identify 4 potential networks identified by cluster above. cluster 1 is the light cycan module genes
# then there is cluster 3 and cluster 11 that seem like discrete networks. Finally the clusters 9, 12, and 7 seem
# interconnected and I think they should be graphed as a single network.
# To plot up these neworks we will need the adjacency matrix filtered for the DE genes. THis is already saved
# (/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_de.non_shaking.RData) and we will need
# the cluster assignments which we will save now

if (transcriptome == "ensembl") {
    save(gene_clusters, file="/home/humebc/projects/ru/nextflow_ru/run_2_nf/mutant.de.non_shaking.adjacency.dendro.clusters.ensembl.RData")
}
if (transcriptome == "NCBI") {
    save(gene_clusters, file="/home/humebc/projects/ru/nextflow_ru/run_2_nf/mutant.de.non_shaking.adjacency.dendro.clusters.RData")
}


###########################

# Axenic DE

###########################

# Now that we've identified the 4 networks in the mutant DE genes
# let's have a look at the axenic de genes
# Load the DE genes from the WT ax vs WT co-cultured
if (pval == 0.05){
  if (transcriptome == "ensembl") {
    load("/home/humebc/projects/ru/nextflow_ru/run_2_nf/run_2_res_ax_vs_co_05.filtered.ensembl.RData")
    load("/home/humebc/projects/ru/nextflow_ru/run_2_nf/geneTraitSignificance.non_shaking.ensembl.RData")
    res_ax_vs_co.filtered = res_ax_vs_co.filtered_05
}
if (transcriptome == "NCBI") {
    load("/home/humebc/projects/ru/nextflow_ru/run_2_nf/run_2_res_ax_vs_co_05.filtered.RData")
    load("/home/humebc/projects/ru/nextflow_ru/run_2_nf/geneTraitSignificance.non_shaking.RData")
    res_ax_vs_co.filtered = res_ax_vs_co.filtered_05
}
}
if (pval == 0.01){
  if (transcriptome == "ensembl") {
    load("/home/humebc/projects/ru/nextflow_ru/run_2_nf/run_2_res_ax_vs_co.filtered.ensembl.RData")
    load("/home/humebc/projects/ru/nextflow_ru/run_2_nf/geneTraitSignificance.non_shaking.ensembl.RData")
    
}
if (transcriptome == "NCBI") {
    load("/home/humebc/projects/ru/nextflow_ru/run_2_nf/run_2_res_ax_vs_co.filtered.RData")
    load("/home/humebc/projects/ru/nextflow_ru/run_2_nf/geneTraitSignificance.non_shaking.RData")
    
}
}

head(res_ax_vs_co.filtered)
head(geneTraitSignificance)
dim(res_ax_vs_co.filtered)
# > dim(res_ax_vs_co.filtered)
# [1] 538   6

# Ensembl
# > dim(res_ax_vs_co.filtered)
# [1] 560   6

res_ax_vs_co.filtered.nomissing = res_ax_vs_co.filtered[rownames(res_ax_vs_co.filtered) %in% rownames(geneTraitSignificance),]
dim(res_ax_vs_co.filtered.nomissing)
# > dim(res_ax_vs_co.filtered.nomissing)
# [1] 538   6

# Ensembl
# > dim(res_ax_vs_co.filtered.nomissing)
# [1] 560   6

# Looks like there's no mismatch between the genes that were part of the network analysis and the DE analysis

de_gene_names_axenic = rownames(res_ax_vs_co.filtered.nomissing)
adjacency_de_axenic = adjacency[de_gene_names_axenic, de_gene_names_axenic]
if (transcriptome == "ensembl") {
    save(adjacency_de_axenic, file="/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_de_axenic.non_shaking.ensembl.RData")
}
if (transcriptome == "NCBI") {
    save(adjacency_de_axenic, file="/home/humebc/projects/ru/nextflow_ru/run_2_nf/adjacency_de_axenic.non_shaking.RData")
}


# Initial look at the heatmap to get an idea of the number of clusters we should be looking to identify
# What happens if we annotate the padj values
annot_df = data.frame(padj=-log10(res_ax_vs_co.filtered[rownames(adjacency_de_axenic),]$padj), abs_lg2fc=abs(res_ax_vs_co.filtered[rownames(adjacency_de_axenic),]$log2FoldChange), row.names=rownames(adjacency_de_axenic))
if (transcriptome == "ensembl") {
    png("/home/humebc/projects/ru/nextflow_ru/run_2_nf/temp_plot.ensembl.png", width=40, height=30, units="cm", res=300)
}
if (transcriptome == "NCBI") {
    png("/home/humebc/projects/ru/nextflow_ru/run_2_nf/temp_plot.png", width=40, height=30, units="cm", res=300)
}

pheatmap(as.dist(1-adjacency_de_axenic), annotation_row=annot_df)
dev.off()


# To identify the other modules we should create a tree and then perform hierarchical clustering on that tree
hc_axenic <- hclust(dist(1-adjacency_de_axenic), method = "complete")
if (transcriptome == "ensembl") {
    png("/home/humebc/projects/ru/nextflow_ru/run_2_nf/temp_plot.ensembl.png", width=40, height=30, units="cm", res=300)
}
if (transcriptome == "NCBI") {
    png("/home/humebc/projects/ru/nextflow_ru/run_2_nf/temp_plot.png", width=40, height=30, units="cm", res=300)
}
as.dendrogram(hc_axenic) %>% plot(horiz = TRUE)
abline(v=4.8, col = "red")
dev.off()
# Gives us 19 clusters
# > table(gene_clusters_axenic)
# gene_clusters_axenic
#  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 
# 19 12 22 62 20 18  8 35 32 71 33 20 30 21 16 21 30 15 17 14 22

# Ensembl
# gives us 20
# > table(gene_clusters_axenic)
# gene_clusters_axenic
#  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
# 23 17 24 47 17 55 84 33 20 20 18 22 25 30 25 30 15 17 11 27

gene_clusters_axenic = cutree(tree = hc_axenic, k = 20)
gene_clusters_df_axenic = data.frame(gene_clusters = gene_clusters_axenic, row.names=names(gene_clusters_axenic))
gene_clusters_df_axenic = gene_clusters_df_axenic %>% mutate(gene_clusters=as.factor(gene_clusters))

# row_annotations = row_annotations  %>% 
# mutate(cluster_1=as.factor(gene_clusters == 1)) %>% 
# mutate(cluster_3=as.factor(gene_clusters == 3)) %>% 
# mutate(cluster_9_12_7=as.factor(gene_clusters == 9 | gene_clusters == 12 | gene_clusters == 7)) %>% 
# mutate(cluster_11=as.factor(gene_clusters == 11))
if (transcriptome == "ensembl") {
    png("/home/humebc/projects/ru/nextflow_ru/run_2_nf/axenic.de.non_shaking.adjacency.dendro.potential_networks.ensembl.png", width=40, height=30, units="cm", res=300)
}
if (transcriptome == "NCBI") {
    png("/home/humebc/projects/ru/nextflow_ru/run_2_nf/axenic.de.non_shaking.adjacency.dendro.potential_networks.png", width=40, height=30, units="cm", res=300)
}

pheatmap(as.dist(1-adjacency_de_axenic), annotation_row=annot_df, annotation_col=gene_clusters_df_axenic, main="DE genes for axenic trait, clustered by adjacency matrix")
dev.off()

# I think there are something like 

# I'm not going to draw up the networks for the idnetified networks because the vizualisations
# are not very informative.

if (transcriptome == "ensembl") {
    save.image("/home/humebc/projects/ru/nextflow_ru/run_2_nf/run_2_wgcna.ensembl.RData")
}
if (transcriptome == "NCBI") {
    save.image("/home/humebc/projects/ru/nextflow_ru/run_2_nf/run_2_wgcna.RData")
}



# In this script we will work on making a visualization of a network
# Specifically we want to plot up the network of the interactions
# between the genes that we found to be DE using the calculated adjacency matric
# Both the objects have been saved from other analyses.
# The DE genes were calculated using run_1_nf.r
# The adjacency network was calculated using run_1_wgcna.follow.up.1.r


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
library(igraph)
library("scales")

# variable to choose which transcriptome to work with
# values can be either "NCBI" or "ensembl"
transcriptome = "ensembl"

# load in adjacency
if (transcriptome == "ensembl"){
  load("/home/humebc/projects/ru/nextflow_ru/run_1_nf/adjacency.mat.subset.ensembl.RData")
}
if (transcriptome == "NCBI"){
  load("/home/humebc/projects/ru/nextflow_ru/run_1_nf/adjacency.mat.subset.RData")
}


# load de_genes
if (transcriptome == "ensembl"){
  load("/home/humebc/projects/ru/nextflow_ru/run_1_nf/de_genes_05.subset.ensembl.RData")
}
if (transcriptome == "NCBI"){
  load("/home/humebc/projects/ru/nextflow_ru/run_1_nf/de_genes.subset.RData")
}


# The idea is to visualize the networks based on the adjacency matrix
# We can either use a hard threshold to make the weighted adjacency matrix
# unweighted, or we can aim to plot the edges with thickness.

# We also need to think about which genes to plot. Do we plot those genes with a high
# axenic gene-wise correlation, or do we plot the DE genes, or do we plot a combination,
# or do we plot the genes from the module.

# Currently we are plotting the DE expressed genes.

# To start with, let's just aim to plot up the network with some differing
# hard thresholds

# subset the adjacency matrix to only the de_genes
adjacency_sub = adjacency[rownames(de_genes_05), rownames(de_genes_05)]

# network <- graph_from_adjacency_matrix(ifelse(adjacency_sub > 0.1, 1, 0), mode="undirected")


# This is the threshold for showing a connection
if (transcriptome == "ensembl"){
  h_thresh = 0.1
}
if (transcriptome == "NCBI"){
  h_thresh = 0.2
}

adjacency_sub[adjacency_sub < h_thresh] <- 0
net <- graph_from_adjacency_matrix(adjacency_sub, mode="undirected", weighted=TRUE, diag=FALSE)
E(net)$width <- E(net)$weight * 10
V(net)$size <- 7
if (transcriptome == "ensembl"){
  V(net)$color = ifelse(names(V(net)) == "Phatr3_J43365", "red", "white")
}
if (transcriptome == "NCBI"){
  V(net)$color = ifelse(names(V(net)) == "PHATRDRAFT_43365", "red", "white")
}

plot(net)
# It would be good to have this as a dataframe so that we can plot it up using ggplot
# it would also be good to have some idea of the GS.axenic significance score
# we could use the node size to represent that.

# Produce dataframes of the vertice and edge information 
edge_df = as_data_frame(net, what="edges")
head(edge_df)
vert_df = as_data_frame(net, what="vertices")
head(vert_df)

# Load the Gene Significance df
if (transcriptome == "ensembl"){
  load("/home/humebc/projects/ru/nextflow_ru/run_1_nf/geneTraitSignificance_adj.subset.ensembl.RData")
}
if (transcriptome == "NCBI"){
  load("/home/humebc/projects/ru/nextflow_ru/run_1_nf/geneTraitSignificance_adj.subset.RData")
}
head(geneTraitSignificance_adj)

coords = layout_nicely(graph=net)
net_df_coords = data.frame(x=coords[,1], y=coords[,2], name=vert_df$name, size=vert_df$size, color=vert_df$color)
rownames(net_df_coords) = net_df_coords$name
net_df_coords$GS = geneTraitSignificance_adj[rownames(net_df_coords), "GS.axenic"]

# To do the edges we need to create a df that has the following for every edge:
# from.x, from.y, to.x, to.y, weight, width
# we will pull this information out as vectors and then combine in a dataframe for plotting
from.x = net_df_coords$x[sapply(edge_df$from, match, net_df_coords$name)]
from.y = net_df_coords$y[sapply(edge_df$from, match, net_df_coords$name)]
to.x = net_df_coords$x[sapply(edge_df$to, match, net_df_coords$name)]
to.y = net_df_coords$y[sapply(edge_df$to, match, net_df_coords$name)]
net_edge_df = data.frame(from=edge_df$from, to=edge_df$to, from.x=from.x, from.y=from.y, to.x=to.x, to.y=to.y, weight=edge_df$weight, width=edge_df$width)

# It is not possible to have a separate legend for each of the size attributes (i.e. the one used for the points and the one used for the lines)
# https://stackoverflow.com/questions/14647794/using-multiple-size-scales-in-a-ggplot

if (transcriptome == "ensembl"){
    ggplot() +
    geom_segment(data=net_edge_df,aes(x=from.x,xend = to.x, y=from.y,yend = to.y, linewidth=weight), colour="grey") +
    geom_point(data=net_df_coords, color="black", fill=net_df_coords$color, shape=21, size=rescale(net_df_coords$GS, c(10,20)), aes(x=x, y=y)) +
    geom_text(data=net_df_coords, aes(label=name, x=x, y=y)) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
            axis.text = element_blank()) + xlab("") + ylab("") +
    guides(size = guide_legend("adjacency score (weight)")) +
    ggtitle("Network of DE genes (P<0.05) according to their adjacency scores Ensembl (connections >0.1 shown)")
  ggsave("/home/humebc/projects/ru/nextflow_ru/run_1_nf/de_adj_net.subset.ensembl.png", height=50, width=40, units="cm")
}
if (transcriptome == "NCBI"){
    ggplot() + 
    geom_segment(data=net_edge_df,aes(x=from.x,xend = to.x, y=from.y,yend = to.y, linewidth=weight), colour="grey") + 
    geom_point(data=net_df_coords, color="black", fill=net_df_coords$color, shape=21, size=rescale(net_df_coords$GS, c(10,20)), aes(x=x, y=y)) + 
    geom_text(data=net_df_coords, aes(label=name, x=x, y=y)) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
            axis.text = element_blank()) + xlab("") + ylab("") +
    guides(size = guide_legend("adjacency score (weight)")) + xlim(c(-4.5, 4.5)) +
    ggtitle("Network of DE genes (P<0.01) according to their adjacency scores NCBI (connections >0.2 shown)")
  ggsave("/home/humebc/projects/ru/nextflow_ru/run_1_nf/de_adj_net.subset.png", height=50, width=40, units="cm")
}



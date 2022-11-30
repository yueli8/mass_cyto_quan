library(cytofWorkflow)
library(readxl)
library(ggplot2)
library(cowplot)
library(Seurat)
setwd("/home/hp/mass_cyto")

#read *.fcs files
rm(list = ls())
p1='paper_fcs/'
fs1=list.files(p1,'*fcs' )
fs1
fs <- read.flowSet(files = fs1,path = p1)
fs
exprs(fs[[1]])[1:6, 1:30]

#read antigens
panel <- "PBMC8_panel_v3.xlsx"
panel <- read_excel(panel)
head(data.frame(panel))

#read samples
md <- "PBMC8_metadata.xlsx"
md <- read_excel(md)
head(data.frame(md))
table(md[,2:4])

# spot check that all panel columns are in the flowSet object
all(panel$fcs_colname %in% colnames(fs))
#TRUE
# specify levels for conditions & sample IDs to assure desired ordering
md$condition <- factor(md$condition, levels = c("Ref", "BCRXL"))
md$sample_id <- factor(md$sample_id, 
                       levels = md$sample_id[order(md$condition)])

# construct SingleCellExperiment
sce <- prepData(fs, panel, md, features = panel$fcs_colname)

#Diagnostic plots
p <- plotExprs(sce, color_by = "condition")
p$facet$params$ncol <- 6
p
n_cells(sce) 
plotCounts(sce, group_by = "sample_id", color_by = "condition")
#pbMDS(sce, color_by = "condition", label_by = "sample_id")
#plotExprHeatmap(sce, scale = "last",
#               hm_pal = rev(hcl.colors(10, "YlGnBu")))
plotNRS(sce, features = "type", color_by = "condition")
set.seed(1234)
sce <- cluster(sce, features = "type",
               xdim = 10, ydim = 10, maxK = 20, seed = 1234)
#plotExprHeatmap(sce, features = "type", 
#                by = "cluster_id", k = "meta20", 
#               bars = TRUE, perc = TRUE)
plotClusterExprs(sce, k = "meta20", features = "type")

#plotMultiHeatmap(sce, 
#               hm1 = "type", hm2 = "pS6", k = "meta20", 
#              row_anno = FALSE, bars = TRUE, perc = TRUE)
# run t-SNE/UMAP on at most 500/1000 cells per sample
set.seed(1234)
sce <- runDR(sce, "TSNE", cells = 500, features = "type")
sce <- runDR(sce, "UMAP", cells = 1e3, features = "type")
plotDR(sce, "UMAP", facet_by = "sample_id")
plotDR(sce, "TSNE", facet_by = "sample_id")
p1 <- plotDR(sce, "TSNE", color_by = "meta20") + 
  theme(legend.position = "none")
p2 <- plotDR(sce, "UMAP", color_by = "meta20")
lgd <- get_legend(p2)
p2 <- p2 + theme(legend.position = "none")
plot_grid(p1, p2, lgd, nrow = 1, rel_widths = c(5, 5, 2))
# facet by sample
plotDR(sce, "UMAP", color_by = "meta20", facet_by = "sample_id")
plotCodes(sce, k = "meta20")
#plotMultiHeatmap(sce, 
 #                hm1 = "type", hm2 = "pS6", k = "som100", m = "meta20", 
  #               row_anno = FALSE, col_anno = FALSE, bars = TRUE, perc = TRUE)
for (gene in names(sce@rowRanges)) {
    plotDR(sce, "UMAP", color_by =  gene )
     ggsave2(filename = paste0(gene,'_gene_UMAP.pdf'))}


url <- "http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow"
merging_table1 <- "PBMC8_cluster_merging1.xlsx"
download.file(file.path(url, merging_table1), 
              destfile = merging_table1, mode = "wb")
merging_table1 <- read_excel(merging_table1)
head(data.frame(merging_table1))

merging_table1$new_cluster <- factor(merging_table1$new_cluster, 
                                     levels = c("B-cells IgM+", "B-cells IgM-", "CD4 T-cells",
                                                "CD8 T-cells", "DC", "NK cells", "monocytes", "surface-"))

# apply manual merging
sce <- mergeClusters(sce, k = "meta20", 
                     table = merging_table1, id = "merging1")
plotDR(sce, "UMAP", color_by = "merging1")

plotExprHeatmap(sce, features = "type", 
                by = "cluster_id", k = "meta20", m = "merging1")

plotExprHeatmap(sce, features = "type",
                by = "cluster_id", k = "merging1")

merging_table2 <- "PBMC8_cluster_merging2_v3.xlsx"
download.file(file.path(url, merging_table2), 
              destfile = merging_table2, mode = "wb")
merging_table2 <- read_excel(merging_table2)
data.frame(merging_table2)

# convert to factor with merged clusters in desired order
merging_table2$new_cluster <- factor(
  merging_table2$new_cluster, 
  levels = levels(merging_table1$new_cluster))

# apply manual merging
sce <- mergeClusters(sce, k = "meta12", 
                     table = merging_table2, id = "merging2")

plotDR(sce, "UMAP", color_by = "merging2")

# tabular comparison of algorithmic & manual merging
table(manual = cluster_codes(sce)[cluster_ids(sce), "merging2"],
      algorithm = cluster_codes(sce)[cluster_ids(sce), "meta8"] )

plot_grid(labels = c("A", "B"),
          plotDR(sce, "UMAP", color_by = "merging2"),
          plotDR(sce, "UMAP", color_by = "meta8"))

#Differential analysis
FDR_cutoff <- 0.05
plotAbundances(sce, k = "merging1", by = "sample_id")




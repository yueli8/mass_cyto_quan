library(cytofWorkflow)
library(readxl)
library(ggplot2)
library(cowplot)
library(Seurat)
library(readxl)
library(HDCytoData)
url <- "https://zenodo.org/records/10039274/files"
setwd("~/CyTOF/quan/10039274")
md <- "PBMC8_metadata.xlsx"
md <- read_excel(md)
head(data.frame(md))
fs <- Bodenmiller_BCR_XL_flowSet()
panel <- "PBMC8_panel_v3.xlsx"
panel <- read_excel(panel)
head(data.frame(panel))
all(panel$fcs_colname %in% colnames(fs))
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
pbMDS(sce, color_by = "condition", label_by = "sample_id")
plotExprHeatmap(sce, scale = "last",
                hm_pal = rev(hcl.colors(10, "YlGnBu")))
plotNRS(sce, features = "type", color_by = "condition")
#Cell population identification with FlowSOM and ConsensusClusterPlus
set.seed(1234)
sce <- cluster(sce, features = "type",
               xdim = 10, ydim = 10, maxK = 20, seed = 1234)
plotExprHeatmap(sce, features = "type", 
                by = "cluster_id", k = "meta20", 
                bars = TRUE, perc = TRUE)
plotClusterExprs(sce, k = "meta20", features = "type")
plotMultiHeatmap(sce, 
                 hm1 = "type", hm2 = "pS6", k = "meta20", 
                 row_anno = FALSE, bars = TRUE, perc = TRUE)
# run t-SNE/UMAP on at most 500/1000 cells per sample
set.seed(1234)
sce <- runDR(sce, "TSNE", cells = 500, features = "type")
sce <- runDR(sce, "UMAP", cells = 1e3, features = "type")
plotDR(sce, "UMAP", color_by = "CD4")
p1 <- plotDR(sce, "TSNE", color_by = "meta20") + 
  theme(legend.position = "none")
p2 <- plotDR(sce, "UMAP", color_by = "meta20")
lgd <- get_legend(p2)
p2 <- p2 + theme(legend.position = "none")
plot_grid(p1, p2, lgd, nrow = 1, rel_widths = c(5, 5, 2))
# facet by sample
plotDR(sce, "UMAP", color_by = "meta20", facet_by = "sample_id")
# facet by condition
plotDR(sce, "UMAP", color_by = "meta20", facet_by = "condition")
plotCodes(sce, k = "meta20")
plotMultiHeatmap(sce, 
                 hm1 = "type", hm2 = "pS6", k = "som100", m = "meta20", 
                 row_anno = FALSE, col_anno = FALSE, bars = TRUE, perc = TRUE)
merging_table1 <- "PBMC8_cluster_merging1.xlsx"
#download.file(file.path(url, merging_table1), 
    #          destfile = merging_table1, mode = "wb")
merging_table1 <- read_excel(merging_table1)
head(data.frame(merging_table1))
# convert to factor with merged clusters in desired order
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
#download.file(file.path(url, merging_table2), 
   #           destfile = merging_table2, mode = "wb")
merging_table2 <- read_excel(merging_table2)
data.frame(merging_table2)
# convert to factor with merged clusters in desired order
merging_table2$new_cluster <- factor(
  merging_table2$new_cluster, 
  levels = levels(merging_table1$new_cluster))
# apply manual merging
sce <- mergeClusters(sce, k = "meta12", 
                     table = merging_table2, id = "merging2")
# tabular comparison of algorithmic & manual merging
table(manual = cluster_codes(sce)[cluster_ids(sce), "merging2"],
      algorithm = cluster_codes(sce)[cluster_ids(sce), "meta8"] )
plot_grid(labels = c("A", "B"),
          plotDR(sce, "UMAP", color_by = "merging2"),
          plotDR(sce, "UMAP", color_by = "meta8"))
FDR_cutoff <- 0.05
plotAbundances(sce, k = "merging1", by = "sample_id")
plotAbundances(sce, k = "merging1", by = "cluster_id", shape_by = "patient_id")
ei <- metadata(sce)$experiment_info
(da_formula1 <- createFormula(ei, 
                              cols_fixed = "condition", 
                              cols_random = "sample_id"))
(da_formula2 <- createFormula(ei, 
                              cols_fixed = "condition", 
                              cols_random = c("sample_id", "patient_id")))
contrast <- createContrast(c(0, 1))
#R session aborted
#da_res1 <- diffcyt(sce, 
   #              formula = da_formula1, contrast = contrast,
     #           analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",
      #         clustering_to_use = "merging1", verbose = FALSE)


p <- plotPbExprs(sce, k = "merging1", fun = "median",
                 facet_by = "cluster_id", shape_by = "patient_id")
p$facet$params$ncol <- 2
p
ds_formula1 <- createFormula(ei, cols_fixed = "condition")
ds_formula2 <- createFormula(ei, 
                             cols_fixed = "condition", cols_random = "patient_id")
sce <- mergeClusters(sce, k = "meta20", id = "merging_all",
                     table = data.frame(old_cluster = seq_len(20), new_cluster = "all"))
p <- plotPbExprs(sce, features = "state", 
                 fun = "median",
                 shape_by = "patient_id")
p$facet$params$ncol <- 3
p


#$quan dexian
setwd("~/CyTOF")
#read *.fcs files
rm(list = ls())
p1='quan/'
fs1=list.files(p1,'*fcs' )
fs1
fs <- read.flowSet(files = fs1,path = p1)
fs
exprs(fs[[1]])[1:6, 1:30]

setwd("~/CyTOF/quan")
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
plotExprHeatmap(sce, features = "type", 
                by = "cluster_id", k = "meta20", 
               bars = TRUE, perc = TRUE)
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



library(cytofWorkflow)
library(readxl)
library(ggplot2)
library(cowplot)
library(Seurat)
library(readxl)
library(HDCytoData)

setwd("~/CyTOF/CyTOF_zhuhuiyuan")
rm(list = ls())
p1='7_wangjinxing_before_61948/'
fs1=list.files(p1,'*fcs' )
fs1
fs <- read.flowSet(files = fs1,path = p1)
fs
exprs(fs[[1]])[1:6, 1:47]

#read antigens
panel <- "antigen_zhuhuiyuan.xlsx"
panel <- read_excel(panel)
head(data.frame(panel))

#read samples
md <- "samples_zhuhuiyuan.xlsx"
md <- read_excel(md)
head(data.frame(md))
table(md[,2:4])

all(panel$fcs_colname %in% colnames(fs))
md$condition <- factor(md$condition, levels = c("Ref", "processed"))
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
#plotMultiHeatmap(sce, 
#                 hm1 = "type", hm2 = "pS6", k = "meta20", 
#                 row_anno = FALSE, bars = TRUE, perc = TRUE)
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

for (gene in names(sce@rowRanges)) {
  plotDR(sce, "UMAP", color_by =  gene )
  ggsave2(filename = paste0(gene,'_gene_UMAP.pdf'))}





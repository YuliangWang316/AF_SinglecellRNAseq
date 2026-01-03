monocle2.pipeline=function(pbmc,mk.cc)
{
library(monocle)  
library(Seurat)
  
cm.clean@meta.data$label<-Idents(cm.clean)
DimPlot(cm.clean, reduction = "umap", group.by='label', label = T)
cm.cds=as.CellDataSet(cm.clean)

cm.cds <- estimateSizeFactors(cm.cds)
cm.cds <- estimateDispersions(cm.cds)
cm.cds <- detectGenes(cm.cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cm.cds),
                                    num_cells_expressed >= 10))
pData(cm.cds)$Total_mRNAs=Matrix::colSums(exprs(cm.cds))

upper_bound <- 10^(mean(log10(pData(cm.cds)$Total_mRNAs)) +
                     2*sd(log10(pData(cm.cds)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(cm.cds)$Total_mRNAs)) -
                     2*sd(log10(pData(cm.cds)$Total_mRNAs)))
cm.cds <- cm.cds[,pData(cm.cds)$Total_mRNAs > lower_bound &
                   pData(cm.cds)$Total_mRNAs < upper_bound]

#cm.cds <- detectGenes(cm.cds, min_expr = 0.1)

ordering_genes=unique(mk.cm.0327$gene[mk.cm.0327$p_val_adj<0.05])

cm.cds <- setOrderingFilter(cm.cds, express_genes)
cm.cds <- reduceDimension(cm.cds, max_components = 2,num_dim = 6,
                          reduction_method  = 'DDRTree',residualModelFormulaStr = ~orig.ident)
cm.cds <- orderCells(cm.cds)
}

express_genes <- subset(mk.cm.0327, p_val_adj<0.05 & pct.1>0.35)$gene
tmp=grep("^MT-", express_genes, invert = TRUE)
express_genes=express_genes[tmp]
tmp=grep("^RPL", express_genes, invert = TRUE)
express_genes=express_genes[tmp]
tmp=grep("^RPS", express_genes, invert = TRUE)
express_genes=express_genes[tmp]
#trace('project2MST', edit = T, where = asNamespace("monocle"))
#and find the code:
#  if (class(projection) != "matrix") 
#    projection <- as.matrix(projection) :


my4colors <- c("#2076B5", "#2B9F2B", "#D52627", "#FF7E0C","#9366BC","#8B544A")
plot_cell_trajectory(cm.cds, color_by = "label") + scale_color_manual(values = my4colors)
plot_cell_trajectory(cm.cds, color_by = "Pseudotime")
plot_cell_trajectory(cm.cds, color_by = "orig.ident")

library(ggpubr)
# 取出cds@phenoData@data的内容
df <- pData(cm.cds)
#View(df)
ggplot(df, aes(Pseudotime, colour = label, fill = label)) +
  geom_density(bw=0.5, size=1, alpha=1)+
  facet_grid(label ~ .)+
  scale_fill_manual(values = my4colors)+
  scale_color_manual(values = my4colors)+
  theme_classic2()

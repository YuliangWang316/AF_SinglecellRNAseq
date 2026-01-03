library(Seurat)
library(patchwork)
library(dplyr)
library(tidyverse)
library(openxlsx)
library(ComplexHeatmap)
library(grid)
library(RColorBrewer)
library(circlize)

pbmc <- readRDS("/Users/yuliangwang/Desktop/merge.seu.obj.RDS")
Idents(pbmc) <- pbmc$label
new.cluster.ids <- c("Cardiomyocytes","Cardiomyocytes","Cardiomyocytes","Cardiomyocytes","Fibroblasts","Fibroblasts","Fibroblasts","Fibroblasts","Fibroblasts","Fibroblasts","Macrophages/Monocytes","Macrophages/Monocytes",
                     "Macrophages/Monocytes","Pericytes","Vascular endothelial cells","Smooth muscle cells","Lymphatic endothelial cells","Neurons","Endocardium","Lymphocytes","Adipocytes","Epicardium")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc$cluster <- Idents(pbmc)
pbmc$cluster <- factor(pbmc$cluster,levels = c("Adipocytes","Cardiomyocytes","Endocardium","Epicardium","Fibroblasts","Lymphatic endothelial cells","Lymphocytes","Macrophages/Monocytes","Neurons","Pericytes","Smooth muscle cells","Vascular endothelial cells"))
Idents(pbmc) <- pbmc$cluster
pbmc <- ScaleData(pbmc)
pbmc.markers <- FindAllMarkers(pbmc,min.pct = 0,only.pos = T)
wb <- createWorkbook()
addWorksheet(wb, "genes")

# 写入数据
writeData(wb, "genes", pbmc.markers)

# 将基因列格式设置为文本
colStyle <- createStyle(numFmt = "TEXT")
addStyle(wb, sheet = "genes", colStyle, rows = 1:(nrow(pbmc.markers)+1), cols = 1)

# 保存文件
saveWorkbook(wb, "genes_data.xlsx", overwrite = TRUE)
DotPlot(pbmc,features = c("FLT1","PTPRB","ADGRL4","MYH11","LMOD1","ACTA2","RGS5","ABCC9","PDGFRB","NRXN1","NRXN3","SORCS1","CSF1R","MRC1","MSR1",
                          "SKAP1","TRAC","THEMIS","FLT4","PROX1","CCL21","DCN","FBLN1","FBN1","PRG4","BNC1","WT1","LEPR","BMPER","NPR3",
                          "RYR2","CHRM2","KCNJ3","PLIN1","PNPLA2","PNPLA3"),
        col.min = -0.5, col.max = 2.5)  +
  
  # 使用自定义深红色色阶
  scale_color_gradient2(
    low = "#f0f0f0",
    mid = "#ff6666",
    high = "#990000",
    midpoint = 1.0,
    limits = c(-0.5, 2.5),
    breaks = c(-0.5, 0.5, 1.5, 2.5)
  ) +

# 旋转x轴标签
theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
  axis.text.y = element_text(size = 10)
)
FB <- subset(pbmc,idents = c("Fibroblasts"))
Idents(FB) <- FB$label
FB.markers <- FindAllMarkers(FB,min.pct = 0,only.pos = T)
wb <- createWorkbook()
addWorksheet(wb, "genes")

# 写入数据
writeData(wb, "genes", FB.markers)

# 将基因列格式设置为文本
colStyle <- createStyle(numFmt = "TEXT")
addStyle(wb, sheet = "genes", colStyle, rows = 1:(nrow(FB.markers)+1), cols = 1)

# 保存文件
saveWorkbook(wb, "FBgenes_data.xlsx", overwrite = TRUE)
Idents(FB) <- FB$group
FB.dff <- FindAllMarkers(FB,min.pct = 0)
wb <- createWorkbook()
addWorksheet(wb, "genes")

# 写入数据
writeData(wb, "genes", FB.dff)

# 将基因列格式设置为文本
colStyle <- createStyle(numFmt = "TEXT")
addStyle(wb, sheet = "genes", colStyle, rows = 1:(nrow(FB.dff)+1), cols = 1)

# 保存文件
saveWorkbook(wb, "FBgenesdff_VAFvs_Donor_data.xlsx", overwrite = TRUE)

CM <- subset(pbmc,idents = c("Cardiomyocytes"))
Idents(CM) <- CM$group
CM.dff <- FindAllMarkers(CM,min.pct = 0)
wb <- createWorkbook()
addWorksheet(wb, "genes")

# 写入数据
writeData(wb, "genes", CM.dff)

# 将基因列格式设置为文本
colStyle <- createStyle(numFmt = "TEXT")
addStyle(wb, sheet = "genes", colStyle, rows = 1:(nrow(CM.dff)+1), cols = 1)

# 保存文件
saveWorkbook(wb, "CMgenesdff_VAFvs_Donor_data.xlsx", overwrite = TRUE)

EC <- subset(pbmc,idents = c("Vascular endothelial cells"))
Idents(EC) <- EC$group
EC.dff <- FindAllMarkers(EC,min.pct = 0)
wb <- createWorkbook()
addWorksheet(wb, "genes")

# 写入数据
writeData(wb, "genes", EC.dff)

# 将基因列格式设置为文本
colStyle <- createStyle(numFmt = "TEXT")
addStyle(wb, sheet = "genes", colStyle, rows = 1:(nrow(EC.dff)+1), cols = 1)

# 保存文件
saveWorkbook(wb, "ECgenesdff_VAFvs_Donor_data.xlsx", overwrite = TRUE)
MO <- subset(pbmc,idents = c("Macrophages/Monocytes"))
Idents(MO) <- MO$group
MO.dff <- FindAllMarkers(MO,min.pct = 0)
wb <- createWorkbook()
addWorksheet(wb, "genes")

# 写入数据
writeData(wb, "genes", MO.dff)

# 将基因列格式设置为文本
colStyle <- createStyle(numFmt = "TEXT")
addStyle(wb, sheet = "genes", colStyle, rows = 1:(nrow(MO.dff)+1), cols = 1)

# 保存文件
saveWorkbook(wb, "MOgenesdff_VAFvs_Donor_data.xlsx", overwrite = TRUE)
FB.diff <- FB.dff %>%
  filter(cluster == "VAF")
CM.diff <- CM.dff %>%
  filter(cluster == "VAF")
EC.diff <- EC.dff %>%
  filter(cluster == "VAF")
MO.diff <- MO.dff %>%
  filter(cluster == "VAF")

diff_gene=as.data.frame(FB.diff)
rownames(diff_gene) <- diff_gene$gene
gene_list=diff_gene[,c("avg_log2FC","p_val_adj")]
colnames(gene_list)=c("logFC","padj")
gene_list$threshold = as.factor(abs(gene_list$logFC) > 0.5 & gene_list$padj < 0.05)
colored_point<-gene_list[gene_list$threshold == "TRUE",]
gene_list$threshold<-as.character(gene_list$threshold)
Spgenes<-colored_point[rownames(colored_point) == "ITGA8",]
gene_list$threshold<-as.character(gene_list$threshold)

gene_list$threshold[which(gene_list$threshold == "TRUE" & gene_list$logFC >0.5)] <- "UP"
gene_list$threshold[which(gene_list$threshold == "TRUE" & gene_list$logFC <0.5)] <- "Down"
gene_list$threshold[rownames(gene_list) == "ITGA8"]<-"LABEL"
Mycolors<-c("Blue","Gray","Black","Red")
library("ggplot2")
pdf("FB_vocano.pdf",width = 3,height = 3)

g = ggplot(data=gene_list, aes(x=logFC, y=-log10(padj),color=threshold)) + geom_point(alpha=0.6, size=0.3)  + xlim(c(-3, 3)) + ylim(c(0, 300)) +xlab("log2 fold change") + ylab("-log10 p-value") + theme_set(theme_bw()) + theme(panel.grid.major=element_line(colour=NA)) + scale_color_manual(values = Mycolors) #+ geom_text(mapping=aes(label=rownames(Spgenes)),data = Spgenes,check_overlap = TRUE,color="Black")
print(g)
dev.off()

diff_gene=as.data.frame(CM.diff)
gene_list=diff_gene[,c("avg_log2FC","p_val_adj")]
colnames(gene_list)=c("logFC","padj")
gene_list$threshold = as.factor(abs(gene_list$logFC) > 0.5 & gene_list$padj < 0.05)
colored_point<-gene_list[gene_list$threshold == "TRUE",]
gene_list$threshold<-as.character(gene_list$threshold)

gene_list$threshold[which(gene_list$threshold == "TRUE" & gene_list$logFC >0.5)] <- "UP"
gene_list$threshold[which(gene_list$threshold == "TRUE" & gene_list$logFC <0.5)] <- "Down"
Mycolors<-c("Blue","Gray","Red")
library("ggplot2")
pdf("CM_vocano.pdf",width = 3,height = 3)

g = ggplot(data=gene_list, aes(x=logFC, y=-log10(padj),color=threshold)) + geom_point(alpha=0.6, size=0.3)  + xlim(c(-3, 3)) + ylim(c(0, 300)) +xlab("log2 fold change") + ylab("-log10 p-value") + theme_set(theme_bw()) + theme(panel.grid.major=element_line(colour=NA)) +ggtitle("CMs VAF vs Donor")  + scale_color_manual(values = Mycolors,labels = c("FALSE" = "No sig"))
print(g)
dev.off()

diff_gene=as.data.frame(EC.diff)
gene_list=diff_gene[,c("avg_log2FC","p_val_adj")]
colnames(gene_list)=c("logFC","padj")
gene_list$threshold = as.factor(abs(gene_list$logFC) > 0.5 & gene_list$padj < 0.05)
colored_point<-gene_list[gene_list$threshold == "TRUE",]
gene_list$threshold<-as.character(gene_list$threshold)

gene_list$threshold[which(gene_list$threshold == "TRUE" & gene_list$logFC >0.5)] <- "UP"
gene_list$threshold[which(gene_list$threshold == "TRUE" & gene_list$logFC <0.5)] <- "Down"
Mycolors<-c("Blue","Gray","Red")
library("ggplot2")
pdf("EC_vocano.pdf",width = 3,height = 3)

g = ggplot(data=gene_list, aes(x=logFC, y=-log10(padj),color=threshold)) + geom_point(alpha=0.6, size=0.3)  + xlim(c(-3, 3)) + ylim(c(0, 300)) +xlab("log2 fold change") + ylab("-log10 p-value") + theme_set(theme_bw()) + theme(panel.grid.major=element_line(colour=NA)) +ggtitle("ECs VAF vs Donor")  + scale_color_manual(values = Mycolors,labels = c("FALSE" = "No sig"))
print(g)
dev.off()

diff_gene=as.data.frame(MO.diff)
gene_list=diff_gene[,c("avg_log2FC","p_val_adj")]
colnames(gene_list)=c("logFC","padj")
gene_list$threshold = as.factor(abs(gene_list$logFC) > 0.5 & gene_list$padj < 0.05)
colored_point<-gene_list[gene_list$threshold == "TRUE",]
gene_list$threshold<-as.character(gene_list$threshold)

gene_list$threshold[which(gene_list$threshold == "TRUE" & gene_list$logFC >0.5)] <- "UP"
gene_list$threshold[which(gene_list$threshold == "TRUE" & gene_list$logFC <0.5)] <- "Down"
Mycolors<-c("Blue","Gray","Red")
library("ggplot2")
pdf("MO_vocano.pdf",width = 3,height = 3)

g = ggplot(data=gene_list, aes(x=logFC, y=-log10(padj),color=threshold)) + geom_point(alpha=0.6, size=0.3)  + xlim(c(-3, 3)) + ylim(c(0, 300)) +xlab("log2 fold change") + ylab("-log10 p-value") + theme_set(theme_bw()) + theme(panel.grid.major=element_line(colour=NA)) + scale_color_manual(values = Mycolors)
print(g)
dev.off()

#####################
counts_matrix <- GetAssayData(CM, assay = "RNA", slot = "counts")

# 获取分组信息
groups <- CM$sample.ident
unique_groups <- unique(groups)

# 创建pseudobulk矩阵（按group求和）
pseudobulk_matrix <- matrix(0, 
                            nrow = nrow(counts_matrix), 
                            ncol = length(unique_groups),
                            dimnames = list(rownames(counts_matrix), unique_groups))

# 对每个group的细胞计数求和
for(group in unique_groups) {
  group_cells <- which(groups == group)
  if(length(group_cells) > 1) {
    pseudobulk_matrix[, group] <- rowSums(counts_matrix[, group_cells, drop = FALSE])
  } else if(length(group_cells) == 1) {
    pseudobulk_matrix[, group] <- counts_matrix[, group_cells]
  }
}

# 2. 标准化pseudobulk矩阵（CPM或TPM）
# 这里使用CPM（每百万计数）
pseudobulk_cpm <- t(t(pseudobulk_matrix) / colSums(pseudobulk_matrix)) * 1e6

# 或者使用log转换（推荐用于热图）
pseudobulk_log <- log2(pseudobulk_cpm + 1)  # 加1避免log(0)

# 3. 提取特定基因集A的表达数据
# 确保基因在矩阵中存在
available_genes <- intersect(CM.diff$gene, rownames(pseudobulk_log))
missing_genes <- setdiff(CM.diff$gene, rownames(pseudobulk_log))
if(length(missing_genes) > 0) {
  warning(paste("以下基因不在数据中:", paste(missing_genes, collapse = ", ")))
}

# 提取基因表达矩阵
gene_matrix <- pseudobulk_log[available_genes, ]

# 4. 按指定顺序排列样本
desired_order <- c("D20210826", "DLA-3152337", "Z0056", "HX-05", "HX-2", "X2Y5r05")

# 检查哪些样本存在
existing_samples <- intersect(desired_order, colnames(gene_matrix))
missing_samples <- setdiff(desired_order, colnames(gene_matrix))

if(length(missing_samples) > 0) {
  warning(paste("以下样本不存在:", paste(missing_samples, collapse = ", ")))
}

# 按指定顺序排列
gene_matrix_ordered <- gene_matrix[, existing_samples]
data <- gene_matrix_ordered
data_new <- data %>% t() %>% scale() %>% t()
colnames(data_new) <- c("Donor1","Donor2","Donor3","VAF1","VAF2","VAF3")
col_fun = colorRamp2(c(-4, 0, 4), c("#200d86", "white", "#d64339"))
htp1 <- ComplexHeatmap::Heatmap(data_new, show_row_names = F, row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                                column_names_gp = gpar(fontsize = 6), name = "Z-score", col = col_fun,
                                cluster_columns = F,cluster_rows =T,show_column_names = T)
pdf("/Users/yuliangwang/Desktop/CM_heatmap.pdf", width = 4, height = 6)
draw(htp1, 
     heatmap_legend_side = "bottom"  # 注释图例位置
)

dev.off()
###############
counts_matrix <- GetAssayData(FB, assay = "RNA", slot = "counts")

# 获取分组信息
groups <- FB$sample.ident
unique_groups <- unique(groups)

# 创建pseudobulk矩阵（按group求和）
pseudobulk_matrix <- matrix(0, 
                            nrow = nrow(counts_matrix), 
                            ncol = length(unique_groups),
                            dimnames = list(rownames(counts_matrix), unique_groups))

# 对每个group的细胞计数求和
for(group in unique_groups) {
  group_cells <- which(groups == group)
  if(length(group_cells) > 1) {
    pseudobulk_matrix[, group] <- rowSums(counts_matrix[, group_cells, drop = FALSE])
  } else if(length(group_cells) == 1) {
    pseudobulk_matrix[, group] <- counts_matrix[, group_cells]
  }
}

# 2. 标准化pseudobulk矩阵（CPM或TPM）
# 这里使用CPM（每百万计数）
pseudobulk_cpm <- t(t(pseudobulk_matrix) / colSums(pseudobulk_matrix)) * 1e6

# 或者使用log转换（推荐用于热图）
pseudobulk_log <- log2(pseudobulk_cpm + 1)  # 加1避免log(0)

# 3. 提取特定基因集A的表达数据
# 确保基因在矩阵中存在
available_genes <- intersect(FB.diff$gene, rownames(pseudobulk_log))
missing_genes <- setdiff(FB.diff$gene, rownames(pseudobulk_log))
if(length(missing_genes) > 0) {
  warning(paste("以下基因不在数据中:", paste(missing_genes, collapse = ", ")))
}

# 提取基因表达矩阵
gene_matrix <- pseudobulk_log[available_genes, ]

# 4. 按指定顺序排列样本
desired_order <- c("D20210826", "DLA-3152337", "Z0056", "HX-05", "HX-2", "X2Y5r05")

# 检查哪些样本存在
existing_samples <- intersect(desired_order, colnames(gene_matrix))
missing_samples <- setdiff(desired_order, colnames(gene_matrix))

if(length(missing_samples) > 0) {
  warning(paste("以下样本不存在:", paste(missing_samples, collapse = ", ")))
}

# 按指定顺序排列
gene_matrix_ordered <- gene_matrix[, existing_samples]
data <- gene_matrix_ordered
data_new <- data %>% t() %>% scale() %>% t()
colnames(data_new) <- c("Donor1","Donor2","Donor3","VAF1","VAF2","VAF3")
col_fun = colorRamp2(c(-4, 0, 4), c("#200d86", "white", "#d64339"))
htp1 <- ComplexHeatmap::Heatmap(data_new, show_row_names = F, row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                                column_names_gp = gpar(fontsize = 6), name = "Z-score", col = col_fun,
                                cluster_columns = F,cluster_rows =T,show_column_names = T)
pdf("/Users/yuliangwang/Desktop/FB_heatmap.pdf", width = 4, height = 6)
draw(htp1, 
     heatmap_legend_side = "bottom"  # 注释图例位置
)

dev.off()
########################################
counts_matrix <- GetAssayData(EC, assay = "RNA", slot = "counts")

# 获取分组信息
groups <- EC$sample.ident
unique_groups <- unique(groups)

# 创建pseudobulk矩阵（按group求和）
pseudobulk_matrix <- matrix(0, 
                            nrow = nrow(counts_matrix), 
                            ncol = length(unique_groups),
                            dimnames = list(rownames(counts_matrix), unique_groups))

# 对每个group的细胞计数求和
for(group in unique_groups) {
  group_cells <- which(groups == group)
  if(length(group_cells) > 1) {
    pseudobulk_matrix[, group] <- rowSums(counts_matrix[, group_cells, drop = FALSE])
  } else if(length(group_cells) == 1) {
    pseudobulk_matrix[, group] <- counts_matrix[, group_cells]
  }
}

# 2. 标准化pseudobulk矩阵（CPM或TPM）
# 这里使用CPM（每百万计数）
pseudobulk_cpm <- t(t(pseudobulk_matrix) / colSums(pseudobulk_matrix)) * 1e6

# 或者使用log转换（推荐用于热图）
pseudobulk_log <- log2(pseudobulk_cpm + 1)  # 加1避免log(0)

# 3. 提取特定基因集A的表达数据
# 确保基因在矩阵中存在
available_genes <- intersect(EC.diff$gene, rownames(pseudobulk_log))
missing_genes <- setdiff(EC.diff$gene, rownames(pseudobulk_log))
if(length(missing_genes) > 0) {
  warning(paste("以下基因不在数据中:", paste(missing_genes, collapse = ", ")))
}

# 提取基因表达矩阵
gene_matrix <- pseudobulk_log[available_genes, ]

# 4. 按指定顺序排列样本
desired_order <- c("D20210826", "DLA-3152337", "Z0056", "HX-05", "HX-2", "X2Y5r05")

# 检查哪些样本存在
existing_samples <- intersect(desired_order, colnames(gene_matrix))
missing_samples <- setdiff(desired_order, colnames(gene_matrix))

if(length(missing_samples) > 0) {
  warning(paste("以下样本不存在:", paste(missing_samples, collapse = ", ")))
}

# 按指定顺序排列
gene_matrix_ordered <- gene_matrix[, existing_samples]
data <- gene_matrix_ordered
data_new <- data %>% t() %>% scale() %>% t()
colnames(data_new) <- c("Donor1","Donor2","Donor3","VAF1","VAF2","VAF3")
col_fun = colorRamp2(c(-4, 0, 4), c("#200d86", "white", "#d64339"))
htp1 <- ComplexHeatmap::Heatmap(data_new, show_row_names = F, row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                                column_names_gp = gpar(fontsize = 6), name = "Z-score", col = col_fun,
                                cluster_columns = F,cluster_rows =T,show_column_names = T)
pdf("/Users/yuliangwang/Desktop/EC_heatmap.pdf", width = 4, height = 6)
draw(htp1, 
     heatmap_legend_side = "bottom"  # 注释图例位置
)

dev.off()
##################################
counts_matrix <- GetAssayData(MO, assay = "RNA", slot = "counts")

# 获取分组信息
groups <- MO$sample.ident
unique_groups <- unique(groups)

# 创建pseudobulk矩阵（按group求和）
pseudobulk_matrix <- matrix(0, 
                            nrow = nrow(counts_matrix), 
                            ncol = length(unique_groups),
                            dimnames = list(rownames(counts_matrix), unique_groups))

# 对每个group的细胞计数求和
for(group in unique_groups) {
  group_cells <- which(groups == group)
  if(length(group_cells) > 1) {
    pseudobulk_matrix[, group] <- rowSums(counts_matrix[, group_cells, drop = FALSE])
  } else if(length(group_cells) == 1) {
    pseudobulk_matrix[, group] <- counts_matrix[, group_cells]
  }
}

# 2. 标准化pseudobulk矩阵（CPM或TPM）
# 这里使用CPM（每百万计数）
pseudobulk_cpm <- t(t(pseudobulk_matrix) / colSums(pseudobulk_matrix)) * 1e6

# 或者使用log转换（推荐用于热图）
pseudobulk_log <- log2(pseudobulk_cpm + 1)  # 加1避免log(0)

# 3. 提取特定基因集A的表达数据
# 确保基因在矩阵中存在
available_genes <- intersect(MO.diff$gene, rownames(pseudobulk_log))
missing_genes <- setdiff(MO.diff$gene, rownames(pseudobulk_log))
if(length(missing_genes) > 0) {
  warning(paste("以下基因不在数据中:", paste(missing_genes, collapse = ", ")))
}

# 提取基因表达矩阵
gene_matrix <- pseudobulk_log[available_genes, ]

# 4. 按指定顺序排列样本
desired_order <- c("D20210826", "DLA-3152337", "Z0056", "HX-05", "HX-2", "X2Y5r05")

# 检查哪些样本存在
existing_samples <- intersect(desired_order, colnames(gene_matrix))
missing_samples <- setdiff(desired_order, colnames(gene_matrix))

if(length(missing_samples) > 0) {
  warning(paste("以下样本不存在:", paste(missing_samples, collapse = ", ")))
}

# 按指定顺序排列
gene_matrix_ordered <- gene_matrix[, existing_samples]
data <- gene_matrix_ordered
data_new <- data %>% t() %>% scale() %>% t()
colnames(data_new) <- c("Donor1","Donor2","Donor3","VAF1","VAF2","VAF3")
col_fun = colorRamp2(c(-4, 0, 4), c("#200d86", "white", "#d64339"))
htp1 <- ComplexHeatmap::Heatmap(data_new, show_row_names = F, row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                                column_names_gp = gpar(fontsize = 6), name = "Z-score", col = col_fun,
                                cluster_columns = F,cluster_rows =T,show_column_names = T)
pdf("/Users/yuliangwang/Desktop/MO_heatmap.pdf", width = 4, height = 6)
draw(htp1, 
     heatmap_legend_side = "bottom"  # 注释图例位置
)

dev.off()

###############################################
feature = c("SCN7A","COL15A1","CDH19","PLA2G2A","KAZN","C3","ELN","ITGA1","HMCN1","ATP5G3","ATP5J","RPLP1","PCOLCE2","SCARA5","SEMA3C","CCL2","RELB","TNFSF14")
FB$label <- factor(FB$label,levels = c("FB-C0-SCN7A","FB-C1-PLA2G2A","FB-C2-ELN","FB-C3-Ribosome/metabolism","FB-C4-PCOLCE2","FB-C5-CCL2"))
Idents(FB) <- FB$label
DotPlot(FB,features = c("SCN7A","COL15A1","CDH19","PLA2G2A","EYA2","NAV3","ELN","ITGA1","HMCN1","ATP5G3","ATP5J","RPLP1","PCOLCE2","SCARA5","SEMA3C","CCL2","RELB","TNFSF14"),
        col.min = -0.5, col.max = 2.5)  +
  
  # 使用自定义深红色色阶
  scale_color_gradient2(
    low = "#f0f0f0",
    mid = "#ff6666",
    high = "#990000",
    midpoint = 1.0,
    limits = c(-0.5, 2.5),
    breaks = c(-0.5, 0.5, 1.5, 2.5)
  ) +
  
  # 旋转x轴标签
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  )
Idents(FB) <- FB$label
FB_C2_ELN <- subset(FB,,idents = c("FB-C2-ELN"))
Idents(FB_C2_ELN) <- FB_C2_ELN$group
FB_C2_ELN.dff <- FindAllMarkers(FB_C2_ELN,min.pct = 0)
wb <- createWorkbook()
addWorksheet(wb, "genes")

# 写入数据
writeData(wb, "genes", FB_C2_ELN.dff)

# 将基因列格式设置为文本
colStyle <- createStyle(numFmt = "TEXT")
addStyle(wb, sheet = "genes", colStyle, rows = 1:(nrow(FB_C2_ELN.dff)+1), cols = 1)

# 保存文件
saveWorkbook(wb, "FB_C2_ELNgenesdff_VAFvs_Donor_data.xlsx", overwrite = TRUE)
FB_C2_ELN.diff <- FB_C2_ELN.dff %>%
  filter(cluster == "VAF")
diff_gene=as.data.frame(FB_C2_ELN.diff)
rownames(diff_gene) <- diff_gene$gene
gene_list=diff_gene[,c("avg_log2FC","p_val_adj")]
colnames(gene_list)=c("logFC","padj")
gene_list$threshold = as.factor(abs(gene_list$logFC) > 0.5 & gene_list$padj < 0.05)
colored_point<-gene_list[gene_list$threshold == "TRUE",]
gene_list$threshold<-as.character(gene_list$threshold)
Spgenes<-colored_point[rownames(colored_point) == "ITGA8",]
gene_list$threshold<-as.character(gene_list$threshold)

gene_list$threshold[which(gene_list$threshold == "TRUE" & gene_list$logFC >0.5)] <- "UP"
gene_list$threshold[which(gene_list$threshold == "TRUE" & gene_list$logFC <0.5)] <- "Down"
gene_list$threshold[rownames(gene_list) == "ITGA8"]<-"LABEL"
Mycolors<-c("Blue","Gray","Black","Red")
library("ggplot2")
pdf("FB_C2ELN_vocano.pdf",width = 3,height = 3)

g = ggplot(data=gene_list, aes(x=logFC, y=-log10(padj),color=threshold)) + geom_point(alpha=0.6, size=0.3)  + xlim(c(-3, 3)) + ylim(c(0, 140)) +xlab("log2 fold change") + ylab("-log10 p-value") + theme_set(theme_bw()) + theme(panel.grid.major=element_line(colour=NA)) + scale_color_manual(values = Mycolors,labels = c("FALSE" = "No sig","LABEL" = "ITGA8")) +ggtitle("FB-C2-ELN VAF vs Donor")  #+ geom_text(mapping=aes(label=rownames(Spgenes)),data = Spgenes,check_overlap = TRUE,color="Black")
print(g)
dev.off()

FB_C2_ELN.marker <- FB.markers %>%
  filter(cluster == "FB-C2-ELN")
# 提取两个数据集的基因表达数据
marker_data <- as.data.frame(FB_C2_ELN.marker)
diff_data <- as.data.frame(FB_C2_ELN.diff)

# 确保行名为基因名
rownames(marker_data) <- marker_data$gene
rownames(diff_data) <- diff_data$gene

# 找出两个数据集中共同的基因
common_genes <- intersect(rownames(marker_data), rownames(diff_data))

# 创建数据框，包含两个数据集的logFC和p_val_adj
gene_list <- data.frame(
  logFC_marker = marker_data[common_genes, "avg_log2FC"],
  padj_marker = marker_data[common_genes, "p_val_adj"],
  logFC_diff = diff_data[common_genes, "avg_log2FC"],
  padj_diff = diff_data[common_genes, "p_val_adj"],
  row.names = common_genes  # 设置行名
)

# 设置阈值：标记在两个数据集中都显著变化的基因
# 条件：两个数据集中p_val_adj < 0.05 且 logFC > 0.5
gene_list$threshold <- "No sig"

# 判断在两个数据集中都显著上调的基因（只考虑上调，不考虑下调）
up_sig <- (gene_list$padj_marker < 0.05 & gene_list$logFC_marker > 0.5) & 
  (gene_list$padj_diff < 0.05 & gene_list$logFC_diff > 0.5)
gene_list$threshold[up_sig] <- "UP"

# 标记ITGA8
if ("ITGA8" %in% rownames(gene_list)) {
  gene_list["ITGA8", "threshold"] <- "ITGA8"
  cat("ITGA8 found in common genes and will be marked\n")
} else {
  cat("Warning: ITGA8 not found in common genes\n")
}

# 设置颜色 - 移除下调的颜色
Mycolors <- c("No sig" = "Gray", "ITGA8" = "Black", "UP" = "Red")

# 创建图形
library("ggplot2")
pdf("FB_C2ELN_comparison_scatter.pdf", width = 4, height = 4)

g <- ggplot(data = gene_list, aes(x = logFC_marker, y = logFC_diff, color = threshold)) +
  geom_point(alpha = 0.6, size = 0.8) +
  xlab("log2FC (Marker)") +
  ylab("log2FC (Diff)") +
  ggtitle("FB-C2-ELN Comparison") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(colour = NA),
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 8)
  ) +
  scale_color_manual(values = Mycolors) +
  # 添加对角线参考线
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "darkgray", alpha = 0.7) +
  # 添加坐标轴参考线
  geom_hline(yintercept = 0, linetype = "solid", color = "darkgray", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "darkgray", alpha = 0.5) +
  # 添加阈值线（只画正值阈值线，因为我们只考虑上调）
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "blue", alpha = 0.5, size = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "blue", alpha = 0.5, size = 0.5)

# 添加ITGA8标签
if ("ITGA8" %in% rownames(gene_list)) {
  g <- g + 
    geom_text(data = gene_list["ITGA8", , drop = FALSE],
              aes(label = "ITGA8"),
              vjust = -1.2,
              hjust = 0.5,
              size = 3.5,
              color = "black",
              fontface = "bold")
}

print(g)
dev.off()

# 显示统计信息
cat("=================================\n")
cat("统计信息:\n")
cat("=================================\n")
cat("总共同基因数:", nrow(gene_list), "\n")
cat("显著上调基因数(UP):", sum(gene_list$threshold == "UP"), "\n")
cat("ITGA8状态:", ifelse("ITGA8" %in% rownames(gene_list), 
                       paste("Present,", gene_list["ITGA8", "threshold"]), 
                       "Not found"), "\n")

# 显示ITGA8的详细信息
if ("ITGA8" %in% rownames(gene_list)) {
  cat("\nITGA8详细信息:\n")
  cat("=================================\n")
  cat("在Marker数据集:\n")
  cat("  log2FC:", gene_list["ITGA8", "logFC_marker"], "\n")
  cat("  p_val_adj:", gene_list["ITGA8", "padj_marker"], "\n")
  cat("在Diff数据集:\n")
  cat("  log2FC:", gene_list["ITGA8", "logFC_diff"], "\n")
  cat("  p_val_adj:", gene_list["ITGA8", "padj_diff"], "\n")
}

color1 <- c(
"Adipocytes" = "#2272a9",
"Cardiomyocytes" = "#ef7c21",
"Endocardium" = "#339939",
"Epicardium" = "#ca2a29",
"Fibroblasts" = "#8b63a8",
"Lymphatic endothelial cells" = "#8a554a",
"Lymphocytes" = "#ce76ac",
"Macrophages/Monocytes" = "#7f7f7f",
"Neurons" = "#b6b72b",
"Pericytes" = "#22b4c3",
"Smooth muscle cells" = "#acc3e2",
"Vascular endothelial cells" = "#f7b778"
)

color2 <- c(
  "FB-C0-SCN7A" = "#2272a9",
  "FB-C1-PLA2G2A" = "#ef7c21",
  "FB-C2-ELN" = "#339939",
  "FB-C3-Ribosome/metabolism" = "#ca2a29",
  "FB-C4-PCOLCE2" = "#8b63a8",
  "FB-C5-CCL2" = "#8a554a"
)

color3 <- c(
  "CM1" = "#2272a9",
  "CM2" = "#339939",
  "CM3" = "#ca2a29",
  "CM4" = "#ef7c21"
)

color4 <- c(
  "Macrophage-1" = "#a74137",
  "Macrophage-2" = "#306a9c",
  "Monocytes" = "#cd7f41"
)

color5 <- c(
  "Arterial" = "#a74137",
  "Capillary" = "#306a9c",
  "Venous" = "#cd7f41"
)
Idents(pbmc) <- pbmc$orig.ident
new.cluster.ids <- c("SR-2","SR-1","VAF-3","VAF-2","VAF-1","SR-3")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc$sample_id <- Idents(pbmc)
pbmc$sample_id <- factor(pbmc$sample_id,levels = c("VAF-1","VAF-2","VAF-3","SR-1","SR-2","SR-3"))
Idents(pbmc) <- pbmc$sample_id
Cellratio <- prop.table(table(pbmc$cluster, pbmc$sample_id), margin = 2)
Cellratio <- data.frame(Cellratio)

Cellratio <- as.data.frame(Cellratio)
Cellratio$Var2<-factor(Cellratio$Var2,levels = c("VAF-1","VAF-2","VAF-3","SR-1","SR-2","SR-3"))
colnames(Cellratio)[1] <-"celltype"
ggplot(Cellratio, aes(x =Var2, y= Freq, fill = celltype)) +
  geom_col(width = 0.5, color='black')+
  # geom_flow(width=0.5,alpha=0.5, knot.pos=0.6)+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  # coord_flip()+
  scale_fill_manual(values = color1)

Idents(FB) <- FB$orig.ident
new.cluster.ids <- c("SR-2","SR-1","VAF-3","VAF-2","VAF-1","SR-3")
names(new.cluster.ids) <- levels(FB)
FB <- RenameIdents(FB, new.cluster.ids)
FB$sample_id <- Idents(FB)
FB$sample_id <- factor(FB$sample_id,levels = c("VAF-1","VAF-2","VAF-3","SR-1","SR-2","SR-3"))
Idents(FB) <- FB$sample_id
Cellratio <- prop.table(table(FB$label, FB$sample_id), margin = 2)
Cellratio <- data.frame(Cellratio)

Cellratio <- as.data.frame(Cellratio)
Cellratio$Var2<-factor(Cellratio$Var2,levels = c("VAF-1","VAF-2","VAF-3","SR-1","SR-2","SR-3"))
colnames(Cellratio)[1] <-"celltype"
ggplot(Cellratio, aes(x =Var2, y= Freq, fill = celltype)) +
  geom_col(width = 0.5, color='black')+
  # geom_flow(width=0.5,alpha=0.5, knot.pos=0.6)+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio') +
  # coord_flip()+
  scale_fill_manual(values = color2)

Idents(CM) <- CM$orig.ident
new.cluster.ids <- c("SR-2","SR-1","VAF-3","VAF-2","VAF-1","SR-3")
names(new.cluster.ids) <- levels(CM)
CM <- RenameIdents(CM, new.cluster.ids)
CM$sample_id <- Idents(CM)
CM$sample_id <- factor(CM$sample_id,levels = c("VAF-1","VAF-2","VAF-3","SR-1","SR-2","SR-3"))
Idents(CM) <- CM$sample_id
Cellratio <- prop.table(table(CM$label, CM$sample_id), margin = 2)
Cellratio <- data.frame(Cellratio)

Cellratio <- as.data.frame(Cellratio)
Cellratio$Var2<-factor(Cellratio$Var2,levels = c("VAF-1","VAF-2","VAF-3","SR-1","SR-2","SR-3"))
colnames(Cellratio)[1] <-"celltype"
ggplot(Cellratio, aes(x =Var2, y= Freq, fill = celltype)) +
  geom_col(width = 0.5, color='black')+
  # geom_flow(width=0.5,alpha=0.5, knot.pos=0.6)+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio') +
  # coord_flip()+
  scale_fill_manual(values = color3)

Idents(MO) <- MO$orig.ident
new.cluster.ids <- c("SR-2","SR-1","VAF-3","VAF-2","VAF-1","SR-3")
names(new.cluster.ids) <- levels(MO)
MO <- RenameIdents(MO, new.cluster.ids)
MO$sample_id <- Idents(MO)
MO$sample_id <- factor(MO$sample_id,levels = c("VAF-1","VAF-2","VAF-3","SR-1","SR-2","SR-3"))
Idents(MO) <- MO$sample_id
Cellratio <- prop.table(table(MO$label, MO$sample_id), margin = 2)
Cellratio <- data.frame(Cellratio)

Cellratio <- as.data.frame(Cellratio)
Cellratio$Var2<-factor(Cellratio$Var2,levels = c("VAF-1","VAF-2","VAF-3","SR-1","SR-2","SR-3"))
colnames(Cellratio)[1] <-"celltype"
ggplot(Cellratio, aes(x =Var2, y= Freq, fill = celltype)) +
  geom_col(width = 0.5, color='black')+
  # geom_flow(width=0.5,alpha=0.5, knot.pos=0.6)+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio') +
  # coord_flip()+
  scale_fill_manual(values = color4)




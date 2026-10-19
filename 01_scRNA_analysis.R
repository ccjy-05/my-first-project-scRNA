# =========================================================
# 0. 加载包
# =========================================================
library(Seurat)
library(ggplot2)
library(rhdf5)
library(dplyr)
library(tidyr)
library(patchwork)
library(Matrix)
library(scTenifoldNet)
library(scTenifoldKnk)
library(parallel)

# =========================================================
# 1. 读入数据并创建 Seurat 对象
# =========================================================
data_dir <- "E:/RSTUDIO/WORKSPACE/raw_2"

count_matrix <- Read10X(data_dir, gene.column = 1)

WT <- CreateSeuratObject(
  counts       = count_matrix,
  project      = "WT",
  min.cells    = 3,
  min.features = 200
)

# =========================================================
# 2. 质控
# =========================================================
WT[["percent.mt"]] <- PercentageFeatureSet(WT, pattern = "^mt-")
WT[["percent.hb"]] <- PercentageFeatureSet(WT, pattern = "^Hb[^(p)]")
WT[["percent.rb"]] <- PercentageFeatureSet(WT, pattern = "^Rbs|Rpl")

WT <- subset(
  WT,
  subset = nFeature_RNA > 200 &
           nFeature_RNA < 4000 &
           nCount_RNA   < 15000 &
           percent.mt   < 5 &
           percent.hb   < 5 &
           percent.rb   < 20
)

# =========================================================
# 3. 标准化、降维、聚类、UMAP
# =========================================================
sce <- NormalizeData(WT,
                     normalization.method = "LogNormalize",
                     scale.factor = 10000)

sce <- FindVariableFeatures(sce,
                            selection.method = "vst",
                            nfeatures = 2000)

all.genes <- rownames(sce)
sce <- ScaleData(sce, features = all.genes)

sce <- RunPCA(sce, features = VariableFeatures(object = sce))

sce <- FindNeighbors(sce, dims = 1:10)
sce <- FindClusters(sce, resolution = 1.0)
sce <- RunUMAP(sce, dims = 1:10)

DimPlot(sce, reduction = "umap")

# =========================================================
# 4. 细胞类型注释
# =========================================================
sce$CellType <- dplyr::recode(
  as.character(sce$seurat_clusters),
  "0"  = "Elongating Spermatids",
  "1"  = "Elongating Spermatids",
  "2"  = "Round Spermatids",
  "3"  = "Spermatozoa",
  "4"  = "Round Spermatids",
  "5"  = "Round Spermatids",
  "6"  = "Spermatogonia",
  "7"  = "Primary Spermatocytes",
  "8"  = "Spermatogonia",
  "9"  = "Secondary Spermatocytes",
  "10" = "Primary Spermatocytes"
)

# 设定细胞类型顺序
cell_order <- c(
  "Spermatogonia",
  "Primary Spermatocytes",
  "Secondary Spermatocytes",
  "Round Spermatids",
  "Elongating Spermatids",
  "Spermatozoa"
)
sce$CellType <- factor(sce$CellType, levels = cell_order)

# 查看细胞类型分布
table(sce$CellType)

# =========================================================
# 5. Marker 基因
# =========================================================
sce.markers <- FindAllMarkers(
  sce,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25
)

write.csv(sce.markers, "sce_all_markers.csv", row.names = FALSE)



  















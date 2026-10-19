# =========================================================
# 0. 加载包
# =========================================================
library(Seurat)
library(Matrix)
library(scTenifoldNet)
library(scTenifoldKnk)
library(parallel)

# =========================================================
# 1. 提取精原细胞（SPG）
# =========================================================
sc.s <- subset(sce, subset = CellType == "Spermatogonia")

# =========================================================
# 2. 取 raw counts
# =========================================================
count <- GetAssayData(sc.s, slot = "counts")

# =========================================================
# 3. 外部 QC（基于 raw counts）
# =========================================================
count_qc <- scTenifoldNet::scQC(
  X                  = count,
  minLibSize         = 500,
  removeOutlierCells = TRUE,
  minPCT             = 50 / ncol(count),   # 至少在 50 个细胞中表达的比例阈值
  maxMTratio         = 0.1
)

# =========================================================
# 4. CPM 标准化
# =========================================================
count_norm <- scTenifoldNet::cpmNormalization(count_qc)

# =========================================================
# 5. 表达覆盖过滤：至少在 50 个细胞中表达
# =========================================================
gene_detected_cells <- rowSums(count_norm > 0)
count_filt <- count_norm[gene_detected_cells >= 50, , drop = FALSE]

# =========================================================
# 6. 检查敲除靶点是否保留
# =========================================================
print("检查敲除靶点与关键基因是否保留：")
print(c("stmn1a", "tubb2b") %in% rownames(count_filt))

# =========================================================
# 7. 定义虚拟敲除函数，避免重复代码
# =========================================================
run_knk <- function(count_filt,
                    gKO,
                    nc_nNet    = 30,
                    nc_nCells  = 150,
                    seed       = 1234) {

  # 防止细胞数不足
  nc_nCells <- min(nc_nCells, ncol(count_filt))

  set.seed(seed)

  result <- scTenifoldKnk(
    countMatrix    = count_filt,
    qc             = FALSE,        # 外部已完成 QC + CPM
    gKO            = gKO,
    nc_nNet        = nc_nNet,
    nc_nCells      = nc_nCells,
    nc_nComp       = 3,
    nc_scaleScores = TRUE,
    nc_symmetric   = FALSE,
    td_K           = 3,
    td_maxIter     = 1000,
    ma_nDim        = 2,
    nCores         = max(1, parallel::detectCores() - 1)
  )

  return(result)
}

# =========================================================
# 8. 运行 stmn1a 虚拟敲除
# =========================================================
stmn1a_result <- run_knk(
  count_filt = count_filt,
  gKO        = "stmn1a",
  nc_nNet    = 30,
  nc_nCells  = 150,
  seed       = 1234
)

df_res <- stmn1a_result$diffRegulation
write.csv(df_res, "SPG_stmn1a_KO_final.csv", row.names = FALSE)

sig_genes <- subset(df_res, abs(Z) >= 1.2 & p.adj < 0.01)
sig_genes <- sig_genes[order(-abs(sig_genes$Z)), ]
head(sig_genes)

# =========================================================
# 9. 运行 tubb2b 虚拟敲除
# =========================================================
tubb2b_result <- run_knk(
  count_filt = count_filt,
  gKO        = "tubb2b",
  nc_nNet    = 30,
  nc_nCells  = 150,
  seed       = 1234
)

tubb2b_df_res <- tubb2b_result$diffRegulation
write.csv(tubb2b_df_res, "SPG_tubb2b_KO_final.csv", row.names = FALSE)

tubb2b_sig_genes <- subset(tubb2b_df_res, abs(Z) >= 1.2 & p.adj < 0.01)
tubb2b_sig_genes <- tubb2b_sig_genes[order(-abs(tubb2b_sig_genes$Z)), ]
head(tubb2b_sig_genes)

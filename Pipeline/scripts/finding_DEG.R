# 필요한 패키지 설치 및 로드
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("edgeR", quietly = TRUE))
  BiocManager::install("edgeR")

if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")

if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!requireNamespace("ggrepel", quietly = TRUE))
  install.packages("ggrepel")

library(edgeR)
library(pheatmap)
library(ggplot2)
library(ggrepel)

# 파일 읽기 함수
read_counts_file <- function(filepath) {
  df <- read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(df) <- c("geneID", filepath)  # 파일 이름을 열 이름으로 사용
  return(df)
}

# 파일 리스트 (NG 및 MG 그룹 예시)
ng_file_list <- list("OG_SRR3643545.tabular", "OG_SRR3643546.tabular", "OG_SRR3643547.tabular", "OG_SRR3643548.tabular", "OG_SRR3643549.tabular", "OG_SRR3643550.tabular", "OG_SRR3643551.tabular", "OG_SRR3643552.tabular", "OG_SRR3643553.tabular", "OG_SRR3643554.tabular", "OG_SRR3643555.tabular")
mg_file_list <- list("OG_SRR3643535.tabular", "OG_SRR3643536.tabular", "OG_SRR3643537.tabular", "OG_SRR3643538.tabular", "OG_SRR3643539.tabular", "OG_SRR3643540.tabular", "OG_SRR3643541.tabular", "OG_SRR3643542.tabular", "OG_SRR3643543.tabular", "OG_SRR3643544.tabular")

# 파일을 읽어와서 리스트에 저장
ng_count_data_list <- lapply(ng_file_list, read_counts_file)
mg_count_data_list <- lapply(mg_file_list, read_counts_file)

# geneID를 기준으로 데이터 병합
ng_merged_counts <- Reduce(function(x, y) merge(x, y, by = "geneID", all = TRUE), ng_count_data_list)
mg_merged_counts <- Reduce(function(x, y) merge(x, y, by = "geneID", all = TRUE), mg_count_data_list)

# 결측값을 0으로 대체
ng_merged_counts[is.na(ng_merged_counts)] <- 0
mg_merged_counts[is.na(mg_merged_counts)] <- 0

# 그룹 지정
group <- factor(c(rep("NG", length(ng_file_list)), rep("MG", length(mg_file_list))))

# 두 그룹 데이터 병합
merged_counts <- merge(ng_merged_counts, mg_merged_counts, by = "geneID")
rownames(merged_counts) <- merged_counts$geneID
merged_counts <- merged_counts[, -1]  # geneID 열 제거

# DGEList 객체 생성
dge <- DGEList(counts = merged_counts, group = group)

# 정규화
dge <- calcNormFactors(dge)

# 실험 디자인 행렬 생성
design <- model.matrix(~group)

# 분산 추정
dge <- estimateDisp(dge, design)

# 차등 발현 유전자 분석
fit <- glmFit(dge, design)
lrt <- glmLRT(fit)

# 결과 추출
deg_results <- topTags(lrt, n = nrow(dge))$table

# 모든 유전자에 대해 FC와 p-value 저장
all_genes_results <- deg_results[, c("logFC", "PValue")]
write.table(all_genes_results, file = "All_Genes_FC_PValue.txt", sep = "\t", row.names = TRUE, quote = FALSE)

# MDS Plot with ggplot2
mds <- plotMDS(dge, plot = FALSE)
mds_data <- data.frame(Sample = rownames(dge$samples), 
                       MDS1 = mds$x, 
                       MDS2 = mds$y, 
                       Group = group)

png(file = "MDS_plot.png", width = 2000, height = 2000, res = 300)
ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Group, label = Sample)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("NG" = "blue", "MG" = "red")) +
  theme_bw(base_size = 16) +
  ggtitle("MDS Plot") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right", legend.box = "horizontal") +
  guides(color = guide_legend(override.aes = list(size = 5)))
dev.off()

# 볼케이노 플롯
deg_results$negLogPValue <- -log10(deg_results$PValue)  # -log10(p-value) 계산

# 조건을 만족하는 데이터 포인트 색칠
deg_results$Significant <- "Not Significant"
deg_results$Significant[abs(deg_results$logFC) > 1 & deg_results$PValue < 0.05] <- "Significant"

# Significant genes 출력 및 저장
significant_genes <- deg_results[deg_results$Significant == "Significant", ]
write.table(significant_genes, file = "Significant_DEG_results.txt", sep = "\t", row.names = TRUE, quote = FALSE)

# 볼케이노 플롯 생성
png(file = "volcano_plot.png", width = 2000, height = 2000, res = 300)
ggplot(deg_results, aes(x = logFC, y = negLogPValue, color = Significant)) +
  geom_point(alpha = 0.4, size = 1.75) +
  scale_color_manual(values = c("black", "red")) +
  theme_bw(base_size = 16) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value") +
  ggtitle("Volcano Plot") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(-1, 1), col = "blue", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = "dashed")
dev.off()



suppressPackageStartupMessages({
  library(regioneR)
  library(regioneReloaded)
  library(GenomicRanges)
  library(rtracklayer)
  library(GenomeInfoDb)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(BSgenome.Hsapiens.UCSC.hg19.masked)
})

# コマンドラインで受け取る引数
args <- commandArgs(trailingOnly = TRUE)
trait_category <- args[1]
trait_name     <- args[2]
trait_id       <- args[3]
ld_r2          <- args[4]
ld_pop         <- args[5]
bed_tag        <- args[6]
CORES          <- args[7]

# パーミュテストのrep数
K_global <- 50000
# local Z scoreのrep数
K_local  <- 10000
# local Z scoreの計算する領域の長さ (bp)
LZ_window <- 1000000
# local Z scoreの計算で刻む領域の長さ (bp)
LZ_step   <- 10000

#各種データの読み込み
## 参照ゲノムとマスクの読み込み
hg     <- BSgenome.Hsapiens.UCSC.hg19
hgmask <- BSgenome.Hsapiens.UCSC.hg19.masked

## SNPデータtsvを読み込み
trait_tsv <- sprintf("example_inputSNP_%s_%s_LD_%s_%s_%s.tsv", 
                     trait_category, ld_pop, ld_r2, trait_name, trait_id)

snps <- read.table(
  trait_tsv,
  sep = "\t",
  header = FALSE,
  col.names = c("seqnames", "start", "id"),
  stringsAsFactors = FALSE
)
###SNP は 1 塩基長の領域として扱う
snps$end <- snps$start
gr_snps <- makeGRangesFromDataFrame(snps, keep.extra.columns = TRUE)

## 遺伝子領域のbedを読み込み
orf_ded <- sprintf("gene_region_ORF.bed")
gr_orf  <- rtracklayer::import(orf_bed)
### BED の name カラムに遺伝子ID (ENSG***_xx) が入っている前提
gene_col <- "name"

### 染色体の表記を統一
seqlevelsStyle(gr_snps) <- "UCSC"
seqlevelsStyle(gr_orf)  <- "UCSC"

###常染色体かつ、観測SNPがある染色体に限定
tbl      <- table(seqnames(gr_snps))
keep_chr <- names(tbl)[tbl > 0]
auto     <- paste0("chr", 1:22)
keep_chr <- intersect(keep_chr, auto)
gr_snps <- keepSeqlevels(gr_snps, keep_chr, pruning.mode = "coarse")
gr_orf  <- keepSeqlevels(gr_orf,  keep_chr, pruning.mode = "coarse")

## 染色体長を取得
seqlengths(gr_snps) <- seqlengths(hg)[seqlevels(gr_snps)]
seqlengths(gr_orf)  <- seqlengths(hg)[seqlevels(gr_orf)]

## パーミューテーションで使っていい領域をフィルタするためのマスク作成
mask <- getMask(hgmask)
mask <- keepSeqlevels(mask, keep_chr, pruning.mode = "coarse")
seqlengths(mask) <- seqlengths(hgmask)[seqlevels(mask)]

# 遺伝子領域のbedから観測SNPが存在するものを取り出す
ov <- findOverlaps(gr_snps, gr_orf, ignore.strand = TRUE)
hit_genes <- unique(as.character(mcols(gr_orf)[subjectHits(ov), gene_col]))
gr_orf_hit <- gr_orf[as.character(mcols(gr_orf)[[gene_col]]) %in% hit_genes]
orflist <- split(gr_orf_hit, mcols(gr_orf_hit)[[gene_col]])

# FDRのフィルタリングを通った遺伝子にさらに絞り込む
## フィルタリングを通った遺伝子のリスト取得
sig_file <- sprintf("significant_gene_GB_result_rep_%s_%s_LD_%s_%s_%s_.csv", K_global, ld_pop, ld_r2, trait_name, trait_id )

sig_annot <- read.csv(sig_file, stringsAsFactors = FALSE)
### そもそもフィルタリングを通った遺伝子が遺伝子がなければ終了
if (nrow(sig_annot) == 0) {
  quit(save = "no")
}
## 有望遺伝子IDだけ抽出
genes_sig_ids <- sig_annot$gene
## orflist から有望遺伝子に対応する要素だけ取り出す
orflist_sig <- orflist[genes_sig_ids]

# local Z scoreを計算
set.seed(42)
mlz <- multiLocalZscore(
  A        = gr_snps,
  Blist    = orflist_sig,
  sampling = FALSE,
  ranFUN   = "circularRandomizeRegions",
  evFUN    = "numOverlaps",
  ntimes   = K_local,
  adj_pv_method = "BH",
  genome   = "hg19",
  window   = LZ_window,
  step     = LZ_step,
  per.chromosome = TRUE,
  mask     = mask,
  mc.cores = CORES
)

cat("Completed multiLocalZscore\n")

# 結果の取得
me_lz <- getMultiEvaluation(mlz)
## 遺伝子ごとの概要
lz_resume <- me_lz$resumeTable
## 計算結果の生データ
lz_shifts <- me_lz$shifts
genes     <- names(lz_shifts)
ngene     <- length(genes)

## 代表として最初の遺伝子のシフト数 N を取得
first_v <- lz_shifts[[1]]
if (is.matrix(first_v)) {
  N <- ncol(first_v)
} else {
  N <- length(first_v)
}

## 行: shift index, 列: gene ID
mat <- matrix(NA_real_, nrow = N, ncol = ngene)
colnames(mat) <- genes

for (j in seq_along(genes)) {
  v <- lz_shifts[[j]]

  # v が matrix なら 1 行目だけ使う（通常1行しかない）
  if (is.matrix(v)) {
    z <- as.numeric(v[1, ])
  } else {
    z <- as.numeric(v)
  }

  # 長さチェック（普通は N と一致）
  if (length(z) != N) {
    warning("gene ", genes[j], ": unexpected localZ vector length ",
            length(z), " vs expected ", N, ".")
  }

  # 可能な範囲で代入
  mat[seq_len(min(N, length(z))), j] <- z
}

## shift の bp ラベル作成（実際のデータはないので、ウィンドウサイズとステップから作成）
shift_bp <- seq(-LZ_window, LZ_window, length.out = N)

## データフレーム化
lz_matrix_df <- data.frame(
  shift_bp = shift_bp,
  mat,
  check.names = FALSE  # gene ID の "." 置換を防ぐ
)



## 有望遺伝子情報（gene_name や global perm の z/FDR）とマージし出力
lz_resume2 <- merge(
  lz_resume,
  sig_annot,
  by.x = "name",  # or "RS2" 実際の列名に合わせる
  by.y = "gene",
  all.x = TRUE
)

out_lz_resume <- sprintf(
  "localZscore_resume_GB_%s_%s_LD_%s_%s_%s.csv",
  trait_category, trait_name, trait_id, ld_pop, ld_r2
)
write.csv(lz_resume2, out_lz_resume, row.names = FALSE)
cat("Save local Z score result:", out_lz_resume, "\n")

## Z scoreの計算結果の生データも出力
out_lz_matrix <- sprintf(
  "localZscore_matrix_GB_%s_%s_LD_%s_%s_%s.csv",
  trait_category, trait_name, trait_id, ld_pop, ld_r2
)

write.csv(lz_matrix_df, out_lz_matrix, row.names = FALSE)
cat("Save local Zscore matrix (shift × gene ID):", out_lz_matrix, "\n")

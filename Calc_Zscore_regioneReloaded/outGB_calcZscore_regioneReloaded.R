
suppressPackageStartupMessages({
  library(regioneR)
  library(regioneReloaded)
  library(GenomicRanges)
  library(rtracklayer)
  library(GenomeInfoDb)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(BSgenome.Hsapiens.UCSC.hg19.masked)
})

#引数設定
args <- commandArgs(trailingOnly = TRUE)

trait_category <- args[1]
trait_name     <- args[2]
trait_id       <- args[3]
ld_r2          <- args[4]
ld_pop         <- args[5]
up_frg         <- args[6]
down_frg       <- args[7]
CORES          <- args[8]


K_global  <- 50000
K_local   <- 10000
LZ_window <- 1000000
LZ_step   <- 10000

#使用するゲノムとマスクの読み込み
hg     <- BSgenome.Hsapiens.UCSC.hg19
hgmask <- BSgenome.Hsapiens.UCSC.hg19.masked

#SNP情報の読み込み
trait_tsv <- sprintf("example_inputSNP_%s_%s_LD_%s_%s_%s.tsv", 
                     trait_category, ld_pop, ld_r2, trait_name, trait_id)

snps <- read.table(
  trait_tsv,
  sep = "\t",
  header = FALSE,
  col.names = c("seqnames", "start", "id"),
  stringsAsFactors = FALSE
)
snps$end <- snps$start
gr_snps <- makeGRangesFromDataFrame(snps, keep.extra.columns = TRUE)

#遺伝子領域の読み込み
##拡大領域
exp_bed <- sprintf("gene_region_u%sk_d%sk.bed",
                   up_frg, down_frg)
gr_exp <- rtracklayer::import(exp_bed)

##ORF
orf_bed <- sprintf("gene_region_ORF.bed")
gr_orf  <- rtracklayer::import(orf_bed)

##遺伝子IDの列指定
gene_col <- "name"

#染色体の調整
##染色体の表記統一
seqlevelsStyle(gr_snps) <- "UCSC"
seqlevelsStyle(gr_exp)  <- "UCSC"
seqlevelsStyle(gr_orf)  <- "UCSC"
##常染色体かつSNPの乗っている遺伝子に限定
tbl      <- table(seqnames(gr_snps))
keep_chr <- intersect(names(tbl)[tbl > 0], paste0("chr", 1:22))
gr_snps <- keepSeqlevels(gr_snps, keep_chr, pruning.mode = "coarse")
gr_exp  <- keepSeqlevels(gr_exp,  keep_chr, pruning.mode = "coarse")
gr_orf  <- keepSeqlevels(gr_orf,  keep_chr, pruning.mode = "coarse")
##染色体の長さ取得
seqlengths(gr_snps) <- seqlengths(hg)[seqlevels(gr_snps)]
seqlengths(gr_exp)  <- seqlengths(hg)[seqlevels(gr_exp)]
seqlengths(gr_orf)  <- seqlengths(hg)[seqlevels(gr_orf)]

#マスクの作成
mask <- getMask(hgmask)
mask <- keepSeqlevels(mask, keep_chr, pruning.mode = "coarse")
seqlengths(mask) <- seqlengths(hgmask)[seqlevels(mask)]

#拡大領域が染色体の端点からはみ出す場合は切る
gr_exp <- IRanges::trim(gr_exp)
gr_orf <- IRanges::trim(gr_orf)

#遺伝子領域のbedから観測SNPが存在するものを取り出す
ov <- findOverlaps(gr_snps, gr_exp, ignore.strand = TRUE)
hit_genes <- unique(as.character(mcols(gr_exp)[subjectHits(ov), gene_col]))

gr_exp_hit <- gr_exp[as.character(mcols(gr_exp)[[gene_col]]) %in% hit_genes]
gr_orf_hit <- gr_orf[as.character(mcols(gr_orf)[[gene_col]]) %in% hit_genes]

exp_list <- split(gr_exp_hit, mcols(gr_exp_hit)[[gene_col]])
orf_list <- split(gr_orf_hit, mcols(gr_orf_hit)[[gene_col]])

exp_list <- lapply(exp_list, function(x) IRanges::trim(x))
orf_list <- lapply(orf_list, function(x) IRanges::trim(x))

#regioneRの解析と同様に拡大領域からORF領域を除いたものを遺伝子領域とする
reg_list <- lapply(names(exp_list), function(gid) {
  exp_g <- exp_list[[gid]]
  orf_g <- orf_list[[gid]]

  exp_g <- IRanges::trim(exp_g)

  if (is.null(orf_g) || length(orf_g) == 0) {
    reg_g <- exp_g
  } else {
    orf_g <- IRanges::trim(orf_g)
    reg_g <- setdiff(exp_g, orf_g)
  }

  reg_g <- IRanges::trim(reg_g)
  reg_g
})
names(reg_list) <- names(exp_list)

##もし遺伝子領域の長さが0になったら落とす
reg_list <- reg_list[lengths(reg_list) > 0]

#FDRのフィルタを通った遺伝子リストの読み込み
sig_file <- sprintf(
  "significant_gene_exp_outGB_u%sk_d%sk_result_rep_%s_%s_LD_%s_%s_%s_.csv",
  up_frg, down_frg, K_global, ld_pop, ld_r2, trait_name, trait_id
)
sig_annot <- read.csv(sig_file, stringsAsFactors = FALSE)
##該当の遺伝子がなければ落とす
if (nrow(sig_annot) == 0) {
  cat("No significant genes. Quit.\n")
  quit(save = "no")
}
genes_sig_ids <- sig_annot$gene

##multiLocalZscoreに渡せる形に変換
reg_list_sig <- reg_list[genes_sig_ids]
reg_list_sig <- reg_list_sig[lengths(reg_list_sig) > 0]
if (length(reg_list_sig) == 0) {
  cat("No regulatory-only regions remain after filtering. Quit.\n")
  quit(save = "no")
}

#local Z scoreを計算
set.seed(42)
mlz <- multiLocalZscore(
  A        = gr_snps,
  Blist    = reg_list_sig,
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

cat("Completed multiLocalZscore \n")

#結果の取得
me_lz <- getMultiEvaluation(mlz)
##概要の方
lz_resume <- me_lz$resumeTable
##計算結果の生データ
lz_shifts <- me_lz$shifts

genes <- names(lz_shifts)
ngene <- length(genes)

# 代表として最初の遺伝子のシフト数 N を取得
first_v <- lz_shifts[[1]]
N <- if (is.matrix(first_v)) ncol(first_v) else length(first_v)

# 行: shift index, 列: gene ID
mat <- matrix(NA_real_, nrow = N, ncol = ngene)
colnames(mat) <- genes

for (j in seq_along(genes)) {
  v <- lz_shifts[[j]]
  z <- if (is.matrix(v)) as.numeric(v[1, ]) else as.numeric(v)
  mat[seq_len(min(N, length(z))), j] <- z
}

#シフト幅の作成
shift_bp <- seq(-LZ_window, LZ_window, length.out = N)

#データフレーム化
lz_matrix_df <- data.frame(
  shift_bp = shift_bp,
  mat,
  check.names = FALSE
)

##概要の方はフィルタ通った方の遺伝子の結果とマージ
lz_resume2 <- merge(
  lz_resume,
  sig_annot,
  by.x = "name",
  by.y = "gene",
  all.x = TRUE
)

#出力
out_lz_resume <- sprintf(
  "localZscore_resume_outGB_%s_%s_%s_%s_%s_LD_%s_%s.csv",
  up_frg, down_frg, trait_category, trait_name, trait_id, ld_pop, ld_r2
)
write.csv(lz_resume2, out_lz_resume, row.names = FALSE)
cat("Save local Z score result:", out_lz_resume, "\n")

out_lz_matrix <- sprintf(
  "localZscore_matrix_outGB_%s_%s_%s_%s_%s_LD_%s_%s.csv",
  up_frg, down_frg, trait_category, trait_name, trait_id, ld_pop, ld_r2
)
write.csv(lz_matrix_df, out_lz_matrix, row.names = FALSE)
cat("Save local Zscore matrix (shift × gene ID):", out_lz_matrix, "\n")


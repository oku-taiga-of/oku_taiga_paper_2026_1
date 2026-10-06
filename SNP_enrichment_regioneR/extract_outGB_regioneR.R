#!/usr/bin/env Rscript

## loading packages
suppressPackageStartupMessages({
  library(regioneR)
  library(parallel)
  library(GenomicRanges)
  library(rtracklayer)
  library(GenomeInfoDb)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(BSgenome.Hsapiens.UCSC.hg19.masked)
  library(S4Vectors)
  library(IRanges)
})

## 引数参照
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
  stop(
    paste0(
      "Usage:\n",
      "Rscript extract_outGB_regioneR.R ",
      "<trait_category> <trait_name> <trait_id> <ld_r2> <ld_pop> <bed_tag> <CORES> [out_dir]\n",
      "\nExample:\n",
      "Rscript extract_outGB_regioneR.R Power APLM 0004980 0.1 EUR up50_down10 4 results_exp_outGB\n"
    )
  )
}

trait_category <- args[1]
trait_name     <- args[2]
trait_id       <- args[3]
ld_r2          <- args[4]
ld_pop         <- args[5]
bed_tag        <- args[6]

## 使用するコア数
CORES <- as.integer(args[7])
if (is.na(CORES) || CORES < 1L) {
  stop("CORES must be a positive integer.")
}

## 出力先
out_dir <- if (length(args) >= 8) args[8] else "."
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## 入力SNP tsv のフルパス（再解析指示_02: 窓修正後の {cat}_{window}/ を指す）。省略時は従来パス
input_tsv_arg <- if (length(args) >= 9 && nzchar(args[9])) args[9] else ""
## 有意判定の FDR 閾値（再解析指示_03 v5: FDR_BH のみ・p_emp 条件は除去）。既定 0.05
sig_fdr <- if (length(args) >= 10 && nzchar(args[10])) as.numeric(args[10]) else 0.05

## ------------------------------
## 遺伝子領域 BED の置き場所
## ------------------------------
## 解決順: 第11引数 → 環境変数 GENE_BED_DIR → このスクリプトと同じディレクトリ。
## 公開パッケージでは BED をスクリプトと同じ場所に置けばそのまま動く。
script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  m  <- grep("^--file=", ca, value = TRUE)
  if (length(m)) normalizePath(dirname(sub("^--file=", "", m[1]))) else getwd()
}
bed_dir <- if (length(args) >= 11 && nzchar(args[11])) args[11] else Sys.getenv("GENE_BED_DIR", "")
if (!nzchar(bed_dir)) bed_dir <- script_dir()
if (!dir.exists(bed_dir)) stop(paste("BED directory not found:", bed_dir))

## 参照ゲノムデータ
hg     <- BSgenome.Hsapiens.UCSC.hg19
hgmask <- BSgenome.Hsapiens.UCSC.hg19.masked

## 設定
## K（permutation 回数）は原則 50000。第12引数で下げられるが、これは**動作確認用**。
## 出力ファイル名に K が入るため、K を変えた結果は本番の結果と取り違えようがない。
K <- 50000
if (length(args) >= 12 && nzchar(args[12])) {
  k_arg <- suppressWarnings(as.integer(args[12]))
  if (is.na(k_arg) || k_arg < 1L) stop("K must be a positive integer.")
  if (k_arg != 50000L) {
    message(sprintf(
      "[NOTE] K = %d (not the 50000 used in the manuscript). Use this only for a quick trial run; the value of K appears in the output file names.",
      k_arg))
  }
  K <- k_arg
}
genome_ref <- hgmask
use_mask <- TRUE

## ------------------------------
## 補助関数
## ------------------------------

## BED 4列目が "ENSG000001234|GENE_NAME" 形式であることを想定し、
## 検定単位には gene_id、出力注釈には gene_name を使う。
## gene_name がBEDにない場合は空欄にする。
add_gene_id_name_from_bed_name <- function(gr, gene_col = "name") {
  if (!gene_col %in% names(mcols(gr))) {
    stop("BEDに name 列がありません。BED 4列目に gene_id または gene_id|gene_name が必要です。")
  }

  raw_name <- as.character(mcols(gr)[[gene_col]])

  gene_id <- sub("\\|.*$", "", raw_name)

  gene_name <- ifelse(
    grepl("\\|", raw_name),
    sub("^[^|]*\\|", "", raw_name),
    ""
  )

  ## "ENSG|" のようなケースは空欄扱い
  gene_name[is.na(gene_name)] <- ""

  mcols(gr)$gene_id   <- gene_id
  mcols(gr)$gene_name <- gene_name

  gr
}

make_gene_name_map <- function(gr) {
  df <- data.frame(
    gene_id   = as.character(mcols(gr)$gene_id),
    gene_name = as.character(mcols(gr)$gene_name),
    stringsAsFactors = FALSE
  )

  df <- unique(df)

  aggregate(
    gene_name ~ gene_id,
    data = df,
    FUN = function(x) paste(unique(x[x != "" & !is.na(x)]), collapse = "/")
  )
}

## ------------------------------
## 入力ファイル
## ------------------------------

## TSV（SNPデータ）をGRangesへ変換
## 注意: このTSVの start は VCF由来の 1-based POS を想定する
trait_tsv <- if (nzchar(input_tsv_arg)) {
  input_tsv_arg
} else {
  ## 第9引数が無いときは、命名規約からファイル名を組み立てて
  ## 環境変数 SNP_TSV_DIR（未設定ならカレントディレクトリ）の下を探す。
  file.path(
    Sys.getenv("SNP_TSV_DIR", "."),
    sprintf("vcf_%s_LD_r%s_%s_%s_withID.tsv", ld_pop, ld_r2, trait_name, trait_id)
  )
}

if (!file.exists(trait_tsv)) {
  stop(paste("SNP TSV not found:", trait_tsv))
}

snps <- read.table(
  trait_tsv,
  sep = "\t",
  header = FALSE,
  col.names = c("seqnames", "start", "id")
)

snps$end <- snps$start
gr_snps <- makeGRangesFromDataFrame(snps, keep.extra.columns = TRUE)

## 拡大遺伝子領域のBED読み込み
## GB版と同じく bed_tag を使って protein_genes_{bed_tag}.bed を読む
exp_bed <- file.path(bed_dir, sprintf("protein_genes_%s.bed", bed_tag))

if (!file.exists(exp_bed)) {
  stop(paste0("Expanded BED not found: ", exp_bed,
              "\n  (BED の置き場所は第11引数 / 環境変数 GENE_BED_DIR / スクリプトと同じディレクトリ の順に解決します)"))
}

gr_exp <- rtracklayer::import(exp_bed)

## gene body領域のBED読み込み
## 拡大領域から差し引く基準として、常に up0_down0 を使う
gb_bed <- file.path(bed_dir, "protein_genes_up0_down0.bed")

if (!file.exists(gb_bed)) {
  stop(paste0("Gene body BED not found: ", gb_bed,
              "\n  (outGB は拡大領域から gene body を差し引くため、常に up0_down0 の BED が要ります)"))
}

gr_gb <- rtracklayer::import(gb_bed)

## BED 4列目から gene_id / gene_name を作る
gr_exp <- add_gene_id_name_from_bed_name(gr_exp, gene_col = "name")
gr_gb  <- add_gene_id_name_from_bed_name(gr_gb,  gene_col = "name")

## 出力用gene_nameは拡大領域BED由来を使う
gene_name_map <- make_gene_name_map(gr_exp)

## ------------------------------
## 染色体処理
## ------------------------------

## 染色体表記を統一
seqlevelsStyle(gr_snps) <- "UCSC"
seqlevelsStyle(gr_exp)  <- "UCSC"
seqlevelsStyle(gr_gb)   <- "UCSC"

## 常染色体かつSNPが乗っている常染色体に限定
auto <- paste0("chr", 1:22)
keep_chr <- intersect(names(table(seqnames(gr_snps)))[table(seqnames(gr_snps)) > 0], auto)

gr_snps <- keepSeqlevels(gr_snps, keep_chr, pruning.mode = "coarse")
gr_exp  <- keepSeqlevels(gr_exp,  keep_chr, pruning.mode = "coarse")
gr_gb   <- keepSeqlevels(gr_gb,   keep_chr, pruning.mode = "coarse")

## 染色体長を取得
seqlengths(gr_snps) <- seqlengths(hg)[seqlevels(gr_snps)]
seqlengths(gr_exp)  <- seqlengths(hg)[seqlevels(gr_exp)]
seqlengths(gr_gb)   <- seqlengths(hg)[seqlevels(gr_gb)]

## はみ出しを全てトリム
gr_snps <- IRanges::trim(gr_snps)
gr_snps <- gr_snps[width(gr_snps) > 0]

gr_exp <- IRanges::trim(gr_exp)
gr_exp <- gr_exp[width(gr_exp) > 0]

gr_gb <- IRanges::trim(gr_gb)
gr_gb <- gr_gb[width(gr_gb) > 0]

## 使用可能な領域のマスクを作成
mask <- getMask(hgmask)
seqlevelsStyle(mask) <- "UCSC"
mask <- keepSeqlevels(mask, keep_chr, pruning.mode = "coarse")
seqlengths(mask) <- seqlengths(hgmask)[seqlevels(mask)]
mask <- IRanges::trim(mask)
mask <- mask[width(mask) > 0]

## ------------------------------
## 観測値の取得
## ------------------------------

gene_col <- "gene_id"

## 拡大領域とgene body領域の観測SNP数をそれぞれ観測
ovF <- findOverlaps(gr_snps, gr_exp, ignore.strand = TRUE)
ovG <- findOverlaps(gr_snps, gr_gb,  ignore.strand = TRUE)

tabF <- table(as.character(mcols(gr_exp)[subjectHits(ovF), gene_col]))
tabG <- table(as.character(mcols(gr_gb)[subjectHits(ovG),  gene_col]))

## data.frameにして geneでマージ（名前ベースで明示的に揃える）
dfF <- data.frame(
  gene = names(tabF),
  full = as.integer(tabF),
  stringsAsFactors = FALSE
)

dfG <- data.frame(
  gene = names(tabG),
  gb   = as.integer(tabG),
  stringsAsFactors = FALSE
)

df_counts <- merge(dfF, dfG, by = "gene", all = TRUE)

## NAを0に
df_counts$full[is.na(df_counts$full)] <- 0L
df_counts$gb[is.na(df_counts$gb)]     <- 0L

## 拡大領域のみの観測数 = FULL − gene body
## 負値は念のため0に丸める
df_counts$outside_gb <- df_counts$full - df_counts$gb
df_counts$outside_gb[df_counts$outside_gb < 0L] <- 0L

## 名前付きベクトルとして取り出す
obs_outside_all <- df_counts$outside_gb
names(obs_outside_all) <- df_counts$gene

## 拡大領域のみの観測SNP数が1以上の遺伝子に限定
## これは計算量削減と、上側検定でobs=0の遺伝子はenrichment候補にならないため
min_obs <- 1L
idx_hit <- obs_outside_all >= min_obs

genes <- names(obs_outside_all)[idx_hit]
obs_vec <- obs_outside_all[genes]

if (length(genes) == 0L) {
  stop("観測でヒットした遺伝子がありません（expanded outside gene body）。min_obs を見直してください。")
}

## gr_expもヒット遺伝子にサブセット
sel <- as.character(mcols(gr_exp)[[gene_col]]) %in% genes
gr_exp_hit <- gr_exp[sel]
gr_exp_hit <- keepSeqlevels(gr_exp_hit, keep_chr, pruning.mode = "coarse")
seqlengths(gr_exp_hit) <- seqlengths(hg)[seqlevels(gr_exp_hit)]
gr_exp_hit <- IRanges::trim(gr_exp_hit)
gr_exp_hit <- gr_exp_hit[width(gr_exp_hit) > 0]

## gene body側もgene_idで使うため、そのまま保持
gr_gb <- keepSeqlevels(gr_gb, keep_chr, pruning.mode = "coarse")
seqlengths(gr_gb) <- seqlengths(hg)[seqlevels(gr_gb)]
gr_gb <- IRanges::trim(gr_gb)
gr_gb <- gr_gb[width(gr_gb) > 0]

## ------------------------------
## permutation
## ------------------------------

## ここから並列Permutation（巨大行列は作らない）
RNGkind("L'Ecuyer-CMRG")
set.seed(42)

## 集計用ベクトル（遺伝子ごと）
sum_x  <- numeric(length(genes))
sum_x2 <- numeric(length(genes))
ge_cnt <- integer(length(genes))

names(sum_x) <- names(sum_x2) <- names(ge_cnt) <- genes

## 1回分のPermutation → 遺伝子ごとのカウントを返す関数
perm_once <- function(i) {
  rp <- circularRandomizeRegions(
    gr_snps,
    genome = genome_ref,
    per.chromosome = TRUE,
    mask = if (use_mask) mask else NULL
  )

  ## FULL側
  ovpF <- findOverlaps(rp, gr_exp_hit, ignore.strand = TRUE)
  tF <- table(as.character(mcols(gr_exp_hit)[subjectHits(ovpF), gene_col]))

  ## gene body側
  ovpG <- findOverlaps(rp, gr_gb, ignore.strand = TRUE)
  tG <- table(as.character(mcols(gr_gb)[subjectHits(ovpG), gene_col]))

  ## outside-gene-bodyカウント用ベクトル
  v_out <- integer(length(genes))
  names(v_out) <- genes

  ## FULLカウントを追加（genesに含まれるものだけ）
  if (length(tF)) {
    gF <- intersect(names(tF), genes)
    if (length(gF)) {
      v_out[gF] <- v_out[gF] + as.integer(tF[gF])
    }
  }

  ## gene bodyカウントを減算（genesに含まれるものだけ）
  if (length(tG)) {
    gG <- intersect(names(tG), genes)
    if (length(gG)) {
      v_out[gG] <- v_out[gG] - as.integer(tG[gG])
    }
  }

  ## 念のためマイナスは0に丸める
  v_out[v_out < 0L] <- 0L

  v_out
}

## 並列でK回走らせ、逐次で合算
CHUNK <- 5000
nchunk <- ceiling(K / CHUNK)

for (b in seq_len(nchunk)) {
  this_B <- if (b < nchunk) CHUNK else (K - CHUNK * (nchunk - 1))
  idxs <- (1:this_B) + (b - 1) * CHUNK

  res_list <- mclapply(idxs, perm_once, mc.cores = CORES)

  for (v_out in res_list) {
    sum_x  <- sum_x  + v_out
    sum_x2 <- sum_x2 + v_out * v_out
    ge_cnt <- ge_cnt + as.integer(v_out >= obs_vec)
  }

  cat(sprintf("chunk %d/%d done (%s perms)\n",
              b, nchunk, format(this_B, big.mark = ",")))
}

## ------------------------------
## 統計量
## ------------------------------

mean_hat <- sum_x / K
var_hat <- (sum_x2 - (sum_x^2) / K) / pmax(1, K - 1)
var_hat[var_hat < 0] <- 0
sd_hat <- sqrt(var_hat)

p_emp <- (ge_cnt + 1) / (K + 1)
z_raw <- ifelse(sd_hat > 0, (obs_vec - mean_hat) / sd_hat, NA_real_)

## 多重検定補正
p_vec <- p_emp
ok <- is.finite(p_vec) & p_vec >= 0 & p_vec <= 1
p_in <- p_vec[ok]

q_bh <- rep(NA_real_, length(p_vec))
p_holm <- rep(NA_real_, length(p_vec))

q_bh[ok] <- p.adjust(p_in, method = "BH")
p_holm[ok] <- p.adjust(p_in, method = "holm")

names(q_bh) <- names(p_vec)
names(p_holm) <- names(p_vec)

## 結果データフレーム
res_gene <- data.frame(
  gene   = genes,
  obs    = as.integer(obs_vec[genes]),
  mean   = mean_hat[genes],
  sd     = sd_hat[genes],
  z      = z_raw[genes],
  p_emp  = p_emp[genes],
  FDR_BH = q_bh[genes],
  HOLM   = p_holm[genes],
  row.names = NULL
)

## BED由来のgene_nameを付与
res_gene$.__ord <- seq_len(nrow(res_gene))
res_gene <- merge(
  res_gene,
  gene_name_map,
  by.x = "gene",
  by.y = "gene_id",
  all.x = TRUE,
  sort = FALSE
)
res_gene <- res_gene[order(res_gene$.__ord), ]
res_gene$.__ord <- NULL
res_gene$gene_name[is.na(res_gene$gene_name)] <- ""

res_gene <- res_gene[, c("gene", "gene_name",
                         setdiff(names(res_gene), c("gene", "gene_name")))]

## ログ
cat(sprintf("# family size tested hit genes m = %d\n", length(genes)))

## 全遺伝子の結果出力
res_gene <- res_gene[order(res_gene$p_emp, -res_gene$z, res_gene$gene), ]

out_csv <- file.path(
  out_dir,
  sprintf(
    "all_gene_exp_outGB_%s_result_rep_%s_%s_LD_%s_%s_%s.csv",
    bed_tag, K, ld_pop, ld_r2, trait_name, trait_id
  )
)

write.csv(res_gene, out_csv, row.names = FALSE)
cat("保存しました（検定対象遺伝子）:", out_csv, "\n")

## 可能性のある遺伝子の出力
## 再解析指示_03 v5: 有意判定は FDR_BH のみ（p_emp≤0.01 の二重条件は除去確定）
thr_obs <- 1
hits <- subset(res_gene, obs >= thr_obs & FDR_BH <= sig_fdr)
hits <- hits[order(hits$p_emp, -hits$z, hits$gene), ]

outfile <- file.path(
  out_dir,
  sprintf(
    "significant_gene_exp_outGB_%s_result_rep_%s_%s_LD_%s_%s_%s.csv",
    bed_tag, K, ld_pop, ld_r2, trait_name, trait_id
  )
)

write.csv(hits, outfile, row.names = FALSE)
cat(sprintf("保存しました(条件を満たした遺伝子): %s（%d genes）\n", outfile, nrow(hits)))

cat("[OK] Finished\n")
cat("[IN] trait_tsv:", trait_tsv, "\n")
cat("[IN] expanded_bed:", exp_bed, "\n")
cat("[IN] gene_body_bed:", gb_bed, "\n")
cat("[IN] bed_tag:", bed_tag, "\n")
cat("[IN] K:", K, "\n")
cat("[IN] CORES:", CORES, "\n")
cat("[IN] out_dir:", out_dir, "\n")
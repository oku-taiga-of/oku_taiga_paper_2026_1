

##loading packages
suppressPackageStartupMessages({
  library(regioneR)
  library(parallel)
  library(GenomicRanges)
  library(rtracklayer)
  library(GenomeInfoDb)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(BSgenome.Hsapiens.UCSC.hg19.masked)
  library(S4Vectors) 
})

##引数参照
args <- commandArgs(trailingOnly = TRUE)
trait_category <- args[1]
trait_name <- args[2]
trait_id <- args[3]
ld_r2 <- args[4]
ld_pop <- args[5]
up_frg <- args[6]
down_frg <- args[7]
CORES <- as.integer(args[8])


##参照ゲノムデータ
hg     <- BSgenome.Hsapiens.UCSC.hg19 
hgmask <- BSgenome.Hsapiens.UCSC.hg19.masked

#TSV（SNPデータ）をGRangesへ変換
trait_tsv <- sprintf("example_inputSNP_%s_%s_LD_%s_%s_%s.tsv", 
                     trait_category, ld_pop, ld_r2, trait_name, trait_id)
snps <- read.table(trait_tsv, sep="\t", header=FALSE,
                   col.names=c("seqnames","start","id"))
snps$end <- snps$start
gr_snps <- makeGRangesFromDataFrame(snps, keep.extra.columns=TRUE)

#拡大遺伝子領域のbed読み込み
exp_ded <- sprintf("gene_region_u%sk_d%sk.bed",
                   up_frg, down_frg)
gr_exp    <- rtracklayer::import(exp_ded)
#ORFのみの領域のbed読み込み
orf_bed <- "gene_region_ORF.bed" 
gr_orf <- rtracklayer::import(orf_bed)


#染色体処理
##染色体表記を統一
seqlevelsStyle(gr_snps)   <- "UCSC"
seqlevelsStyle(gr_exp)    <- "UCSC"
seqlevelsStyle(gr_orf)    <- "UCSC"

##常染色体かつSNPが乗っている常染色体に限定
auto <- paste0("chr", 1:22)
keep_chr <- intersect(names(table(seqnames(gr_snps)))[table(seqnames(gr_snps)) > 0], auto)
gr_snps <- keepSeqlevels(gr_snps, keep_chr, pruning.mode="coarse")
gr_exp  <- keepSeqlevels(gr_exp,  keep_chr, pruning.mode="coarse")
gr_orf <- keepSeqlevels(gr_orf, keep_chr, pruning.mode="coarse")

##マスク作成
mask <- getMask(hgmask)
seqlevelsStyle(mask) <- "UCSC"
mask <- keepSeqlevels(mask, keep_chr, pruning.mode="coarse")
## 染色体長を取得
seqlengths(gr_snps) <- seqlengths(hg)[seqlevels(gr_snps)]
seqlengths(gr_exp)  <- seqlengths(hg)[seqlevels(gr_exp)]
seqlengths(gr_orf) <- seqlengths(hg)[seqlevels(gr_orf)]
seqlengths(mask)    <- seqlengths(hgmask)[seqlevels(mask)]
##はみ出しを全てトリム
gr_snps <- IRanges::trim(gr_snps); gr_snps <- gr_snps[width(gr_snps) > 0]
gr_exp  <- IRanges::trim(gr_exp);  gr_exp  <- gr_exp [width(gr_exp ) > 0]
gr_orf <- IRanges::trim(gr_orf); gr_orf <- gr_orf[width(gr_orf) > 0]
mask    <- IRanges::trim(mask);    mask    <- mask   [width(mask   ) > 0]

#検定
##設定
K <- 50000
genome_ref <- hgmask
use_mask <- TRUE

##遺伝子IDの列名を取得
gene_col <- "name"
mcols(gr_orf)$name <- as.character(mcols(gr_orf)$name)
mcols(gr_exp)$name <- as.character(mcols(gr_exp)$name)

##ORFと拡大領域の観測SNP数をそれぞれ観測
ovF  <- findOverlaps(gr_snps, gr_exp, ignore.strand=TRUE)
ovO  <- findOverlaps(gr_snps, gr_orf, ignore.strand=TRUE)

tabF <- table(as.character(mcols(gr_exp)[subjectHits(ovF), gene_col]))
tabO <- table(as.character(mcols(gr_orf)[subjectHits(ovO), gene_col]))

## data.frameにして geneでマージ（名前ベースで明示的に揃える）
dfF <- data.frame(
  gene = names(tabF),
  full = as.integer(tabF),
  stringsAsFactors = FALSE
)
dfO <- data.frame(
  gene = names(tabO),
  orf  = as.integer(tabO),
  stringsAsFactors = FALSE
)

df_counts <- merge(dfF, dfO, by = "gene", all = TRUE)

## NA を 0 に
df_counts$full[is.na(df_counts$full)] <- 0L
df_counts$orf [is.na(df_counts$orf)]  <- 0L

## フランキング観測数 = FULL − ORF（負値は0に丸める）
df_counts$flank <- df_counts$full - df_counts$orf
df_counts$flank[df_counts$flank < 0L] <- 0L

## 名前付きベクトルとして取り出す
obs_flank_all <- df_counts$flank
names(obs_flank_all) <- df_counts$gene


##拡大部分の観測SNP数が1以上の遺伝子に限定
min_obs <- 1L
idx_hit <- (obs_flank_all >= min_obs)
genes   <- names(obs_flank_all)[idx_hit]
obs_vec <- obs_flank_all[genes]
if (length(genes) == 0L) {
  stop("There are no gene regions where the input SNP overlaps")
}


##gr_exp もヒット遺伝子にサブセット
sel <- as.character(mcols(gr_exp)[[gene_col]]) %in% genes
gr_exp_hit <- gr_exp[sel]
gr_exp_hit <- keepSeqlevels(gr_exp_hit, keep_chr, pruning.mode="coarse")
seqlengths(gr_exp_hit) <- seqlengths(hg)[seqlevels(gr_exp_hit)]
gr_exp_hit <- IRanges::trim(gr_exp_hit)

##permutation
## ここから並列Permutation（巨大行列は作らない）
RNGkind("L'Ecuyer-CMRG"); set.seed(42)  

## 集計用ベクトル（遺伝子ごと）
sum_x  <- numeric(length(genes))
sum_x2 <- numeric(length(genes))
ge_cnt <- integer(length(genes))
names(sum_x) <- names(sum_x2) <- names(ge_cnt) <- genes

## 1回分のPermutation → 遺伝子ごとのカウント（名前付き整数ベクトル）を返す関数
perm_once <- function(i) {
  rp <- circularRandomizeRegions(
    gr_snps,
    genome = genome_ref,
    per.chromosome = TRUE,
    mask = if (use_mask) mask else NULL
  )
  
  ## FULL側
  ovpF <- findOverlaps(rp, gr_exp, ignore.strand=TRUE)
  tF <- table(as.character(mcols(gr_exp)[subjectHits(ovpF), gene_col]))
  ## ORF側
  ovpO <- findOverlaps(rp, gr_orf, ignore.strand=TRUE)
  tO <- table(as.character(mcols(gr_orf)[subjectHits(ovpO), gene_col]))
  
  ## フランキングカウント用ベクトル（必ず genes 長・names）を用意
  vFlank <- integer(length(genes))
  names(vFlank) <- genes
  
  ## FULLカウントを追加（genes に含まれるものだけ）
  if (length(tF)) {
    gF <- intersect(names(tF), genes)
    if (length(gF)) {
      vFlank[gF] <- vFlank[gF] + as.integer(tF[gF])
    }
  }
  
  ## ORFカウントを減算（genes に含まれるものだけ）
  if (length(tO)) {
    gO <- intersect(names(tO), genes)
    if (length(gO)) {
      vFlank[gO] <- vFlank[gO] - as.integer(tO[gO])
    }
  }
  
  ## 念のためマイナスは0に丸める
  vFlank[vFlank < 0L] <- 0L
  
  vFlank
}

## 並列でK回走らせ、逐次で合算（巨大オブジェクトを作らない）
## chunkごとにまとめて返し、親プロセスで加算する方式（コードは簡潔なまま）
CHUNK  <- 5000
nchunk <- ceiling(K / CHUNK)

for (b in seq_len(nchunk)) {
  this_B <- if (b < nchunk) CHUNK else (K - CHUNK * (nchunk - 1))
  idxs   <- (1:this_B) + (b - 1) * CHUNK
  
  res_list <- mclapply(idxs, perm_once, mc.cores = CORES)
  
  for (vFlank in res_list) {
    ## ★ ΣX, ΣX² を更新 → sd, z を計算可能に
    sum_x  <- sum_x  + vFlank
    sum_x2 <- sum_x2 + vFlank * vFlank
    ge_cnt <- ge_cnt + as.integer(vFlank >= obs_vec)
  }
  cat(sprintf("chunk %d/%d done (%s perms)\n",
              b, nchunk, format(this_B, big.mark = ",")))
}

##統計量
mean_hat <- sum_x / K
var_hat  <- (sum_x2 - (sum_x^2) / K) / pmax(1, K - 1)
var_hat[var_hat < 0] <- 0
sd_hat   <- sqrt(var_hat)

p_emp <- (ge_cnt + 1) / (K + 1)
z     <- ifelse(sd_hat > 0, (obs_vec - mean_hat) / sd_hat, NA_real_)

#多重検定補正
p_vec <- p_emp
ok    <- is.finite(p_vec) & p_vec >= 0 & p_vec <= 1
p_in  <- p_vec[ok]

q_bh   <- rep(NA_real_, length(p_vec))
p_holm <- rep(NA_real_, length(p_vec))

q_bh[ok]   <- p.adjust(p_in, method = "BH")
p_holm[ok] <- p.adjust(p_in, method = "holm")

names(q_bh)   <- names(p_vec)
names(p_holm) <- names(p_vec)

##結果データフレーム
res_gene <- data.frame(
  gene   = genes,
  obs    = as.integer(obs_vec[genes]),
  mean   = mean_hat[genes],
  sd     = sd_hat[genes],
  z      = z[genes],
  p_emp  = p_emp[genes],
  FDR_BH = q_bh[genes],
  HOLM   = p_holm[genes],
  row.names = NULL
)


## ログ（family size と最小p/qの確認）
cat(sprintf("# family size (hit genes) m = %d\n", length(genes)))

#全遺伝子の結果出力
res_gene <- res_gene[order(res_gene$p_emp, -res_gene$z, res_gene$gene), ]
out_csv <- sprintf(
  "all_gene_exp_outGB_u%sk_d%sk_result_rep_%s_%s_LD_%s_%s_%s_.csv",
  up_frg, down_frg, K, ld_pop, ld_r2, trait_name, trait_id
)
write.csv(res_gene, out_csv, row.names = FALSE)
cat("Saved the results (for all genes analyzed)", out_csv, "\n")


#可能性のある遺伝子の出力
##しきい値
thr_obs <- 1
thr_p   <- 0.01
thr_q <- 0.1
##条件でフィルタ：obs≥1 & p_emp≤0.01 & FDR≤0.05
hits <- subset(
  res_gene,
  p_emp <= thr_p & FDR_BH <= thr_q
)
##結果の並び替え
hits <- hits[order(hits$p_emp, -hits$z, hits$gene), ]

map_csv  <- "gene_id_name_map.csv"
gene_map <- read.csv(map_csv, stringsAsFactors = FALSE)

gene_map_collapsed <- aggregate(
  gene_name ~ gene_id,
  data = gene_map,
  FUN = function(x) paste(unique(na.omit(x)), collapse = "/")
)

hits$.__ord <- seq_len(nrow(hits))
hits_annot <- merge(
  hits,
  gene_map_collapsed,
  by.x = "gene",
  by.y = "gene_id",
  all.x = TRUE,
  sort  = FALSE
)
hits_annot <- hits_annot[order(hits_annot$.__ord), ]
hits_annot$.__ord <- NULL

hits_annot <- hits_annot[, c("gene", "gene_name",
                             setdiff(names(hits_annot), c("gene", "gene_name")))]

outfile <- sprintf(
  "significant_gene_exp_outGB_u%sk_d%sk_result_rep_%s_%s_LD_%s_%s_%s_.csv",
  up_frg, down_frg, K, ld_pop, ld_r2, trait_name, trait_id
)
write.csv(hits_annot, outfile, row.names = FALSE)
cat(sprintf("Saved the results (all genes that became significant): %s（%d genes）\n",
            outfile, nrow(hits_annot)))

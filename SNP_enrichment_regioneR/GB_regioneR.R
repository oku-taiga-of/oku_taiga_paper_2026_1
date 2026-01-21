

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
bed_tag        <- args[6]
CORES <- as.integer(args[7])

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

#遺伝子領域のbed読み込み
orf_ded <- sprintf("gene_region_ORF.bed")

gr_orf    <- rtracklayer::import(orf_ded)

#染色体処理
##染色体表記を統一
seqlevelsStyle(gr_snps)   <- "UCSC"
seqlevelsStyle(gr_orf)    <- "UCSC"
##SNPが乗っている常染色体に限定
tbl <- table(seqnames(gr_snps))
keep_chr <- names(tbl)[tbl > 0]
auto <- paste0("chr", 1:22)
keep_chr <- intersect(keep_chr, auto) 
##常染色体に限定
gr_snps  <- keepSeqlevels(gr_snps,  keep_chr, pruning.mode="coarse")
gr_orf   <- keepSeqlevels(gr_orf,   keep_chr, pruning.mode="coarse")

##染色体長を取得
seqlengths(gr_snps) <- seqlengths(hg)[seqlevels(gr_snps)]
seqlengths(gr_orf)  <- seqlengths(hg)[seqlevels(gr_orf)]

##使用可能な領域のマスクを作成
mask <- getMask(hgmask) 
mask <- keepSeqlevels(mask, keep_chr, pruning.mode="coarse") 
seqlengths(mask) <- seqlengths(hgmask)[seqlevels(mask)]  

#検定
##設定
K <- 50000
genome_ref <- hgmask
use_mask <- TRUE

##遺伝子IDの列名を取得
gene_col <- "name"

##観測値の取得
ov  <- findOverlaps(gr_snps, gr_orf, ignore.strand=TRUE)
genes_all <- unique(as.character(mcols(gr_orf)[[gene_col]]))
obs_counts <- table(as.character(mcols(gr_orf)[subjectHits(ov), gene_col]))
obs_vec_all <- setNames(integer(length(genes_all)), genes_all)
obs_vec_all[names(obs_counts)] <- as.integer(obs_counts)

##少なくともSNPが1つ載っている遺伝子に限定
idx_hit <- (obs_vec_all >= 1L)
min_obs <- 1L
idx_hit <- idx_hit & (obs_vec_all >= min_obs)

genes  <- names(obs_vec_all)[idx_hit]
obs_vec <- obs_vec_all[genes]
if (length(genes) == 0L) {
  stop("There are no gene regions where the input SNP overlaps")
}

##gr_orf もヒット遺伝子にサブセット
sel <- as.character(mcols(gr_orf)[[gene_col]]) %in% genes
gr_orf_hit <- gr_orf[sel]
gr_orf_hit <- keepSeqlevels(gr_orf_hit, keep_chr, pruning.mode="coarse")
seqlengths(gr_orf_hit) <- seqlengths(hg)[seqlevels(gr_orf_hit)]
gr_orf_hit <- IRanges::trim(gr_orf_hit)

##permutation
## ここから並列Permutation（巨大行列は作らない）
RNGkind("L'Ecuyer-CMRG"); set.seed(42)  


## 集計用ベクトル（遺伝子ごと）
sum_x   <- numeric(length(genes))         # 合計（平均用）
sum_x2  <- numeric(length(genes))         # 二乗和（分散用）
ge_cnt  <- integer(length(genes))         # ≥ obs の回数（p_emp用）
names(sum_x) <- names(sum_x2) <- names(ge_cnt) <- genes

## 1回分のPermutation → 遺伝子ごとのカウント（名前付き整数ベクトル）を返す関数
perm_once <- function(i) {
  rp <- circularRandomizeRegions(
    gr_snps,
    genome = genome_ref,
    per.chromosome = TRUE,
    mask = if (use_mask) mask else NULL
  )
  ovp <- findOverlaps(rp, gr_orf_hit, ignore.strand=TRUE)
  tab <- table(as.character(mcols(gr_orf_hit)[subjectHits(ovp), gene_col]))
  stats::setNames(as.integer(tab), names(tab))
}

## 並列でK回走らせ、逐次で合算（巨大オブジェクトを作らない）
## chunkごとにまとめて返し、親プロセスで加算する方式（コードは簡潔なまま）
CHUNK <- 5000                                # 好きな粒度。5k×10回=5万
nchunk <- ceiling(K / CHUNK)

for (b in seq_len(nchunk)) {
  this_B <- if (b < nchunk) CHUNK else (K - CHUNK*(nchunk-1))
  idxs <- (1:this_B) + (b-1)*CHUNK
  
  res_list <- mclapply(idxs, perm_once, mc.cores = CORES)
  
  ## 親側で合算
  for (cnt in res_list) {
    # 全遺伝子ゼロ埋め → 出現分だけ上書き
    v_full <- integer(length(genes)); names(v_full) <- genes
    if (length(cnt) > 0) {
      v_full[names(cnt)] <- as.integer(cnt)
    }
    # 平均・分散用の合算
    sum_x  <- sum_x  + v_full
    sum_x2 <- sum_x2 + v_full * v_full
    # ≥ obs の回数判定
    ge_cnt <- ge_cnt + as.integer(v_full >= obs_vec) 
  }
  cat(sprintf("chunk %d/%d done (%s perms)\n",
              b, nchunk, format(this_B, big.mark=",")))
}


##統計量
mean_hat <- sum_x / K
var_hat <- (sum_x2 - (sum_x^2) / K) / pmax(1, K - 1)
var_hat[var_hat < 0] <- 0  # 数値誤差で負になったら0に丸める
sd_hat <- sqrt(var_hat)
p_emp <- (ge_cnt + 1) / (K + 1)
z_raw <- ifelse(sd_hat > 0, (obs_vec - mean_hat) / sd_hat, NA_real_)
##多重検定補正
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
  gene  = genes,
  obs   = as.integer(obs_vec[genes]),
  mean  = mean_hat[genes],
  sd    = sd_hat[genes],
  z     = z_raw[genes],
  p_emp = p_emp[genes],
  FDR_BH   = q_bh[genes],
  HOLM  = p_holm[genes],
  row.names = NULL
)

#全遺伝子の結果出力
res_gene <- res_gene[order(res_gene$p_emp, -res_gene$z, res_gene$gene), ]
out_csv <- sprintf("all_gene_GB_result_rep_%s_%s_LD_%s_%s_%s_.csv", K, ld_pop, ld_r2, trait_name, trait_id )
write.csv(res_gene, out_csv, row.names = FALSE)
cat("Saved the results (for all genes analyzed)", out_csv, "\n")

#可能性のある遺伝子の出力
##しきい値
thr_obs <- 1
thr_p   <- 0.01
thr_fdr <- 0.1
##条件でフィルタ：obs≥1 & p_emp≤0.01 & FDR≤0.1
hits <- subset(res_gene, obs >= thr_obs & p_emp <= thr_p & FDR_BH <= thr_fdr)
##見やすい並び（p小さい→z大きい→gene名）
hits <- hits[order(hits$p_emp, -hits$z, hits$gene), ]
##ファイル名に条件やK/CORESを埋め込む（記録しやすい）
outfile <- sprintf("significant_gene_GB_result_rep_%s_%s_LD_%s_%s_%s_.csv", K, ld_pop, ld_r2, trait_name, trait_id )
##遺伝子名の紐づけCSVを読み込む
map_csv <- "gene_id_name_map.csv"
gene_map <- read.csv(map_csv, stringsAsFactors = FALSE)

##gene_id ごとに gene_name を "/" で結合（重複対応）
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

##geneの直後にgene_nameを配置
hits_annot <- hits_annot[, c("gene","gene_name", setdiff(names(hits_annot), c("gene","gene_name")))]

##csvとして出力
write.csv(hits_annot, outfile, row.names = FALSE)
cat(sprintf("Saved the results (all genes that became significant): %s（%d genes）\n", outfile, nrow(hits_annot)))

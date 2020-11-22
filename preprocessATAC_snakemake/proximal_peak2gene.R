### Proximal peaks 2 genes ###
suppressPackageStartupMessages({
  library(argparse)
  library(Signac)
  library(GenomicRanges)
  library(stringr)
  library(Matrix)
  # library(cellatacUtils)
  library(ensembldb)
  library(EnsDb.Hsapiens.v86)
  library(EnsDb.Mmusculus.v79)
})

parser <- ArgumentParser()
parser$add_argument("atac_feats", type="character",
                    help = "Path to csv file contaning ATAC peak ids")
parser$add_argument("rna_feats", type="character",
                    help = "Path to csv file contaning RNA gene ids")
# parser$add_argument("outdir",
#                     type="character",
#                     help = "Path to output directory")
parser$add_argument("--genome", default = "hg38",
                    type="character",
                    help = "Reference genome ID")
parser$add_argument("--prox_window", default = 50000,
                    type="integer",
                    help = "Window around genes to detect proximal peaks")
args <- parser$parse_args()

#' Find peaks close to features of interest
#'
#' @param peaks_gr GenomicRanges object containing peaks
#' @param features_gr GenomicRanges object containing features (e.g. genes)
#' @param d distance to include peak, in bps (default 50000)
#' @param feat_anno column in `features_gr@elementMetadata` containing annotation to name features (if NULL converts Granges to string)
#'
#' @return Sparse adjacency matrix indicating hits
peak2feature <- function(peaks_gr, features_gr, d=50000, feat_anno=NULL){
  seqlevelsStyle(features_gr) <- seqlevelsStyle(peaks_gr)
  
  ## Find peaks overlapping the search range around the features
  ext_gr <- Signac::Extend(features_gr, upstream = d, downstream = d)
  ovs <- findOverlaps(peaks_gr, ext_gr)
  
  ## Define identifiers for peaks and features
  all_peaks <- GRangesToString(peaks_gr, sep = c(":", '-'))
  if (is.null(feat_anno)) {
    all_feats <- GRangesToString(features_gr, sep = c(":", '-'))
  } else {
    all_feats <- features_gr@elementMetadata[[feat_anno]]
  }
  
  ## Build adjacency matrix for hits
  adj_mat <- Matrix(data=0, nrow = length(all_peaks), ncol=length(all_feats))
  for (i in unique(subjectHits(ovs))) {
    # if (length(adj_mat[queryHits(ovs[subjectHits(ovs)==i]),i]) > 0) {
    adj_mat[queryHits(ovs[subjectHits(ovs)==i]),i] <- 1
    # }
  }
  colnames(adj_mat) <- all_feats
  rownames(adj_mat) <- all_peaks
  
  adj_mat
  
}

atac.feats.file <- args$atac_feats
rna.feats.file <- args$rna_feats
# outdir <- args$outdir
prox_win <- args$prox_window
genome <- args$genome


## Define reference gene annotation

if (genome=="hg38") {
  EnsDb.genome = EnsDb.Hsapiens.v86
} else if (genome=="mm10"){
  EnsDb.genome = EnsDb.Mmusculus.v79 
} else {
  stop("Unrecognized genome ID. Supported genomes: hg38, mm10")
}

genes_gr <- genes(EnsDb.genome)

## Load peaks
peaks <- read.csv(atac.feats.file)[["peak_id"]]
peaks_gr <- StringToGRanges(peaks, sep=c(":", "-"))

## Load genes
feats <- read.csv(rna.feats.file)[["X0"]]
feats_gr <- genes_gr[genes_gr$gene_name %in% feats]

## Compute peak2gene mat
adj_mat <- peak2feature(peaks_gr, feats_gr, feat_anno = "gene_name", d=prox_win)

## Write adjacency matrix
prefix <- str_remove(atac.feats.file, "ATAC.features.csv")

writeMM(adj_mat, paste0(prefix, "peak2gene.mmtx"))
write(x = rownames(adj_mat), file = paste0(prefix, "peak2gene.peaks.tsv"))
write(x = colnames(adj_mat), file = paste0(prefix, "peak2gene.genes.tsv"))

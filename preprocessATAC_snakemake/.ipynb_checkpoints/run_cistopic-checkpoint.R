### Run cisTopic ### 
suppressPackageStartupMessages({
    library(cisTopic)
    library(Matrix)
    library(readr)
    library(tibble)
    library(argparse)
    library(data.table)
    library(ensembldb)
    library(EnsDb.Hsapiens.v86)
    library(EnsDb.Mmusculus.v79)
})

parser <- ArgumentParser()
parser$add_argument("count_mat", type="character",
                    help = "Path to tsv file for input count matrix")
# parser$add_argument("outdir", type="character",
#                     help = "Path to output directory")
parser$add_argument("--n_cores", default = 20,
                    type="integer",
                    help = "number of cores")
parser$add_argument("--seed", default = 2020,
                    type="integer",
                    help = "random seed")
parser$add_argument("--genome", default = 'hg38',
                    type="character",
                    help = "Reference genome (for gene annotations)")
args <- parser$parse_args()

countpath <- args$count_mat
# outdir <- args$outdir
n_cores <- args$n_cores
seed <- args$seed

# countpath <- "/nfs/team205/ed6/data/paired_RNA_ATAC_datasets/Chen_2019/P0_brain/GSE126074_P0_BrainCortex_SNAREseq_ATAC_cisTopic.tsv"
print("Starting to load ...")
count.matrix <- fread(countpath, sep = "\t", header=TRUE)
count.matrix = as.data.frame(count.matrix)
count.matrix = column_to_rownames(count.matrix, "peak_id")

## Initialize object ##
cisTopicObject <- createcisTopicObject(count.matrix, project.name='ATAC_cisTopic')

## Run models (updated to cisTopic v3) ##
topics_vec <- c(30:50, 60, 70, 80, 90, 100)
print("Starting to model...")
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=topics_vec, seed=seed, nCores=n_cores, iterations = 500)
print("Done w modelling!")
saveRDS(cisTopicObject, paste0(tools::file_path_sans_ext(countpath), ".trainedCistopic.RDS"))

## Select model ##
cistopic_model <- selectModel(cisTopicObject, type="derivative")

## Save topic matrix ##
cellassign <- modelMatSelection(cistopic_model, 'cell', 'Probability')
cistopic_df <- as.data.frame(t(cellassign)) 
cistopic_df <- rownames_to_column(cistopic_df, "cell")
write_csv(cistopic_df, paste0(tools::file_path_sans_ext(countpath), ".topics.csv"))

## Save predictive distribution ##
p <- predictiveDistribution(cistopic_model)

writeMM(p, paste0(tools::file_path_sans_ext(countpath), ".predDist.mmtx"))
write(x = rownames(p), file = paste0(tools::file_path_sans_ext(countpath), ".predDist.peaks.tsv"))
write(x = colnames(p), file = paste0(tools::file_path_sans_ext(countpath), ".predDist.cells.tsv"))

## Calculate gene scores ##
# Following the method from Bravo-Gonzales et al. 2020
peak2feature <- function(peaks_gr, features_gr, d=50000, feat_anno=NULL){
  seqlevelsStyle(features_gr) <- seqlevelsStyle(peaks_gr)
  
  # Find peaks overlapping the search range around the features
  ext_gr <- Signac::Extend(features_gr, upstream = d, downstream = d)
  ovs <- findOverlaps(peaks_gr, ext_gr)
  
  # Define identifiers for peaks and features
  all_peaks <- GRangesToString(peaks_gr, sep = c(":", '-'))
  if (is.null(feat_anno)) {
    all_feats <- GRangesToString(features_gr, sep = c(":", '-'))
  } else {
    all_feats <- features_gr@elementMetadata[[feat_anno]]
  }
  
  # Build adjacency matrix for hits
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

# Define reference gene annotation
if (genome=="hg38") {
  EnsDb.genome = EnsDb.Hsapiens.v86
} else if (genome=="mm10"){
  EnsDb.genome = EnsDb.Mmusculus.v79 
} else {
  stop("Unrecognized genome ID. Supported genomes: hg38, mm10")
}

genes_gr <- genes(EnsDb.genome)
peaks_gr <- StringToGRanges(rownames(p), sep=c(":", "-"))

adj_mat <- peak2feature(peaks_gr, genes_gr, feat_anno = "gene_name", d=5000)
keep.genes <- which(colSums(adj_mat) > 2) # Filter genes with at least 3 hits
gene_mat <- sapply(keep.genes, function(g) apply(p[which(adj_mat[,g]==1),],2, median)) # This is slow... might optimize in the future
gene_mat <- t(gene_mat)
                   
writeMM(gene_mat, paste0(tools::file_path_sans_ext(countpath), ".geneScores.mmtx"))
write(x = rownames(gene_mat), file = paste0(tools::file_path_sans_ext(countpath), ".geneScores.genes.tsv"))
write(x = colnames(gene_mat), file = paste0(tools::file_path_sans_ext(countpath), ".geneScores.cells.tsv"))

                
                   

                   
                   


                   
                   
                   
                   
### Run cisTopic ### 
library(cisTopic)
library(Matrix)
library(readr)
library(tibble)
library(argparse)

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
args <- parser$parse_args()

countpath <- args$count_mat
n_cores <- args$n_cores
seed <- args$seed

# countpath <- "/nfs/team205/ed6/data/paired_RNA_ATAC_datasets/Chen_2019/P0_brain/GSE126074_P0_BrainCortex_SNAREseq_ATAC_cisTopic.tsv"

count.matrix <- read.table(countpath, sep = "\t", header=TRUE, row.names = "peak_id")

# Transpose matrix (to calculate probability over the population)
peaks <- rownames(count.matrix)
cells <- colnames(count.matrix)

count.matrix.t <- t(count.matrix)

## Dummy feature names to make cisTopic work
rownames(count.matrix.t) <- paste0("chr1:1000000-", 1000000 + 1:nrow(count.matrix.t))

# Initialize object
cisTopicObject <- createcisTopicObject(count.matrix.t, project.name='ATAC_cisTopic')
# Run models (updated to cisTopic v3)
topics_vec <- c(30:50, 60, 70, 80, 90, 100)
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=topics_vec, seed=seed, nCores=n_cores, iterations = 500)
saveRDS(cisTopicObject, paste0(tools::file_path_sans_ext(countpath), ".trainedTransposedCistopic.RDS"))

# # Select model
# cistopic_model <- selectModel(cisTopicObject, type="derivative")
# # cellassign <- modelMatSelection(cistopic_model, 'cell', 'Probability')
# # cistopic_df <- as.data.frame(t(cellassign)) 
# # cistopic_df <- rownames_to_column(cistopic_df, "cell")

# # Calculate accessibility probability based on model
# p <- predictiveDistribution(cistopic_model)
# prob.matrix <- t(p)
# # rownames(prob.matrix) <- peaks
# # colnames(prob.matrix) <- cells

# writeMM(prob.matrix, paste0(tools::file_path_sans_ext(countpath), ".probabilityT.mmtx"))
# write(x = peaks, file = paste0(tools::file_path_sans_ext(countpath), ".probabilityT.peaks.tsv"))
# write(x = cells, file = paste0(tools::file_path_sans_ext(countpath), ".probabilityT.cells.tsv"))

# # write_csv(cistopic_df, paste0(tools::file_path_sans_ext(countpath), ".topics.csv"))

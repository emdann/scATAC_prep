## Pre-processing of scATAC-seq data

This repo contains a collection of scripts and notebooks that I use for initial pre-processing and dimensionality reduction of 10X scATAC-seq data. These handle steps that come **after** peak calling with the [cellatac pipeline](https://github.com/cellgeni/cellatac) (CellGenIT can run this for you). 

### Input 
After running `cellatac` you will have a results folder. From that we will need the contents of `results/peak_matrix`

### Setting up

Make a new folder to store all your results

```
mkdir my_scATAC_dir
outdir=./my_scATAC_dir
cellatac_outs=results/peak_matrix
```

Go the the folder w the scripts
```
cd scATAC_preprocess/preprocessATAC_snakemake
```

Install required R packages (other packages are installed directly in the notebooks)
```
install.packages(c("data.table","Signac", "tidyverse", "argparse"))

if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}
BiocManager::install(c("ensembldb", "EnsDb.Mmusculus.v79", "EnsDb.Hsapiens.v86","GenomicRanges"))

```

Install required python packages:
```
pip install scanpy anndata numpy pandas anndata2ri rpy2 scipy
```


#### 1. Annotate peaks 

This script uses functionality in [`GenomicRanges`](https://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html) to compute some basic stats on the peaks identified by `cellatac`. We will use these for filtering and in later stages of the analysis. 

In terminal:
```
Rscript annotate_peaks.R $cellatac_outs/peaks.txt $outdir --genome hg38
```

#### 2. Make anndata and filter peaks

See notebook: `cellatac2anndata.ipynb`

Here we make an anndata object, we filter out a lot of peaks and start saving different layers (loading the big matrix takes some time, be patient).

The notebook saves an `.h5ad` abject storing the raw counts and peak annotations for the dataset.

#### 3. Run cisTopic 

See notebook: `add_cistopic.ipynb`

_cisTopic_ uses Latent Dirichlet Allocation (LDA) to perform dimensionality reduction on a binary (or sort of binary) matrix of peak x cell accessibility. LDA is a robust Bayesian method used in text mining to group documents addressing similar topics and related words into topics. _cisTopic_ treats cells as documents and peaks as words to group cells where the same "regulatory topics" are active. Read all about it [here](https://github.com/aertslab/cisTopic).

Here we use the LDA model implemented in cisTopic for 3 purposes

1. Dimensionality reduction: the topic x cell matrix obtained with LDA can be used to construct a KNN graphs, UMAPs, clustering etc... (basically an alternative to PCA in scRNA-seq pipelines)
2. Data de-noising: scATAC-seq data from 10X Genomics is very sparse, de-noising based on the LDA model helps in exploratory data analysis and when examining accessibility patterns of single-peaks
3. Calculating accessibility scores per gene: the denoised peak x cell matrix can be used to aggregate signal over genes, useful to examine markers 

References for uses of cisTopic: 
- [Bravo Gonzales-Blas et al. (2017) cisTopic: cis-regulatory topic modeling on single-cell ATAC-seq data, Nat Methods](https://www.nature.com/articles/s41592-019-0367-1)
- [Bravo Gonzales-Blas et al. (2020) Identification of genomic enhancers through spatial integration of single‐cell transcriptomics and epigenomics. Mol. Syst. Biol.](https://www.embopress.org/doi/full/10.15252/msb.20209438) 
 

#### 4. Match genes to proximal peaks

See (short) notebook: `peak2genes.ipynb`

A range of downstream analyses require to associate genes to peaks in their proximity (i.e. overlapping the gene body or less than _n_ kbs away). The script `proximal_peak2gene.R` uses functionality in the `GenomicRanges` R package to create a sparse matrix of peak x gene assignment, where 1s indicate that a peak is in the proximity of the gene (N.B. A peak might be proximal to multiple genes). The default proximity window is 50 kbs.

## FAQs

**Why cisTopic and not x,y,z...?** 



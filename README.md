## Pre-processing of scATAC-seq data

This repo contains a collection of scripts and notebooks that I use for initial pre-processing and dimensionality reduction of 10X scATAC-seq data. These handle steps that come **after** peak calling with the [cellatac pipeline](https://github.com/cellgeni/cellatac) (CellGenIT can run this for you). 

### Input 
After running `cellatac` you will have a results folder. From that we will need the contents of `results/peak_matrix`

- A file containing the peak ids (peak_matrix/peaks.txt)
- A peak x cell count matrix

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

Install a bunch of R packages. In R:
```
install.packages("Signac")

if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}
BiocManager::install(c("ensembldb", "EnsDb.Mmusculus.v79", "EnsDb.Hsapiens.v86","GenomicRanges"))

```


#### Annotate peaks (R) 

This script calculates some basic stats on the peaks identified by `cellatac` that we will use for filtering and in later stages of the analysis. 

In terminal:
```
Rscript annotate_peaks.R $cellatac_outs/peaks.txt $outdir --genome hg38
```

#### Make anndata 

See notebook: `cellatac2anndata.ipynb`

Here we make an anndata object, we filter out a lot of peaks and start saving different layers (loading the big matrix takes some time, be patient).

The notebook saves (A) An `.h5ad` abject storing the raw counts and peak annotations for the dataset; (B) The binary counts matrix in `.mtx`, that is used to run cisTopic (lots of time-consuming I/O here, could and should be improved)

#### Run cisTopic 


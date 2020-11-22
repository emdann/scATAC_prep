## Scripts to pre-process scATAC-seq data

#### Input 
After running `cellatac` you will have a results folder. From that we will need the contents of `results/peak_matrix`

- A file containing the peak ids (peak_matrix/peaks.txt)
- A peak x cell count matrix

#### Setting up

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

#### Annotate peaks (R) 

This script calculates some basic stats on the peaks identified by `cellatac` that we will use for filtering and in later stages of the analysis

In terminal:
```
Rscript annotate_peaks.R $cellatac_outs/peaks.txt $outdir --genome hg38
```

#### Make anndata (python)

Here we make an anndata object, we filter out a lot of peaks and start saving different layers (loading the big matrix takes time, be patient)

```

```

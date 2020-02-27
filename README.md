## hicrep
R package to evaluate the reproducibility of Hi-C data
(Genome Research 2017. doi: 10.1101/gr.220640.117.)

Tao Yang  <xadmyangt@gmail.com>

## Introduction

Hi-C data analysis and interpretation are still in their early stages. In particular, there has been a lack of sound statistical metric to evaluate the quality of Hi-C data. When biological replicates are not available, investigators often rely on either visual inspection of Hi-C interaction heatmap or examining the ratio of long-range interaction read pairs over the total sequenced reads, neither of which are supported by robust statistics. When two or more biological replicates are available, it is a common practice to compute either Pearson or Spearman correlation coefficients between the two Hi-C data matrices and use them as a metric for quality control. However, these kind of over-simplified approaches are problematic and may lead to wrong conclusions, because they do not take into consideration the unique characteristics of Hi-C data, such as distance-dependence and domain structures. As a result, two un-related biological samples can have a strong Pearson correlation coefficient, while two visually similar replicates can have poor Spearman correlation coefficient. It is also not uncommon to observe higher Pearson and Spearman correlations between unrelated samples than those between real biological replicates. 

We develop a novel framework, `hicrep`, for assessing the reproducibility of Hi-C data. It first minimizes the effect of noise and biases by smoothing Hi-C matrix, and then addresses the distance-dependence effect by stratifying Hi-C data according to their genomic distance. We further adopt a stratum-adjusted correlation coefficient (SCC) as the measurement of Hi-C data reproducibility. The value of `SCC` ranges from -1 to 1, and it can be used to compare the degrees of differences in reproducibility. Our framework can also infer confidence intervals for `SCC`, and further estimate the statistical significance of the difference in reproducibility measurement for different data sets.


## Citation

Cite our paper:

HiCRep: assessing the reproducibility of Hi-C data using a 
stratum-adjusted correlation coefficient. Tao Yang, Feipeng Zhang, Galip
Gurkan Yardimci, Fan Song, Ross C Hardison, William Stafford Noble, 
Feng Yue, Qunhua Li. Genome Research 2017. doi: 10.1101/gr.220640.117.

## Rationale of method

This is a 2-step method (Figure1). In Hi-C data it is often difficult to achieve sufficient coverage. When samples are not sufficiently sequenced, the local variation introduced by under-sampling can make it difficult to capture large domain structures. To reduce local variation, we first smooth the contact map before assessing reproducibility. Although a smoothing filter will reduce the individual spatial resolution, it can improve the contiguity of the regions with elevated interaction, consequently enhancing the domain structures. We use a 2D moving window average filter to smooth the Hi-C contact map. This choice is made for the simplicity and fast computation of mean filter, and the rectangular shape of Hi-C compartments.
 
In the second step, we stratify the Hi-C reads by the distance of contacting loci, calculate the Pearson correlations within each stratum, and summarize the stratum-specific correlation coefficients into an aggregated statistic. We name it as Stratum-adjusted Correlation Coefficient (SCC). For the methodology details, please refer to our manuscript.

Figure1. `hicrep` pipeline schematic representation
                          
![Figure1. `hicrep` pipeline schematic representation](https://github.com/MonkeyLB/hicrep/blob/master/vignettes/hicrep-pipeline.JPG)


## System requirements and installation 

Hardware requirement: Linux, MacOS, Windows systems

Software requirement: R version > 3.3.0

Dependencies: `stat`, `strawr`, `rhdf5`, `rmarkdonw`, `testhat`

Installation: HiCRep can be installed in three ways by running the following commands in R 

1. Install from Github:
```
devtools::install_github(“MonkeyLB/hicrep”)
```
2. Install from Bioconductor:
```
if (!requireNamespace("BiocManager", quietly = TRUE))  install.packages("BiocManager")
BiocManager::install("hicrep")
```
3. Install from source:
Download the source package `hicrep_1.11.0.tar.gz` from Github.

(NOTE: version number can change, please find the newest one.)

Link: https://github.com/qunhualilab/hicrep/raw/master/hicrep_1.11.0.tar.gz
```
> install.packages("/PATH/TO/SOURCE/hicrep_1.11.0.tar", repo = NULL, type = "source")
```

## Input file format comaptibility

HiCRep takes a pair of intra-chromosomal contact matrices as the input. Each matrix contains the contacts from a single chromosome. HiCRep now can take the contact matrix in multiple formats.

1. Squared matrix format: This is the default format for HiCRep. It can read the files in this format directly.
2. 3-column matrix indexed by genome coordinates: this format can be converted to the squared matrix format using the `bed2mat()` function in HiCRep.

```
# hic.bin.txt is a 3-column contact matrix indexed by bin numbers
bed <- read.delim(“hic.bin.txt”)
mat <- bed2mat(bed, resol = 40000) # resol specifies the bin size in base pair
```

3. 3-column matrix indexed by bin number: this format can also be converted to the squared matrix format using the `bed2mat()` function by set `resol = “NONE”`:

```
# hic.coordinate.txt is a 3-column contact matrix indexed by bin numbers
bed <- read.delim(“hic.coordinate.txt”)
mat <- bed2mat(bed, resol = “NONE”)
```

4. `.hic` file: the `hic2mat()` function directly read in `.hic` format file and convert it to squared matrix:

```
# To work with .hic file, the package 'strawr' is needed. Install it:
remotes::install_github("aidenlab/straw/R")

# test.hic is example .hic format file
mat <- hic2mat("test.hic", chromosome1 = “18”, chromosome2 = “18”, resol =   25000, method = "NONE") 

# Note: the resolution has to be precomputed and stored in the .hic file

# method = “NONE” will extract the counts. Users may extract normalized counts too.
```

5. `.cool` file: HiCRep provides a function that reads and converts a `.cool` file to a squared matrix for a user-specified chromosome:

```
# Read a .cool file, extract intrachromosomal interactions on chr18 and convert it into a squared matrix
mat <- cool2matrix(“hic.cool”, chr = ‘chr18’)
```

## Compute reproducibility score using HiCRep
HiCRep computes the reproducibility score for a pair of chromosomal matrices using the function get.scc(). The code below demonstrates the analysis for a example data: 

### 1. Compute one chromosome
```
# Read the HiC contact maps, which are in the squared matrix format
mat1 = read.table(“Example.HiC.Rep1.chr18.matrix”)
mat2 = read.table(“Example.HiC.Rep2.chr18.matrix”)
# Compute SCC for one chromosome
scc.out = get.scc(mat1, mat2, resol = 40000, h = 5, lbr = 0, ubr = 5000000)
```

### 2. Compute all chromosomes

```
all.scc <- list()
for (i in paste(“chr”, c(as.character(1:19), “X”, “Y”)){
    mat1.chr = read.table(paste0(“Example.HiC.Rep1.”, i, “.matrix”))
    mat2.chr = read.table(paste0(“Example.HiC.Rep2.”, i, “.matrix”))
    all.scc[[i]] = get.scc(mat1.chr, mat2.chr, 40000, 5, lbr = 0, ubr = 5000000)
}
```

## HiCRep parameters

HiCRep has several user-adjustable parameters: 

`resol` is the resolution of the Hi-C data, i.e. bin size. For example, resol= 40000 for data with 40kb resolution. 
lbr and ubr are the lower and upper bounds of the genomic distance between interaction loci.  Using these bounds, users can limit the reproducibility assessment to the interactions within a given range of genomic distance. In this example, the interactions whose interacting loci are no more than 5M apart are considered.

`h` is the smoothing parameter that controls the smoothing level. A larger value of `h` leads to a higher level of smoothing.  When `h=0`, no smoothing is applied. The optimal parameter choice can be obtained by running the function `htrain()` on a pair of reasonably deeply sequenced data. Users may use these values for their analyses if training is not feasible. Below is an example code of htrain() using training samples mat1 and mat2.  It will select the best h in the range of 0-10. 

```
h_value <- htrain(mat1, mat2, resol = 40000, lbr = 0, ubr = 5000000, range = 0:10)
```

Training the smoothing parameter `h` will increase the computational time. However, it is not necessary to train `h` every time. In general, for a given resolution, the value of `h` trained from a pair of deeply sequenced biological replicates can be used for other datasets with the same resolution. Here we provide the `h` values trained from two hESC replicates in Dixon et al 2015 (GEO accession: GSE52457), which were sequenced at 500 million reads, for a range of resolutions. Users may directly use these values.

```
Resolution ----- h
10kb ----- 20
25kb ----- 10
40kb ----- 5
100kb ----- 3
500kb ----- 1 or 2
1Mb ----- 0 or 1
```

For a given HiC dataset, a higher resolution matrix usually requires more smoothing, i.e. a higher `h` value, to enhance its domain structures, due to the increasing level of sparsity in the data. To compare reproducibility between samples with the same resolution, the same smoothing parameter should be used. 

## hicrep
R package to evaluate the reproducibility of Hi-C data
(Genome Research 2017. doi: 10.1101/gr.220640.117.)

Tao Yang  <xadmyangt@gmail.com>

## Introduction

Hi-C data analysis and interpretation are still in their early stages. In particular, there has been a lack of sound statistical metric to evaluate the quality of Hi-C data. When biological replicates are not available, investigators often rely on either visual inspection of Hi-C interaction heatmap or examining the ratio of long-range interaction read pairs over the total sequenced reads, neither of which are supported by robust statistics. When two or more biological replicates are available, it is a common practice to compute either Pearson or Spearman correlation coefficients between the two Hi-C data matrices and use them as a metric for quality control. However, these kind of over-simplified approaches are problematic and may lead to wrong conclusions, because they do not take into consideration of the unique characteristics of Hi-C data, such as distance-dependence and domain structures. As a result, two un-related biological samples can have a strong Pearson correlation coefficient, while two visually similar replicates can have poor Spearman correlation coefficient. It is also not uncommon to observe higher Pearson and Spearman correlations between unrelated samples than those between real biological replicates. 

we develop a novel framework, `hicrep`, for assessing the reproducibility of Hi-C data. It first minimizes the effect of noise and biases by smoothing Hi-C matrix, and then addresses the distance-dependence effect by stratifying Hi-C data according to their genomic distance. We further adopt a stratum-adjusted correlation coefficient (SCC) as the measurement of Hi-C data reproducibility. The value of SCC ranges from -1 to 1, and it can be used to compare the degrees of differences in reproducibility. Our framework can also infer confidence intervals for SCC, and further estimate the statistical significance of the difference in reproducibility measurement for different data sets. 

Please read the package vignette for a guide through of `hicrep` analysis: 
https://github.com/MonkeyLB/hicrep/blob/master/vignettes/hicrep-vigenette.Rmd


## Citation

Cite our paper:

HiCRep: assessing the reproducibility of Hi-C data using a 
stratum-adjusted correlation coefficient. Tao Yang, Feipeng Zhang, Galip
Gurkan Yardimci, Fan Song, Ross C Hardison, William Stafford Noble, 
Feng Yue, Qunhua Li. Genome Research 2017. doi: 10.1101/gr.220640.117.


## Installation
## Installation

Download the source package (hicrep_xxx.tar.gz) from Github.
Or install it from Bioconductor:
```
## try http:// if https:// URLs are not supported
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("hicrep")
```

## Rationale of method

This is a 2-step method (Figure1). In Hi-C data it is often difficult to achieve sufficient coverage. When samples are not sufficiently sequenced, the local variation introduced by under-sampling can make it difficult to capture large domain structures. To reduce local variation, we first smooth the contact map before assessing reproducibility. Although a smoothing filter will reduce the individual spatial resolution, it can improve the contiguity of the regions with elevated interaction, consequently enhancing the domain structures. We use a 2D moving window average filter to smooth the Hi-C contact map. This choice is made for the simplicity and fast computation of mean filter, and the rectangular shape of Hi-C compartments.
 
In the second step, we stratify the Hi-C reads by the distance of contacting loci, calculate the Pearson correlations within each stratum, and then summarize the stratum-specific correlation coefficients into an aggregated statistic. We name it as Stratum-adjusted Correlation Coefficient (SCC). For the methodology details, please refer to our manuscript.

Figure1. `hicrep` pipeline schematic representation
                          
![Figure1. `hicrep` pipeline schematic representation](https://github.com/MonkeyLB/hicrep/blob/master/vignettes/hicrep-pipeline.JPG)

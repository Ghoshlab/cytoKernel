cytoKernel
====

### `cytoKernel`: An R package for differential expression using kernel-based tests for high-dimensional biological data.

`cytoKernel` computes the 
feature-wise p values and their corresponding 
adjusted p values in high-dimensional single cell experiments.

### Method

>Ghosh, T., Baxter, R.M., Seal, S., Lui, V.G., Rudra, P., Vu, T., Hsieh, E.W. and Ghosh, D., 2024. cytoKernel: Robust kernel >embeddings for assessing differential expression of single cell data. bioRxiv.
> <https://doi.org/10.1101/2024.08.16.608287>

### cytoKernel Manuscript source code

>Ghosh, T. (2024). cytoKernel: cytoKernel manuscript source code.
> Zenodo. <https://doi.org/10.5281/zenodo.14003780>


### Installing cytoKernel (development version)

The R-package **cytoKernel** can be installed from GitHub using the R package
[remotes](https://github.com/r-lib/remotes):
  
  Use to install the development version of **cytoKernel** from GitHub:
```r
  if (!require("remotes")) install.packages("remotes")
remotes::install_github("Ghoshlab/cytoKernel@devel")
```

# cytoKernel Analysis Workflow

This README provides a step-by-step guide to preprocess single-cell data and perform analyses using the `cytoKernel` package alongside other complementary tools.

---
  
  ## Prerequisites
  
  The following libraries are required to run the script:
  
  ```r
library(cytoKernel)
suppressMessages({
  library(scater)
  library(muscat)
  library(sctransform)
  library(BASiCS)
  library(Linnorm)
})
library(muscData)
library(scRNAseq)
library(SingleCellExperiment)
```

### Load and Prepare Data

We begin by loading the Bacher T Cell Data using the scRNAseq package. The data is converted into a SingleCellExperiment (SCE) object with appropriate metadata. 

```r
data2 <- BacherTCellData()
rownames(data2)

coldata2 <- colData(data2)
sce <- SingleCellExperiment(
  assay = list(counts = data2@assays@data$counts),
  colData = list(
    ind = factor(coldata2$sample),
    stim = factor(coldata2$severity),
    cluster = factor(coldata2$seurat_clusters),
    cell = factor(coldata2$new_cluster_names)
  )
)

rownames(sce) <- rownames(data2)
```

### Filter and Preprocess Data


- Select cells based on conditions:
  - Include only mild-moderate and severe samples.
- Remove unassigned cells and ensure a minimum of 20 non-zero counts per gene.

```r

sel <- sce$stim == c("mild-moderate", "severe")
sce <- sce[, sel]
sce <- sce[, !is.na(sce$cell)]

sce$stim <- factor(sce$stim)
sce$cluster_id <- sce$cell
sce$sample_id <- paste(sce$stim, sce$ind)
sce$group_id <- sce$stim
colData(sce) <- colData(sce)[, -c(1:4)]
sce <- sce[rowSums(counts(sce) > 0) >= 20, ]
```

### Subset Data
- Subsample up to 200 cells per unique sample to create a subset of the Bacher COVID data for illustrative example.
```r
subsample_cells <- function(sce, n = 200) {
  sample_ids <- unique(sce$sample_id)
  subsampled_list <- list()
  
  for (sample in sample_ids) {
    sample_indices <- which(sce$sample_id == sample)
    n_to_sample <- min(length(sample_indices), n)
    sampled_indices <- sample(sample_indices, n_to_sample)
    subsampled_list[[sample]] <- sce[, sampled_indices]
  }
  
  do.call(cbind, subsampled_list)
}

subsampled_sce <- subsample_cells(sce, n = 200)
```

### Normalize Counts
- Normalize counts and compute counts per million (CPM).
```r
Bacher_COVID_subset <- subsampled_sce[1:100, ]
Bacher_COVID_subset <- computeLibraryFactors(Bacher_COVID_subset)
Bacher_COVID_subset <- logNormCounts(Bacher_COVID_subset)
assays(Bacher_COVID_subset)$cpm <- calculateCPM(Bacher_COVID_subset)
```

### Prepare Data for muscat single cell experiment format

- Prepare the SingleCellExperiment object for use with the muscat package.
```r
Bacher_COVID_subset <- prepSCE(
  Bacher_COVID_subset, 
  "cluster_id", "sample_id", "group_id", TRUE
)
```

### Create metadata
- Extract experimental metadata.
```r
samples <- metadata(Bacher_COVID_subset)$experiment_info$sample_id
patient_id <- factor(sub(".* ", "", samples))
group <- metadata(Bacher_COVID_subset)$experiment_info$group_id
design <- model.matrix(~ group)
rownames(design) <- samples

group_factor <- as.numeric(as.factor(group)) - 1
```

### Run CytoKernel Analysis
- Run without covariates:
  
 ```r
  
Bacher_mms_cytoKernel <- cytoKernel::cytoKSCE_clusters_Fpsrf(
  Bacher_COVID_subset,
  group_factor = group_factor,
  covars = NULL,
  residuals = TRUE,
  initial_permutations = 100,
  thresholds = c(0.1, 0.01, 0.001),
  perm_increments = c(500, 2000, 10000),
  random_seed = 999,
  residual_tolerance = 0.5,
  name_assays_expression = "logcounts",
  name_cluster = "cluster_id",
  name_sample = "sample_id",
  n_cores = 8,
  FDR_method = "BH"
)
```

- Run with (dummy) covariates:
  
```r
Bacher_mms_cytoKernel_design_covar <- cytoKernel::cytoKSCE_clusters_Fpsrf(
  Bacher_COVID_subset,
  group_factor = group_factor,
  covars = data.frame(design[, 1]),
  residuals = TRUE,
  initial_permutations = 100,
  thresholds = c(0.1, 0.01, 0.001),
  perm_increments = c(500, 2000, 10000),
  random_seed = 999,
  residual_tolerance = 0.5,
  name_assays_expression = "logcounts",
  name_cluster = "cluster_id",
  name_sample = "sample_id",
  n_cores = 8,
  FDR_method = "BH"
)
```

### The output column of the cytoKSCE_clusters_Fpsrf() function includes:

- gene  (feature/gene names) 
- cluster_id (cell subpopulation)   
- p_val (p value unadjusted)
- p_adj.loc (adjusted p value within clusters [locally])   
- p_adj.glb (adjusted p value [globally])  

    

### `cytoKernel`: An R/Bioconductor package for differential expression using kernel-based score test for high-dimensional biological data.

`cytoKernel` computes the 
feature-wise p values and their corresponding 
adjusted p values in high-dimensional biological experiments.

### Methods

> Liu D, Lin X, Ghosh D."Semiparametric regression of
> multi-dimensional genetic pathway data: least-squares 
> kernel machines and linear mixed models". Biometrics. 2007;63(4):1079-1088.
> <doi:10.1111/j.1541-0420.2007.00799.x1028-1039>

> Zhan X, Patterson AD, Ghosh D. "Kernel approaches for   
> differential expression analysis of mass spectrometry-based
> metabolomics data". BMC Bioinformatics. 2015;16:77. Published 2015 Mar 11. 
> <doi:10.1186/s12859-015-0506-3>

> Liu D, Ghosh D, Lin X. "Estimation and testing for the effect > of a genetic pathway on a disease outcome using logistic
> kernel machine regression via logistic mixed models". BMC 
> Bioinf. 2008; 9(1):292.
> <https://doi.org/10.1186/1471-2105-9-292>


### Installing cytoKernel

The R-package **cytoKernel** can be installed from GitHub using the R package
[devtools](https://github.com/hadley/devtools):

Use to install the latest version of **cytoKernel** from GitHub:

    if (!require("devtools")) install.packages("devtools")
    devtools::install_github("Ghoshlab/cytoKernel")
    
It can also be installed using Bioconductor: 

```s
# install BiocManager from CRAN (if not already installed)
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

# install cytoKernel package
BiocManager::install("cytoKernel")
```    

After installation, the package can be loaded into R.
```s
library(cytoKernel)
```
    

### Using cytoKernel

The main function in the **cytoKernel** package is `CytoK()`. The `CytoK()`
function needs two required objects and three optional objects: (1)
object: a data frame or a matrix or a Summarized Experiment with one
assay object with observations (e.g., cluster-marker combinations or genes) on the rows.
and samples as the columns (e.g. let’s call it `dataSE`).
(2) group_factor: a binary categorical response variable
that represents the group condition for each sample. For example if the samples represent two different groups or conditions (e.g., before stimulation and after stimulation), provide CytoK() with a phenotype representing which columns in the object are different groups. (e.g. let’s call it `groupSamples`).
(3) lowerRho (optional) a positive value that represents the lower bound of the kernel parameter. Default is 2.
(4) upperRho (optional) a positive value that represents the upper bound of the kernel parameter. Default is 12.
(5) gridRho (optional) a positive value that represents the number of grid points in the interval of upper and bound of the kernel parameter. Default is 4.
(6) alpha (optional) level of significance to control the False Discovery
rate (FDR). Default is 0.05.
(7) featureVars (optional) Vector of the columns which identify features. If a SummarizedExperiment is used for data, row variables will be used. Default is NULL.

To run the `CytoK()` function,

    CytoKOutput <- CytoK(object = dataSE,
    group_factor =groupSamples, lowerRho=2, upperRho=12, gridRho=4,
    alpha = 0.05, featureVars = NULL)

Individual slots can be extracted using accessor methods:

    CytoKFeatures(CytoKOutput) # extracts the data.frame with shrunken effect size, shrunken effect size sd, unadjusted p value and adjusted p value for each feature
    
     CytoKFeaturesOrdered(CytoKOutput) # extracts the data.frame with shrunken effect size, shrunken effect size sd, unadjusted p value and adjusted p value for each feature ordered by unadjusted p value from low to high
     
     CytoKDEfeatures(CytoKOutput) # extracts the percent of differentially expressed features
    
    CytoKData(CytoKOutput) # extracts the original data object
    
    CytoKalpha(CytoKOutput) # extracts the specified level of significance
    
    CytoKFeatureVars(CytoKOutput) # extracts the value of featureVars
    
The heatmap of the expressed matrix of features on rows ordered by the adjusted p values from low to high can be directly
plotted using the `plotCytoK()` function.

     plotCytoK(object = CytoKOutput,
     group_factor =groupSamples,topK=K,...)

For more details, see `vignettes`.

Bug reports
===========

Report bugs as issues on the [GitHub repository new
issue](https://github.com/Ghoshlab/cytoKernel/issues/new)

Contributors
============

-   [Tusharkanti Ghosh](https://github.com/tghosh30)
-   [Victor Lui]()
-   [Pratyaydipta Rudra]()
-   [Souvik Seal]()
-   [Thao Vu]()
-   [Elena Hsieh]()
-   [Debashis Ghosh](https://github.com/ghoshd)


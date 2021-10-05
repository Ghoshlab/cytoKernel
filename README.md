cytoKernel
====

### `cytoKernel`: An R/Bioconductor package for differential expression using kernel-based score test for high-dimensional biological data.

`cytoKernel` computes the 
feature-wise p values and their corresponding 
adjusted p values in high-dimensional biological experiments.

### Method

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
function needs two required object and three optional objects: (1)
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


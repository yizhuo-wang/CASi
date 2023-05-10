# CASi: Cross-timepoint Analysis of Single-cell RNA sequencing data

<table>
<tr>
<td>
  CASi is a comprehensive R package for analyzing multi-timepoint scRNA-seq data, which provides users with: (1) cross-time points cell annotation, (2) detection of potentially novel cell types emerged over time, (3) visualization of cell population evolution, and (4) identification of temporal differentially expressed genes (tDEGs).
</td>
</tr>
</table>


![Fig1_big](https://github.com/yizhuo-wang/CASi/assets/88061326/75620fc2-b1f6-40d4-81f5-0fe0a08a3017)


## Installation

Before the installation, please make sure to install JAGS-4.x.y.exe from http://www.sourceforge.net/projects/mcmc-jags/files. 

```
require(devtools)
devtools::install_github("yizhuo-wang/CASi")
```

Alternatively, you can download the zip file and install it from using:

```
install.packages("file_name_and_path", repos = NULL, type="source")
```

OR, if you have RStudio installed, you can download the zip file and open the CASi.Rproj to install the whole package.

## Usage

Before using CASi, please make sure to set up TensorFlow from https://tensorflow.rstudio.com/install/, and that the dependent package "keras" is available to use.

```
require(keras)
require(CASi)
```

## Input

- Gene expression matrix of the initial data point/t0.
- Gene expression matrix of the following data points/t1.
- (optional) Cell labels of the initial data point.


## Vignettes

A thorough demonstration can be found at the CASi.html file under the Vignettes folder.

## To-do
- Publish on Bioconductor.


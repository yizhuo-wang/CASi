# CASi: Cross-timepoint Analysis of Single-cell RNA sequencing data

<table>
<tr>
<td>
  CASi is a comprehensive R package for analyzing multi-timepoint scRNA-seq data, which provides users with: (1) cross-time points cell annotation, (2) detection of potentially novel cell types emerged over time, (3) visualization of cell population evolution, and (4) identification of temporal differentially expressed genes (tDEGs).
</td>
</tr>
</table>

![Fig1_big](https://user-images.githubusercontent.com/88061326/236935983-1fdb3534-6670-47e3-93ee-338de92799aa.png)



## Site

This site is under reconstruction.


## Installation

Before the installation, please make sure to install JAGS-4.x.y.exe from http://www.sourceforge.net/projects/mcmc-jags/files. 

```
require(devtools)
devtools::install_github("yizhuo-wang/CASi")
```

Alternatively, you can download the zip file and install it from local.

```
install.packages("file_name_and_path", repos = NULL, type="source")
```

## Usage

Before using CASi, please make sure to set up TensorFlow from https://tensorflow.rstudio.com/install/, and that the dependent package "keras" is available to use.

```
require(keras)
require(CASi)
```

## Development

- Fork the repo
- Create a new branch (`git checkout -b improve-feature`)
- Make the appropriate changes in the files
- Add changes to reflect the changes made
- Commit your changes (`git commit -am 'Improve feature'`)
- Push to the branch (`git push origin improve-feature`)
- Create a Pull Request 

### Input




## Vignettes

A thorough demonstration can be found at .

## To-do
- Publish on Bioconductor.


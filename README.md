# tissueResolver
tissueResolver is a package for converting bulk RNA-seq datasets into virtual tissues using information from similar single cell datasets by assigning weights to true single cells, maintaining their molecular integrity. Virtual tissues can be analyzed in a similar way as conventional single cell datasets.

<div style="text-align: center">
    <img src="man/figures/schematics.png" width=500 alt="tR Pipeline" />
</div>

For details, refer to: [Oxford Academic](https://doi.org/10.1093/bioinformatics/btae709)

## Installation
Use devtools to install the package:

``` R
library(devtools)
devtools::install_github("spang-lab/tissueResolver")
```

## Usage and Documentation
We provide a detailed [vignette reproducing the results of our paper](https://github.com/spang-lab/tissueResolver-docs).


## Citation
**Virtual tissue expression analysis** by\
Jakob Simeth, Paul Hüttl, Marian Schön, Zahra Nozari, Michael Huttner, Tobias Schmidt, Michael Altenbuchinger, Rainer Spang, 2024.

For the article, refer to: [Oxford Academic](https://doi.org/10.1093/bioinformatics/btae709)
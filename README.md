# Benchmarking Study of Covariate-Adjusted FDR Methods 

This repository includes all code used to perform analyses and generate figures for the manuscript, [_A practical guide to methods controlling false discoveries in computational biology_](https://www.biorxiv.org/content/early/2018/10/31/458786). 

## Notes

Before attempting to run the code in this repository, please read the following notes.

### Cloning

We **do not** recommond cloning this entire repository. The repository includes all historical commits and is incredibly large. Instead, we recommend either making a shallow clone of the repository using the following code.

```bash
 git clone --depth=1 https://github.com/pkimes/benchmark-fdr.git
```

### SummarizedBenchmark

The code in this repository was run using an older version of the `SummarizedBenchmark` package and will not run with newer versions of the package. An older version of the package should be installed from GitHub using the following R code.

```r
devtools::install_github("areyesq89/SummarizedBenchmark", ref = "fdrbenchmark")
```

### Rendered reports

Final rendered `Rmd` files for simulations and case studies included in the manuscripts are available in [this repository](https://github.com/pkimes/benchmark-fdr-html), and are rendered [here](http://www.pkimes.com/benchmark-fdr-html/) via GitHub pages.

## Citation

Korthauer K&dagger;, Kimes PK&dagger;, Duvallet C&Dagger;, Reyes A&Dagger;, Subramanian A&Dagger;, Teng M, Shukla C, Alm EJ, Hicks SC\*.
_A practical guide to methods controlling false discoveries in computational biology._ bioRxiv. doi: https://doi.org/10.1101/458786

&dagger;: co-first authors, &Dagger;: co-second authors, \*: corresponding author

## License/Copyright
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)  
All code in this project is made available under a MIT license.  
All non-code text in this project is made available under a CC BY-NC-ND 4.0 license.

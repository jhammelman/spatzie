# spatzie: Identification of enriched motif pairs from chromatin interaction data

[![license: GPL-3](https://img.shields.io/badge/license-GPL--3-blue)](https://opensource.org/licenses/GPL-3.0) [![DOI](https://img.shields.io/badge/DOI-10.1093%2Fnar%2Fgkac036-blue.svg)](https://doi.org/10.1093/nar/gkac036) [![BioC](https://img.shields.io/badge/BioC-1.8.0-brightgreen.svg)](https://doi.org/doi:10.18129/B9.bioc.spatzie) [![platforms](https://bioconductor.org/shields/availability/release/spatzie.svg)](https://bioconductor.org/packages/release/bioc/html/spatzie.html#archives)  [![Coverage Status](https://coveralls.io/repos/github/jhammelman/spatzie/badge.svg?branch=master)](https://coveralls.io/github/jhammelman/spatzie?branch=master)

https://spatzie.mit.edu

Given a database of DNA sequence motifs representing transcription factors and enhancer promoter interaction data, spatzie performs statistical analysis to identify co-enriched transcription factors.

## Installation

The *spatzie* package is part of Bioconductor since release 3.14. To install it on your system, enter:

```
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("spatzie")
```

Alternatively, the latest version can be installed directly from this repository:

```
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("jhammelman/spatzie", build_vignettes = TRUE)
```

Note: For most use cases it is not necessary to install the *spatzie* package locally, as a substantial part of its functionality is offered as an online service at https://spatzie.mit.edu.

## Usage

For interaction data aligned to the most recent human or mouse genome assemblies (`hg38`, `hg19`, `mm10`, or `mm9`), the most common spatzie use cases are covered by the function `find_ep_coenrichment`, which is prominently featured in one of the vignettes:
```
vignette("single-call", package = "spatzie")
```

The functionality displayed in the vignette above is also available online at [spatzie.mit.edu](https://spatzie.mit.edu).

If more flexibility is required, e.g., different genome assemblies, locally cached promoter annotations, non-standard ways to filter interactions, this vignette is a good starting point:
```
vignette("individual-steps", package = "spatzie")
```

## Build status

| Platform | Status |
|------|------|
| Travis CI | [![Travis build status](https://travis-ci.com/jhammelman/spatzie.svg?branch=master)](https://travis-ci.com/jhammelman/spatzie) |
| Bioconductor 3.18 (release) | [![BioC release](https://bioconductor.org/shields/build/release/bioc/spatzie.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/spatzie/) |
| Bioconductor 3.19 (devel) | [![BioC devel](https://bioconductor.org/shields/build/devel/bioc/spatzie.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/spatzie/) |


## Citation

If you use *spatzie* in your research, please cite:

**spatzie: An R package for identifying significant transcription factor motif co-enrichment from enhancer-promoter interactions**  
Jennifer Hammelman, Konstantin Krismer, and David K. Gifford  
*Nucleic Acids Research*, Volume 50, Issue 9, 20 May 2022, Page e52; DOI: https://doi.org/10.1093/nar/gkac036

## Funding

The development of this method was supported by National Institutes of Health (NIH) grants 1R01HG008754 and 1R01NS109217, and a National Science Foundation Graduate Research Fellowship (1122374).

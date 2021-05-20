# spatzie: Identification of enriched motif pairs from chromatin interaction data

[![license: GPL-3](https://img.shields.io/badge/license-GPL--3-blue)](https://opensource.org/licenses/GPL-3.0) [![Coverage Status](https://coveralls.io/repos/github/jhammelman/spatzie/badge.svg?branch=master)](https://coveralls.io/github/jhammelman/spatzie?branch=master)

Given a database of DNA sequence motifs representing transcription factors and enhancer promoter interaction data, spatzie performs statistical analysis to identify co-enriched transcription factors.

## Installation

The *spatzie* package can be installed directly from this repository:

```
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("jhammelman/spatzie")
```

<!---
## Usage

TODO refer to vignettes
-->

## Build status

| Platform | Status |
|------|------|
| Travis CI | [![Travis build status](https://travis-ci.com/jhammelman/spatzie.svg?branch=master)](https://travis-ci.com/jhammelman/spatzie) |

## Citation

If you use *spatzie* in your research, please cite:

**spatzie: An R package for identifying significant transcription factor motif co-enrichment from enhancer-promoter interactions**  
Jennifer Hammelman, Konstantin Krismer, and David K. Gifford  
Journal name, Volume *TODO*, Issue *TODO*, *TODO*; DOI: https://doi.org/*TODO*

## Funding

The development of this method was supported by National Institutes of Health (NIH) grants *TODO*.

# spatzie: Identification of enriched motif pairs from chromatin interaction data

[![license: GPL-3](https://img.shields.io/badge/license-GPL--3-blue)](https://opensource.org/licenses/GPL-3.0)

Given a database of DNA sequence motifs representing transcription factors and enhancer promoter interaction data, spatzie performs statistical analysis to identify co-enriched transcription factors.

## Installation

The *spatzie* package can be installed directly from this repository:

```
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("jhammelman/spatzie")
```

## Citation

If you use *spatzie* in your research, please cite:

**spatzie paper**  
Jennifer Hammelman, Konstantin Krismer, and David Gifford  
Journal name, Volume *TODO*, Issue *TODO*, *TODO*; DOI: https://doi.org/*TODO*

## Funding

The development of this method was supported by National Institutes of Health (NIH) grants *TODO*.

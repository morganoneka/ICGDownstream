# ICGDownstream
The Giotto package provides some amazing functionality for spatial transcriptomics data. The [Interaction Changed Genes](https://rubd.github.io/Giotto_site/reference/findICG.html) function identifies genes that are influenced by cell proximity. This script provides a few functions for more easily parsing and visualizing the output from this function.

## Functionality:
- Read in output and subset statistically significant results
- Identify the following:
  - What genes are influenced in the largest number of interactions
  - What interactions a gene of interested is affected by
  - What genes are over- and under-expressed for each interaction studied
  - How differential expression differs for a given interaction between two samples

## How to use:
- Run `findICG` on your Giotto object(s)
- Save the output as a `.csv` file.
- Import this script into yours using `source(ICGDownstream)`

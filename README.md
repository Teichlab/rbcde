# Rank-biserial correlation

RBCDE is a Python implementation of the rank-biserial correlation coefficient [(Cureton, 1956)](https://link.springer.com/article/10.1007/BF02289138), which can be used as an effect size equivalent of the Wilcoxon test [(Kerby, 2014)](https://journals.sagepub.com/doi/full/10.2466/11.IT.3.1), which in turn was deemed to perform well on single cell data problems [(Soneson, 2018)](https://www.nature.com/articles/nmeth.4612). Using effect size analyses is recommended for problems with large population sizes [(Sullivan, 2012)](https://www.jgme.org/doi/full/10.4300/JGME-D-12-00156.1). The package comes with both a [scanpy](https://scanpy.readthedocs.io/en/latest/)-compatible version and a standalone function that ingests a data matrix and an assignment vector.

## Citation

Stay tuned!

## Installation

RBCDE depends on numpy, scipy and pandas. The package is available on pip, and can be easily installed as follows:

	pip3 install rbcde

## Usage and Documentation

RBCDE can slot into a scanpy workflow and accept an object with `log(CPM/100 + 1)` data stored as a layer or `.raw`, and the desired clustering/grouping vector as an `.obs` column:

	import rbcde
	rbcde.RBC(adata)
	degs, plot_dict = rbcde.filter_markers(adata)

`rbcde.RBC()`'s `clus_key` argument controls which `.obs` column is used for the grouping, and a combination of `layer` and `use_raw` can instruct the function to retrieve expression data from `.X`, `.layers` or `.raw`. `rbcde.filter_markers()` takes the computed coefficient values and thresholds them into a list of per-cluster markers. The thresholding can be controlled via the `thresh` argument, with a range of literature [critical values](https://en.wikipedia.org/wiki/Effect_size#Pearson_r_or_correlation_coefficient) available. A helper dictionary, compatible with the formatting scanpy plotting functions accept in the `var_names` argument.

Analogous functions exist for scanpy-independent data analysis, and can ingest any data matrix with variables as rows and observations as columns. The filtering function does not produce a helper dictionary, only yielding the marker data frame.

	results = rbcde.matrix.RBC(data, clusters, vars)
	degs = rbcde.matrix.filter_markers(results)

An HTML render of the RBCDE function docstrings, detailing all the parameters, can be accessed at [ReadTheDocs](https://rbcde.readthedocs.io/en/latest/).

## Example Notebook

[rbc_demo.ipynb](https://nbviewer.jupyter.org/github/Teichlab/rbcde/blob/master/examples/rbc_demo.ipynb) computes the rank-biserial correlation coefficient for demonstration 10X PBMC data, yielding a similar standard of markers to established approaches while reporting only ~13% of the gene total. This more compact summary does not require any heuristic filtering to obtain. The full marker export yielded by the analysis can be found at `examples/markers.csv`
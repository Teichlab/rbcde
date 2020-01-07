import rbcde.matrix
import pandas as pd
try:
	from scanpy import logging as logg
except ImportError:
	pass

def RBC(adata, clus_key='leiden', layer=None, use_raw=False):
	'''
	Compute the rank-biserial correlation coefficient for each gene in each cluster. The 
	results can be subsequently turned into a marker list via the helper function 
	``rbcde.filter_markers()``. The primary output is stored as part of either `.var` or 
	`.raw.var`, depending on whether `.raw` data is used.
	
	The rank-biserial correlation coefficient 
	`(Cureton, 1956) <https://link.springer.com/article/10.1007/BF02289138>`_ 
	can be used as an effect size equivalent of the Wilcoxon test 
	`(Kerby, 2014) <https://journals.sagepub.com/doi/full/10.2466/11.IT.3.1>`_, which in 
	turn was deemed to perform well on single cell data problems 
	`(Soneson, 2018) <https://www.nature.com/articles/nmeth.4612>`_. Using effect size 
	analyses is recommended for problems with large population sizes 
	`(Sullivan, 2012) <https://www.jgme.org/doi/full/10.4300/JGME-D-12-00156.1>`_.
	
	Input
	-----
	adata : ``AnnData``
		Needs ``log(CPM/100 + 1)`` data stored somewhere in the object (as either sparse or 
		dense), and the desired clustering/grouping vector included in `.obs`.
	clus_key : ``str``, optional (default: "leiden")
		The name of the `.obs` column containing the clustering/grouping.
	layer : ``str`` or ``None``, optional (default: ``None``)
		If specified, take the expression data from the matching ``.layers`` field. Overrides 
		``use_raw`` if provided.
	use_raw : ``bool``, optional (default: ``False``)
		If no ``layer`` was specified and this is set to ``True``, take the data from the 
		``.raw`` field of the object. Store results in ``.raw.var`` to match dimensionality.
	'''
	
	start = logg.info('computing rank-biserial correlation')
	#extract appropriate data form, along with gene names
	#layer trumps use_raw
	raw_prefix = ''
	if layer is not None:
		data = adata.layers[layer]
		genes = adata.var_names
	else:
		#try to get out .raw if possible (and specified)
		if adata.raw is not None and use_raw:
			data = adata.raw.X.copy()
			genes = adata.raw.var_names
			raw_prefix = '.raw'
		else:
			data = adata.X.copy()
			genes = adata.var_names
	#extract cluster list
	clusters = adata.obs[clus_key].copy()
	#think that's everything we need. call matrix version of function
	rbc_out = rbcde.matrix.RBC(data, clusters, genes)
	#append the fact this is RBC output to each column name and stash it in var
	for col in rbc_out.columns:
		if use_raw:
			adata.raw.var['RBC_'+col] = rbc_out[col]
		else:
			adata.var['RBC_'+col] = rbc_out[col]
	logg.info('	finished', time=start,
		deep=(
			'added\n'
			"	'RBC_' columns for each of the clusters to "+raw_prefix+".var"
		),)

def filter_markers(adata, thresh=0.5, use_raw=False):
	'''
	Filter the rank-biserial correlation coefficients computed with ``rbcde.RBC()`` to a 
	list of markers for each cluster, provided as a data frame and a Scanpy plotting compatible 
	``var_names`` cluster marker dictionaty. Returns those two objects, in this order.
	
	Input
	-----
	adata : ``AnnData``
		Needs to have been processed with ``rbcde.RBC()``.
	thresh : ``float``, optional (default: 0.5)
		The threshold value used to call markers. Literature 
		`critical values <https://en.wikipedia.org/wiki/Effect_size#Pearson_r_or_correlation_coefficient>`_ 
		can be used.
	use_raw : ``bool``, optional (default: ``False``)
		Set this to ``True`` if the raw data was used for the computation so that the 
		results can be retrieved from the correct field of the object.
	'''
	
	
	#extract the RBC results embedded in .var and remove the prefix
	if use_raw:
		results = adata.raw.var.loc[:,[i.startswith('RBC_') for i in adata.raw.var.columns]]
	else:
		results = adata.var.loc[:,[i.startswith('RBC_') for i in adata.var.columns]]
	results.columns = [i.replace('RBC_','',1) for i in results.columns]
	#call the matrix version to get a marker data frame
	degs = rbcde.matrix.filter_markers(results, thresh)
	#parse up a plotting cluster marker dictionary
	plot_dict = {}
	for clus in results.columns:
		plot_dict[clus] = degs.loc[degs['cluster']==clus, :].index
		logg.hint(str(len(plot_dict[clus]))+' markers found for cluster '+clus)
	#return both the data frame and the plot-ready form
	return degs, plot_dict
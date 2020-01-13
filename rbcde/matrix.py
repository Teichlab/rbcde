import pandas as pd
import numpy as np
from scipy import stats, sparse

def RBC(data, clusters, vars):
	'''
	Compute the rank-biserial correlation coefficient for each gene in each cluster. The 
	results can be subsequently turned into a marker list via the helper function 
	``rbcde.matrix.filter_markers()``. Returns a data frame with the coefficient value for 
	each gene in each cluster.
	
	The rank-biserial correlation coefficient 
	`(Cureton, 1956) <https://link.springer.com/article/10.1007/BF02289138>`_ 
	can be used as an effect size equivalent of the Wilcoxon test 
	`(Kerby, 2014) <https://journals.sagepub.com/doi/full/10.2466/11.IT.3.1>`_. Using effect 
	size analyses is recommended for problems with large population sizes 
	`(Sullivan, 2012) <https://www.jgme.org/doi/full/10.4300/JGME-D-12-00156.1>`_.
	
	Input
	-----
	data : ``np.array`` or ``scipy.sparse``
		Per cell normalised, if using single cell count data. Variables as rows, 
		observations as columns.
	clusters : ``np.array`` or ``list``
		A vector of cluster/group assignments for each observation.
	vars: ``np.array`` or ``list``
		A vector of variable names, for output generation purposes.
	'''
	
	#set up sparse data as CSC format for quickness of access
	#and a little helper flag to see if we need to densify or not
	issparse = False
	if sparse.issparse(data):
		data = data.tocsc()
		issparse = True
	#prepare some cluster'y stuff
	clusters_categories = np.unique(clusters)
	#store the boolean mask of cluster membership, nU1 (number of cells in each cluster)
	#and nU1*nU2 (multiply that by the number of cells not in that cluster)
	masks = {}
	nU1 = {}
	nU1nU2 = {}
	for clus in clusters_categories:
		masks[clus] = np.array([i==clus for i in clusters])
		nU1[clus] = np.sum(masks[clus])
		nU1nU2[clus] = nU1[clus]*np.sum(~masks[clus])
	#make results array
	results = np.zeros((len(vars),len(clusters_categories)))
	#okay, with all this out of the way, can finally get computing
	for i in np.arange(len(vars)):
		#pre-compute the gene's ranks once
		if issparse:
			expr = data[:,i].todense()
		else:
			expr = data[:,i]
		ranks = stats.rankdata(expr)
		for j, clus in enumerate(clusters_categories):
			#compute the actual RBC coefficient and store it
			#making active use of the fact U1+U2 = nU1*nU2
			U1 = np.sum(ranks[masks[clus]]) - nU1[clus]*(nU1[clus]+1)/2
			U2 = nU1nU2[clus] - U1
			results[i,j] = (U1-U2)/nU1nU2[clus]
	return(pd.DataFrame(results, index=vars, columns=clusters_categories))

def filter_markers(results, thresh=0.5):
	'''
	Filter the rank-biserial correlation coefficients computed with ``rbcde.matrix.RBC()`` to a 
	list of markers for each cluster. Returns a data frame of the computed markers.
	
	Input
	-----
	results : ``pd.DataFrame``
		The output of ``rbcde.matrix.RBC()``.
	thresh : ``float``, optional (default: 0.5)
		The threshold value used to call markers. Literature 
		`critical values <https://en.wikipedia.org/wiki/Effect_size#Pearson_r_or_correlation_coefficient>`_ 
		can be used.
	'''
	
	degs = pd.DataFrame(columns=['cluster','RBC'])
	for clus in results.columns:
		sub = pd.DataFrame(clus, index=results.index, columns=['cluster','RBC'])
		sub['RBC'] = results[clus]
		degs = degs.append(pd.DataFrame(sub[sub['RBC']>thresh]).sort_values(by='RBC', ascending=False))
	return(degs)
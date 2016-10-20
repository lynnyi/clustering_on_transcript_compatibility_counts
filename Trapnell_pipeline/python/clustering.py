# clustering module on TCCs


import os
import getopt
import pickle
import itertools
import multiprocessing as mp
import numpy as np
import random

nclust = 3
path_to_reads = '/home/vasilis/clustering_EQ/Trapnell/read_data/'
ncells = 271
filenames = np.loadtxt('/home/lynnyi/clustering/clustering_on_transcript_compatibility_counts/Trapnell_pipeline/Trapnells_data/Trapnell_filenames.txt',dtype=str)
base_dir = '/home/lynnyi/clustering/clustering_on_transcript_compatibility_counts/Trapnell_pipeline/'


def perform_clustering():
	# Clustering
	with open(base_dir + 'Trapnell_TCC_distribution.dat', 'rb') as infile:
	    X = pickle.load(infile)
	with open(base_dir + '/Trapnell_TCC_pairwise_distance.dat','rb') as infile:
	    D = pickle.load(infile)
	Trap_labels=np.loadtxt(base_dir + '/Trapnells_data/Trapnell_labels.txt',dtype=str)
	assert np.all(np.isclose(D,D.T))
	assert np.all(np.isclose(np.diag(D),np.zeros(np.diag(D).shape)))
	# Clustering is done using scikit-learn's implementation of affinity propagation
	# D is a symmetric N-by-N distance matrix where N is the number of cells
	from sklearn import cluster
	def AffinityProp(D,pref,damp):
	    aff= cluster.AffinityPropagation(affinity='precomputed',
					     preference=pref,damping=damp, verbose=True)
	    labels=aff.fit_predict(D)
	    return labels
	# Jensen-shannon metric used to compute distances. 
	# This code is used in get_pairwise_distances.py and is repeated here for convenience. 
	from scipy.stats import entropy
	def jensen_shannon(p, q):
	    m=0.5*p+0.5*q
	    p = np.transpose(p[p > 0])
	    q = np.transpose(q[q > 0])
	    m = np.transpose(m[m > 0])
	    return np.sqrt(entropy(m)-0.5*entropy(q)-0.5*entropy(p))
	from sklearn.metrics.pairwise import pairwise_distances
	pref = -1.3*np.ones(ncells)
	labels3=AffinityProp(-D,pref,0.95)
	return labels3

# get names of cells within each cluster
def get_names_per_cluster(labels):
	flnames_cluster = []
	for i in range(nclust):
		flnames_cluster.append([])
	for i in range(ncells):
		flnames_cluster[labels[i]].append(filenames[i])
	return flnames_cluster




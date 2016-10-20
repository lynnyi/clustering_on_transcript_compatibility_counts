# chooses subsamples of cells to quantify


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
idxpath = '/home/lynnyi/clustering/clustering_on_transcript_compatibility_counts/Trapnell_pipeline/kallisto_index/Trapnell_index.idx'
nprocesses=30
filenames = np.loadtxt('/home/lynnyi/clustering/clustering_on_transcript_compatibility_counts/Trapnell_pipeline/Trapnells_data/Trapnell_filenames.txt',dtype=str)
kallipso_path = '/home/lynnyi/clustering/clustering_on_transcript_compatibility_counts/modified-kallisto-source/kallisto_pseudo_paired/build/src/kalliPso'


base_dir = '/home/lynnyi/clustering/clustering_on_transcript_compatibility_counts/Trapnell_pipeline/'
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


# get names of cells within each cluster
def get_names_per_cluster():
	flnames_cluster = []
	for i in range(nclust):
		flnames_cluster.append([])
	for i in range(ncells):
		flnames_cluster[labels3[i]].append(filenames[i])
	return flnames_cluster


def make_directories():
	#os.system('mkdir ' + output_path)
	for i in range(nclust):
		i_directory = output_path+str(i)
		os.system('mkdir '+i_directory)


# functions for subsampling # divides total number of cells by ncells to determine the number of subsamples
def make_subsamples(nsubsamples):
	flnames_clust = get_names_per_cluster()
	subsamples = []
	for i in range(nclust):
		subsamples.append([])
		ncells = len(flnames_clust[i]) / nsubsamples
		flnames=flnames_clust[i]
		for j in range(nsubsamples):
			subsamples[i].append(flnames[j*ncells:(j+1)*ncells])
	return subsamples

# make subsamples randomly
def make_subsamples_random(ncells):
        flnames_clust = get_names_per_cluster()
	print(flnames_clust)
        subsamples = []
        for i in range(nclust):
                subsamples.append([])
                nsubsamples = len(flnames_clust[i]) / ncells
                flnames=flnames_clust[i]
                for j in range(nsubsamples):
			this_subsample = random.sample(flnames, ncells)
			subsamples[i].append(this_subsample)
			flnames = [x for x in flnames if x not in this_subsample]
        return subsamples

# make filetuples over all cells for kallisto to run
def make_kallisto_fltuples():
	fltuples = []
	for i in range(ncells):
		out_path = output_path+str(labels3[i])
		fltuples.append([filenames[i], idxpath, out_path])
	return fltuples


# make filetuples for subsamples of cells, returns list of list, with subsample array, subsample[i] corresponding to cluster i's subsamples -- each subsample[i][j] corresponds to filenames chosen for that subsample. first element of fltuple is list of filenames.
def make_kallisto_fltuples_from_subsamples(subsamples):
	fltuples = []
	for i in range(nclust):
		for j in range(len(subsamples[i])):
			out_path = output_path + 'cluster' + str(i)
			out_path = out_path + 'subsample' + str(j)
			fltuples.append([subsamples[i][j], idxpath, out_path])
	return fltuples

# centroids is just list of length 3 with filenames pertaining to each cluster
def make_kallisto_centroid_tupls(centroids):
	fltuples = []
	for i in range(nclust):
		out_path = output_path + str(i)
		fltuples.append([centroids[i], idxpath, out_path])
	return fltuples

# perform quantifications, map to a per_cell basis
def run_kallisto(fltuple):
	flname = fltuple[0]
	idxpath=fltuple[1]
	outpath=fltuple[2]+'/'+flname+'/'
	flname_1=path_to_reads + fltuple[0]+'_1.fastq.gz'
	flname_2=path_to_reads + fltuple[0]+'_2.fastq.gz'
	cmd = kallipso_path + ' quant -b 100 -i ' + idxpath + ' -o ' + outpath + ' ' + flname_1 + ' ' + flname_2
	os.system(cmd)

#input tupl is lsit of cells that kallisto will quantify simultaneously quantifying cluster without bootstraps
def run_kallisto_cluster(tupl):
        flname_clust = tupl[0]
        idxpath = tupl[1]
        outpath = tupl[2]
        cmd = kallipso_path + ' quant -i '
        cmd = cmd + idxpath + ' -o ' + outpath
        for j in range(len(flname_clust)):
                flname_1=path_to_reads + flname_clust[j]+'_1.fastq.gz'
                flname_2=path_to_reads + flname_clust[j]+'_2.fastq.gz'
                cmd = cmd + ' ' + flname_1 + ' ' + flname_2
        os.system(cmd)

# functions for quantifying cluster center
def run_kallisto_cluster_with_bootstraps(tupl):
	flname_clust = tupl[0]
	idxpath = tupl[1]
	outpath = tupl[2]
	cmd = kallipso_path + ' quant -b 100 -i ' + idxpath +  ' -o ' + outpath
	for j in range(len(flname_clust)):
        	flname_1=path_to_reads + flname_clust[j]+'_1.fastq.gz'
	        flname_2=path_to_reads + flname_clust[j]+'_2.fastq.gz'
		cmd = cmd + ' ' + flname_1 + ' ' + flname_2
	print(cmd)
	os.system(cmd)



"""
output_path = './singlecellquant2cluster'
make_directories()
tpls = make_kallisto_fltuples()
pool = mp.Pool(30)
pool.map(run_kallisto,tpls)
"""
"""
for i in [2, 5, 10]:
	ncellspersubsample = i
	output_path = './subsampled_' + str(i) +'cellsper/'
	os.system('mkdir ' + output_path)
	make_directories()
	subsamples = make_subsamples_random(ncellspersubsample)
	print(subsamples)
	tpls = make_kallisto_fltuples_from_subsamples(subsamples)
	pool = mp.Pool(20)
	pool.map(run_kallisto_cluster_with_bootstraps,tpls) 
"""
"""
output_path = './labels'+str(nclust)+'centroids'
make_directories()
centroids = get_names_per_cluster()
centroid_tpl = make_kallisto_centroid_tupls(centroids)
pool = mp.Pool(nclust)
pool.map(run_kallisto_cluster,centroid_tpl)
"""



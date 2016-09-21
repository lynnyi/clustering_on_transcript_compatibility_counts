import os
import getopt
import pickle
import itertools
import multiprocessing as mp
import numpy as np

nclust = 3
path_to_reads = '/home/vasilis/clustering_EQ/Trapnell/read_data/'
ncells = 271
idxpath = '/home/lynnyi/clustering/clustering_on_transcript_compatibility_counts/Trapnell_pipeline/kallisto_index/Trapnell_index.idx'
nprocesses=30
filenames = np.loadtxt('./Trapnells_data/Trapnell_filenames.txt',dtype=str)
kallipso_path = '/home/lynnyi/clustering/clustering_on_transcript_compatibility_counts/modified-kallisto-source/kallisto_pseudo_paired/build/src/kalliPso'


# Clustering
with open('./Trapnell_TCC_distribution.dat', 'rb') as infile:
    X = pickle.load(infile)
with open('./Trapnell_TCC_pairwise_distance.dat','rb') as infile:
    D = pickle.load(infile)
Trap_labels=np.loadtxt('./Trapnells_data/Trapnell_labels.txt',dtype=str)

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

# perform quantifications over all cells
def run_kallisto(fltuple):
	flname = fltuple[0]
	idxpath=fltuple[1]
	outpath=fltuple[2]+'/'+flname+'/'
	flname_1=path_to_reads + fltuple[0]+'_1.fastq.gz'
	flname_2=path_to_reads + fltuple[0]+'_2.fastq.gz'
	cmd = kallipso_path + ' quant -b 100 -i ' + idxpath + ' -o ' + outpath + ' ' + flname_1 + ' ' + flname_2
	os.system(cmd)


def make_directories():
	#os.system('mkdir ' + output_path)
	for i in range(nclust):
		i_directory = output_path+str(i)
		os.system('mkdir '+i_directory)

# make filetuples over all cells for kallisto to run
def make_kallisto_fltuples():
	fltuples = []
	for i in range(ncells):
		out_path = output_path+str(labels3[i])
		fltuples.append([filenames[i], idxpath, out_path])
	return fltuples

# functions for subsampling # not complete
def make_ncell_subsamples(ncells):
        flnames_clust = get_names_per_cluster()
	print(flnames_clust)
        subsamples = []
        for i in range(nclust):
                subsamples.append([])
                nsubsamples = len(flnames_clust[i]) / ncells
                flnames=flnames_clust[i]
                for j in range(nsubsamples):
                        subsamples[i].append(flnames[j*ncells:(j+1)*ncells])
        return subsamples

# functions for subsampling # not complete
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

# make filetuples for subsamples of cells
"""
def make_kallisto_fltuples(subsamples):
	fltuples = []
	uor i in range(nclust):
		for j in range(len(subsamples[i])):
			out_path = output_path + 'cluster' + str(i)
			out_path = out_path + 'subsample' + str(j)
			fltuples.append([subsamples[i][j], idxpath, out_path])
	return fltuples
"""
def make_kallisto_centroid_tupls(centroids):
	fltuples = []
	for i in range(nclust):
		out_path = output_path + str(i)
		fltuples.append([centroids[i], idxpath, out_path])
	return fltuples

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


#quantifying cluster without bootstraps
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

output_path = './singlecellquant2cluster'
make_directories()
tpls = make_kallisto_fltuples()
pool = mp.Pool(30)
pool.map(run_kallisto,tpls)

"""
ncellspersubsample = 20
output_path = './subsampled_20cellsper/'
make_directories()
subsamples = make_ncell_subsamples(ncellspersubsample)
print(subsamples)
subsample_tpl = make_kallisto_fltuples(subsamples)
pool = mp.Pool(20)
pool.map(run_kallisto_cluster_with_bootstraps,subsample_tpl) 

output_path = './labels'+str(nclust)+'centroids'
make_directories()
centroids = get_names_per_cluster()
centroid_tpl = make_kallisto_centroid_tupls(centroids)
pool = mp.Pool(nclust)
pool.map(run_kallisto_cluster,centroid_tpl)
"""



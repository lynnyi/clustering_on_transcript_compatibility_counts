
import numpy as np
import scipy as sp
import scipy.special.digamma as digamma


#load in the TCCs
execfile('sampling.py')
names = get_names_per_cluster()
execfile('tcc.py')
make_dictionary_of_TCCs(names[0])

#d will be the dictionary merged
#num-ECs
num_ECs = len(d.keys())

#num-samples
num_samples = num_samples

#log_p is the sum of the log of the kth component over all the samples
def calculate_logp(dicts)
	log_p = np.array([])
	return log_p

#gradient of likelihood at specific point of alpha
def gradient(alpha):
	alpha_sum = np.sum(alpha)	
	return grad

#given old alpha, update old alpha to give new alpha
def update(old_alpha):
	return new_alpha



# This file does clustering on Trapnell data then
# arranges Trapnell 1st cluster into dictionary
# of TCC counts indexed by TCC numbers
# Attempted to fit to Dirichlet distribution
# by gradient descent.
# Lynn Yi October 2016

import scipy.special as sp



#load in the TCCs
execfile('clustering.py')
labels = perform_clustering()
names = get_names_per_cluster(labels)
execfile('tcc.py')
TCCs = make_dictionary_of_TCCs(names[0])
TCCs = turn_dict_into_list(TCCs)

#d will be the dictionary merged
#num-ECs is 561662 for cluster 0 of trapnell data set
num_ECs = num_ECs

#num-samples
num_samples = num_samples

# todo: add the multiplicative factor
#log_p is the sum of the log of the kth component over all the samples
def calculate_log_p(ec_counts):
	log_ec_counts = [map(np.log, x) for x in ec_counts]
	log_p = np.array([np.sum(x) for x in log_ec_counts]) 
	#i think this part of the paper is incorrect
	#log_p /= num_samples
	return log_p

log_p = calculate_log_p(TCCs)
assert(len(log_p)==num_ECs)

#gradient of likelihood at specific point of alpha
def gradient(alpha):
	digamma_alpha_sum = sp.digamma(np.sum(alpha))
	alpha = np.array(log_p) - np.array([sp.digamma(x) for x in alpha])
	alpha += digamma_alpha_sum
	alpha *= num_samples
	return alpha

#given old alpha, update old alpha to give new alpha
def update(alpha, step_size):
	alpha += step_size*np.array(gradient(alpha))
	return alpha

#given specific alpha, get likelihood over the dirichlet
def log_likelihood(alpha):
        log_gamma_sum = sp.gammaln(np.sum(alpha))
	print('1 alpha sum')
	print(log_gamma_sum)
	sum_log_gamma = np.sum([sp.gammaln(x) for x in alpha])
        print('2 sum log gamma')
	print(sum_log_gamma)
	sum_dot_log_p = np.sum(np.dot([x-1 for x in alpha], log_p))
	print('3 sum_dot_log_p')
	print(sum_dot_log_p)
	likelihood = num_samples * (log_gamma_sum - sum_log_gamma + sum_dot_log_p)
	return likelihood
	
def gradient_descent(num_steps, step_size):
	alpha = np.ones(num_ECs)
	assert(len(log_p) == len(alpha))
	print(log_likelihood(alpha))
	for i in range(num_steps):
		alpha = update(alpha, step_size)
		print(log_likelihood(alpha))
	return alpha





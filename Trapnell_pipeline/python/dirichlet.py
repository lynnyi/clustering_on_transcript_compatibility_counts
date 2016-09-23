
import numpy as np
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

#log_p is the sum of the log of the kth component over all the samples
def calculate_log_p(ec_counts):
	log_ec_counts = [map(np.log, x) for x in ec_counts]
	log_p = map(np.sum, log_ec_counts)
	log_p /= num_samples
	return log_p

log_p = calculate_log_p(TCCs)
assert(len(log_p)==num_ECs)

#gradient of likelihood at specific point of alpha
def gradient(alpha):
	digamma_alpha_sum = sp.digamma(np.sum(alpha))
	alpha = np.array(log_p) - np.array([sp.digamma(x) for x in alpha])
	alpha += digamma_alpha_sum
	alpha *= num_samples
	#alpha = [(x + digamma_alpha_sum) * num_samples for x in alpha] 
	return alpha

#given old alpha, update old alpha to give new alpha
def update(alpha, step_size):
	alpha += step_size*np.array(gradient(alpha))
	return alpha

#given specific alpha, get likelihood over the dirichlet
def likelihood(alpha):
        log_gamma_sum = np.log(sp.gamma(np.sum(alpha)))
	sum_log_gamma = np.sum([np.log(sp.gamma(x)) for x in alpha])
        sum_dot_log_p = np.sum(np.dot([x-1 for x in alpha], log_p))
	likelihood = num_samples * (log_gamma_sum - sum_log_gamma + sum_dot_log_p)
	return likelihood
	
def gradient_descent(num_steps, step_size):
	alpha = np.random.rand(num_ECs)
	assert(len(log_p) == len(alpha))
	for i in range(num_steps):
		alpha = update(alpha, step_size)
		print likelihood(alpha)
	return alpha

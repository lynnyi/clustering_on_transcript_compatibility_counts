
# problem with my method is I don't keep track of cell counts or distributions
# better to stick with .dat file so nothing is removed
# may TCC_version2 is actually to modify the tcc
import numpy as np
import pandas as pd
matplotlib.use('Agg')
import matplotlib.pyplot as plt

base_dir = '/home/lynnyi/clustering/clustering_on_transcript_compatibility_counts/Trapnell_pipeline/'
directory_class_files = base_dir + '/transcript_compatibility_counts/'


def read_classfile_into_dictionary(filename):
	d = {}
	filename = directory_class_files + filename + '.class'
	print('reading:')
	print(filename)
	with open(filename) as f:
		for line in f:
			line = line.split()
			if len(line)!=2:
				print('error')
			ec_num = line[0]
			ec_count = line[1]
			d[ec_num] = int(ec_count)
	return d

def get_unique_keys(list_dict):
	keys = []
	for dic in list_dict:
		keys = keys + dic.keys()
	keys = np.unique(keys)
	return keys

def normalize_median(dictionary):
        values = dictionary.values()
        median = np.median(values)
        for key, value in dictionary.iteritems():
		dictionary[key]=float(value)/median
        return dictionary

def normalize_log_median(dictionary):
	values = dictionary.values()
	median = np.median(values)
	for key, value in dictionary.iteritems():
		dictionary[key]=np.log(float(value)/median+1.0)

def normalize_multinomial(dictionary):
	total = np.sum(dictionary.values())
	for key,value in dictionary.iteritems():
		dictionary[key] = float(value)/total
	return dictionary

def merge_dicts(dicts):
	d = {}
	for dic in dicts:
		for key in dic:
			try:
				d[key].append(dic[key])
			except KeyError:
				d[key] = [dic[key]]
	return d


global num_samples
global num_ECs

def make_dictionary_of_TCCs(filenames):
	x = map(read_classfile_into_dictionary, filenames)
	#print('normalizing counts')
	#x = map(normalize_median, x)
	d = merge_dicts(x)
	global num_samples
	global num_ECs
	num_samples = len(x)
	num_ECs = len(d.keys())
	return d

global df
def print_TCC_to_csv(filenames, filename):
	x = map(read_classfile_into_dictionary, filenames)
	global df
	print('making dataframe')
	df = pd.DataFrame(x)
	print('filling NaNs')
	df = df.fillna(value=0)
	df = df.astype(int)
	#assert(isinstance(df.iloc[[10],[11]], int))
	df.to_csv(filename, header=False)
	return

def turn_dict_into_list(dict_TCC):
	list_TCC = []
	for key in np.sort(dict_TCC.keys()):
		list_TCC.append(dict_TCC[key])
	assert(len(list_TCC) == num_ECs)
	return list_TCC


def plot_mean_var(d):
	mean = []
	var = []
	for key, value in d.iteritems():
		this_mean = np.sum(value) / num_samples
		this_var = np.var(list(value) + list(np.zeros(num_samples - len(value))))
		mean.append((this_mean))
		var.append((this_var))
	print('correl')
	print(np.corrcoef(mean,var))
	fig, ax = plt.subplots()
	ax.scatter(mean,var)
 	plt.savefig('mean_var_tcc_nonorm.png')




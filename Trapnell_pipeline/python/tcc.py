
import numpy as np
import pandas as pd
import os
import matplotlib
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
			d[ec_num] = float(ec_count)
	return d

def get_unique_keys(list_dict):
	keys = []
	for dic in list_dict:
		keys = keys + dic.keys()
	keys = np.unique(keys)
	return keys

def normalize_counts(dictionary):
        values = dictionary.values()
        median = np.median(values)
        if median == 0:
                print('median is 0')
        else:
                for key, value in dictionary.iteritems():
                        dictionary[key]=value/median
        return dictionary

def normalize_counts_with_log(dictionary):
	values = dictionary.values()
	median = np.median(values)
	for key, value in dictionary.iteritems():
		dictionary[key]=(value/median+1.0)
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
	print('normalizing counts')
	x = map(normalize_counts, x)
	d = merge_dicts(x)
	global num_samples
	global num_ECs
	num_samples = len(x)
	num_ECs = len(d.keys())
	return d


def plot_mean_var():
	mean = []
	var = []
	for key, value in d.iteritems():
		this_mean = np.sum(value) / num_samples
		this_var = np.var(list(value) + list(np.zeros(num_samples - len(value))))
	mean.append(np.log(this_mean))
	var.append(np.log(this_var))
	fig, ax = plt.subplots()
	ax.scatter(mean,var)
 	plt.savefig('mean_var_tcc_log_log_plus_one.png')




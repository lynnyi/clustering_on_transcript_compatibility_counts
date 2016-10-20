# reads Trapnell TCC counts into dictionary
# with normalizing options
# then plots mean variance plots
import numpy as np
import pandas as pd
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

base_dir = '/home/lynnyi/clustering/clustering_on_transcript_compatibility_counts/Trapnell_pipeline/'
directory = base_dir + '/subsampled_20cellsper/'
directory = '/home/lynnyi/clustering/clustering_on_transcript_compatibility_counts/Trapnell_pipeline/singlecellquant2/singlecellquant2cluster0/'

def read_counts_into_dictionary_from_file(filename):
	d = {}
	print('reading:')
        print(filename)
        with open(filename) as f:
                title = f.readline()
                for line in f:
                        line = line.split()
                        target_id = line[0]
                        est_counts = line[3]
                        d[target_id] = float(est_counts)
	return d

def read_counts_into_dictionary_from_folder(foldername):
	filename = directory + foldername + '/abundance.tsv'
	return read_counts_into_dictionary_from_file(filename)

def read_counts_into_dictionary_from_folderpath(folderpath):
	filename = folderpath + '/abundance.tsv'
	return read_counts_into_dictionary_from_file(filename)

def read_counts_into_dictionary_from_folders(parentfolder):
	folders = os.listdir(parentfolder)
	dicts = map(read_counts_into_dictionary_from_folder, folders)
	return dicts

def normalize_counts(dictionary):
	values = dictionary.values()
	median = np.median(values)
	if median == 0:
		print('median is 0')
	else:
		for key, value in dictionary.iteritems():
			dictionary[key]=np.log(value/median)
	return dictionary

def normalize_nonzero_counts(dictionary):
	values = np.nonzero(dictionary.values())
	median = np.median(values)
	if median == 0:
		print('median is 0')
	else:
		print('median is')
		print(median)
		for key,value in dictionary.iteritems():
			dictionary[key] = (value/median)
	return dictionary


def read_counts_in_folder_and_plot(directory):
	dictionaries = []
	#foldernames = glob.glob(directory + '/cluster0*')
	#foldernames = os.listdir(directory)
	foldernames = glob.glob(directory+'/*')
	print('number of files: ')
	print(len(foldernames))
	for foldername in foldernames:
		dictionaries.append(read_counts_into_dictionary_from_folderpath(foldername))
	map(normalize_nonzero_counts, dictionaries)
	data = pd.DataFrame(dictionaries)
	#mean-variance curve
	mean = np.array(list(data.mean()))
	var = np.array(list(data.var()))
	fig, ax = plt.subplots()
	ax.scatter(mean,var)
	#plt.xlim(0,5000)
	#plt.ylim(0,5000)
	#plt.gca().set_aspect('equal', adjustable='box')
	plt.savefig('meanvarsinglecells_nonzero_norm_cluster0.png')


"""
ax1 = fig.add_subplot(111)
ax1 = plt.scatter(mean,var)

##plotting
target_ids = list(data.columns.values)
ax1 = fig.add_subplot(221)
ax1 = plt.hist(data[target_ids[0]])
plt.savefig('test.png')
"""

"""
d = read_counts_into_dictionary_from_file(base_dir + '/singlecellquant2/singlecellquant2cluster0/SRR1033209/abundance.tsv')
f = read_counts_into_dictionary('./singlecellquant2/singlecellquant2cluster0/SRR1033210/abundance.tsv')

x = sorted(d)
y = sorted(f)

for i in range(len(x)):
	if (x[i] != y[i]):
		print('problem!')
normalize_counts(d)
normalize_counts(f)
"""



import numpy as np
import pandas as pd
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

directory = '/home/lynnyi/clustering/clustering_on_transcript_compatibility_counts/Trapnell_pipeline/singlecellquant2/singlecellquant2cluster0/'

def read_counts_into_dictionary(foldername):
	d = {}
	filename = directory + foldername + '/abundance.tsv'
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

def normalize_counts(dictionary):
	values = dictionary.values()
	median = np.median(values)
	for key, value in dictionary.iteritems():
		if median != 0:
			dictionary[key]=value/median

dictionaries = []
foldernames = os.listdir(directory)
print('number of files: ')
print(len(foldernames))

for foldername in foldernames:
	dictionaries.append(read_counts_into_dictionary(foldername))
map(normalize_counts, dictionaries)
data = pd.DataFrame(dictionaries)

#mean-variance curve
mean = np.array(list(data.mean()))
var = np.array(list(data.var()))

fig, ax = plt.subplots()
ax.plot(mean,var)
plt.xlim(0,5000)
plt.ylim(0,5000)
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig('meanvar.png')

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
d = read_counts_into_dictionary('./singlecellquant2/singlecellquant2cluster0/SRR1033209/abundance.tsv')
f = read_counts_into_dictionary('./singlecellquant2/singlecellquant2cluster0/SRR1033210/abundance.tsv')

x = sorted(d)
y = sorted(f)

for i in range(len(x)):
	if (x[i] != y[i]):
		print('problem!')
normalize_counts(d)
normalize_counts(f)
"""


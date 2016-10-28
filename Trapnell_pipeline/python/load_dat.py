
import pickle
import numpy as np
import csv
infile=open('./labels.dat','rb')
labels=pickle.load(infile)
labels = labels.tolist()

infile = open('../Trapnell_TCC.dat', 'rb')
load = pickle.load(infile)
print('converting to list')
df = load.toarray()

f = open('./TCC.txt', 'w')
writer = csv.writer(f)
for row in df:
	print(row)
	writer.writerow(row)


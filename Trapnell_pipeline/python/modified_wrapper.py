import os
import getopt
import sys
import numpy as np

num_proc = 30

"""
# Parse the options from command line
try:
    opts, args = getopt.getopt(sys.argv[1:],"i:n:k:t:",["idir=","njobs=","hacked-kallisto-path=","reference-transcriptome="])
except getopt.GetoptError:
    print ("getopterrror")
    print ('usage is : \n python Trapnell_wrapper.py -i input_SRA_dir -k path-to-hacked-kallisto -t path-to-human-reference-transcriptome [-n number-of-processes-to-use]')
    sys.exit(1)

# The four options used in this pipeline are    
SRA_dir=''
num_proc=1
kallipso_path=''
ref_transcriptome=''

for opt,arg in opts:
    if opt in ("-i", "--idir"):
        SRA_dir=arg
    elif opt in ("-n","--njobs"):
        num_proc=int(arg)
    elif opt in ("-k","--hacked-kallisto-path"):
        kallipso_path=arg
    elif opt in ("-t","--reference-transcriptome"):
        ref_transcriptome=arg
    
# Error message if function not called correctly        
if (not SRA_dir) or (not kallipso_path) or (not ref_transcriptome):
    print ('usage is : \n python Trapnell_wrapper.py -i input_SRA_dir -k path-to-hacked-kallisto -t path-to-human-reference-transcriptome [-n number-of-processes-to-use]')
    sys.exit(1)


# Execute process_SRA, which converts SRA to fastq
print('Extracting reads from SRAs...')
os.system('mkdir -p ./reads/')
os.system('rm -f ./reads/*')
os.system('python process_SRA.py -i '+SRA_dir+' -o ./reads/ -n '+str(num_proc))

# This removes files that Trapnell didn't use. (Specific for analysis of this dataset)
print('Removing files that Trapnell throws...')
relevant_fls=np.loadtxt('Files_to_keep.txt',dtype=str)
fls_to_remove=[cells for cells in os.listdir('./reads/') if cells.split('_')[0] not in relevant_fls ]
for fl in fls_to_remove:
    os.system('rm '+'./reads/'+fl)
"""


"""
ref_transcriptome = '~/clustering/Homo_sapiens.GRCh38.rel79.cdna.all.fa'
# Use kallisto to build index. should return same results as hacked kallisto
print('Generating the Kallisto index (with hacked kallisto)...')
os.system('mkdir -p ./kallisto_index')
os.system('rm -f ./kallisto_index/*')
index_path='./kallisto_index/Trapnell_index.idx'
os.system(kallipso_path+' index -i '+index_path+' '+ref_transcriptome)
metadata_cmd=kallipso_path+' metadata '+index_path
os.system(metadata_cmd)
num_ec = sum(1 for line in open('./kallisto_index/Trapnell_index.idx_ecmap.txt'))
print(num_ec)

kallipso_path =  '~/clustering/clustering_on_transcript_compatibility_counts/modified-kallisto-source/kallisto_pseudo_paired/build/src/kalliPso' 
index_path='./kallisto_index/Trapnell_index.idx'
path_to_reads = '/home/vasilis/clustering_EQ/Trapnell/read_data/'

# Use hacked kallisto to generate TCCs (output is .class file that contains counts of ECs)
print('Generating TCC (with hacked kallisto)...')
os.system('mkdir -p ./transcript_compatibility_counts/')
os.system('rm -f ./transcript_compatibility_counts/*')
os.system('python get_pseudoalignments_paired_end.py -i '+ path_to_reads +' -o ./transcript_compatibility_counts/ -k '+kallipso_path+ ' -t '+ index_path +' -n '+ str(num_proc))

"""

num_ec = sum(1 for line in open('./kallisto_index/Trapnell_index.idx_ecmap.txt'))
path_to_tcc = '/home/lynnyi/clustering/clustering_on_transcript_compatibility_counts/Trapnell_pipeline/TCC_2ndhalf/'

# Process the .class file into a matrix. normalize such that each row sums to 1 
print('Generating TCC distribution...')
os.system('python get_tcc_dist.py -i ' + path_to_tcc + ' -m '+str(num_ec)+' -t ./Trapnell_TCC_2ndhalf.dat -d ./Trapnell_TCC_2ndhalf_distribution.dat')

# Process the TCC distribution matrix to get pairwise distributions
print('Generating pairwise distances...')
os.system('python get_pairwise_distances.py ./Trapnell_TCC_2ndhalf_distribution.dat ./Trapnell_TCC_2ndhalf_pairwise_distance.dat '+str(num_proc))

"""

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

pref = -1.3*np.ones(271)
labels3=AffinityProp(-D,pref,0.95)

"""

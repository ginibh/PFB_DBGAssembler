#!/usr/bin/env python3

#from scipy.signal import argrelmin
import numpy as np
import sys
import matplotlib.pyplot as plt

filename = sys.argv[1]
kmer_length = int(sys.argv[2])

fileop = open(filename, 'r')

def kmer_count(filename, kmer_length):
   FILE = open(filename, 'r')
   if kmer_length % 2 == 0: raise ValueError('Kmer length cannot be an even number')
   kmers = {}
   for line in FILE:
    header = line
    seq = next(FILE)
    seq = seq.rstrip()
    plus = next(FILE)
    quality = next(FILE)
    #sliding window
    for i in range(0, len(seq) - kmer_length + 1):
       kmer = seq[i:i + kmer_length]
       if kmer in kmers:  kmers[kmer] += 1
       else:  kmers[kmer] = 1
   return(kmers)

kmers = kmer_count(filename, kmer_length)

kmers_sorted = sorted(kmers.items(), key = lambda x: x[1])

dist_kmer = []
for kmer in kmers:
	if kmers[kmer] > 5:
		dist_kmer.append(kmers[kmer])	

kmer_list = []
for kmer in kmers:
	kmer_list.append(kmers[kmer])
kmer_list = sorted(kmer_list)
kmer_count = {}

for i in range(min(kmer_list), max(kmer_list)+1):
	tempcount = kmer_list.count(i)
	kmer_count[i] = tempcount

minkmer = max(kmer_count.values())*100
value_greater_than_minkmer=[]
for key, value in sorted(kmer_count.items(), key = lambda x: x[0]):
#for i in range(1, max(kmer_count.keys())):
	print(key, value, minkmer)
	if value < minkmer:
		minkmer = value
	elif value > minkmer:
		print(key)
		value_greater_than_minkmer.append(value)
	else:
		print(key-1)
		break
print(value_greater_than_minkmer)
print(kmer_count)

#plt.hist(kmers.values())
#plt.plot()
#plt.savefig('graphforreals.png')
plt.hist(dist_kmer)
#plt.xlim([0,45])
#plt.ylim([0, 100000])
plt.plot()
plt.savefig('graphrealfreqy.png')
#print(sum(kmers.values())/len(kmers))


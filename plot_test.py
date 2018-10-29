#!/usr/bin/env python3
import matplotlib.pyplot as plt
import sys
filename=sys.argv[1]
kmer_length=int(sys.argv[2])

def kmer_count(filename, kmer_length):
	FILE=open(filename,'r')
	if kmer_length % 2 == 0: raise ValueError('Kmer length cannot be an even number')
	kmers={}
	for line in FILE:
	 header=line
	 seq=next(FILE)
	 seq=seq.rstrip()
	 plus=next(FILE)
	 quality=next(FILE)
	 for i in range(0,len(seq)-kmer_length+1):
	 	kmer=seq[i:i+kmer_length]
		if kmer in kmers: kmer[kmer]+=1
		else: kmers[kmer]=1
	return(kmers)
plt.hist(dict.values())
plt.plot()
plt.savefig(".png"

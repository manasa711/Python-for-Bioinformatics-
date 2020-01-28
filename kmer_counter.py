#!/usr/bin/env python3

import sys

def kmer(seq,k):
    count = {}
    list = []
    for i in range(len(seq)-k+1):
        x = seq[i:i+k]
        list.append(x)
    for j in list:
        if not (j in count):
            count[j] = 0
        count[j] += 1
    #print(count)
    for key in sorted(count.keys()):
        print(key+"\t"+str(count[key]))

k = sys.argv[1]

k = int(k)

fastafile = sys.argv[2]

fhand = open(fastafile)

seqs = ""

for line in fhand:
    if line.startswith('>'):
        continue
    line = line.rstrip()
    seqs = seqs + line

#print(seqs)

kmer(seqs,k)

#!/usr/bin/env python2.7
# It compares 2 genome files and make a image after it.


import unittest
from dnaseqlib import *

### Utility classes ###

# Maps integer keys to a set of arbitrary values.
class Multidict:
    # Initializes a new multi-value dictionary, and adds any key-value
    # 2-tuples in the iterable sequence pairs to the data structure.
    def __init__(self, pairs=[]):
        self.D=dict()
        for p in pairs:
            self.put(p[0],p[1])
    def put(self, k, v):
        try:
            self.D[k].append(v)
        except:
            self.D[k]=[v]
    # Gets any values that have been associated with the key k; or, if
    # none have been, returns an empty sequence.
    def get(self, k):
        try:
            return self.D[k]
        except:
            return []
# Given a sequence of nucleotides, return all k-length subsequences
# and their hashes.  (What else do you need to know about each
# subsequence?)
def getExactSubmatches(a, b, k, m):
    print "Making dictionary for seq a"
    Mydict=Multidict(intervalSubsequenceHashes(a,k,m))
    print "dictionary table done....."
    for hsh,(sequence,pos) in subsequenceHashes(b,k):
        for s,p in Mydict.get(hsh):
            if(sequence==s):
                yield (p,pos)
    return

def subsequenceHashes(seq, k):
    try:
        start= ""
        for i in range(k):
            start+=(seq.next())
        start2=RollingHash(start)
        pos=0
        yield (start2.current_hash(),(start,pos))
        while True:
            new = seq.next()
            start2.slide(start[0],new)
            start=start[1:]+new
            pos += 1
            yield (start2.current_hash(),(start,pos))
    except StopIteration:
        return

def intervalSubsequenceHashes(seq, k, m):
    try:
        pos = 0
        while True:
            answerseq = ''
            for i in range(k):
                answerseq += seq.next()
            for i in range(m-k):
                seq.next()
            rollseq = RollingHash(answerseq)
            yield (rollseq.current_hash(), (answerseq, pos))
            pos+=m
    except StopIteration:
        return

compareSequences(getExactSubmatches, 'my-shortmat-dog_8_100.png', (500,500), 'shortmat.fa', 'fdog0.fa', 8, 100)
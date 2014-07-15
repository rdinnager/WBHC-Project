import numpy
import math
import itertools
import sys

fastafile = sys.argv[1]
k = sys.argv[2]

def makeAllMers(k):
    allstr=[]
    for tempstr in itertools.product('ATGC', repeat=k):
        temp=''
        for i in range(0,k):
            temp=temp+tempstr[i]
        allstr.append(temp)
    return allstr

def makeAllMersUpTo(k):
    kmers = []
    ks = range(1,k+1)
    for kk in ks:
        kmers.append(makeAllMers(kk))
    return list(itertools.chain.from_iterable(kmers))

def computeRTDvectors(genome,mers):
    vector=[]
    for strtemp in mers:
        position=-1
        genomesequence=genome
        position=genomesequence.find(strtemp)
        listtemp=[]
        while position!=-1:
            listtemp.append(position)
            genomesequence=genomesequence[position+4:-1]
            position=genomesequence.find(strtemp)
        if len(listtemp)==0:
            vector.append(-1)
            vector.append(-1)
        else:
            vector.append(numpy.mean(listtemp))
            vector.append(math.sqrt(numpy.var(listtemp)))
    return vector

mers = makeAllMersUpTo(k)
mermean = []
mervar = []
for i in mers:
    mermean.append(str(i)+"_mean")
    mervar.append(str(i)+"_var")
pmers=list(itertools.chain.from_iterable(itertools.izip(mermean, mervar)))
print "SeqID,", ', '.join(pmers)
genomesequence=''
vector=[]
first=1
f=open(fastafile)
for line in f:
    line.replace('\n','')
    line.replace('\r','')
    if len(line)>1:
        if line[0]=='>':
            name = line[1:-1]
            if first==1:
                print str(name)+",",
                first=0
            if genomesequence!='':
                vector=computeRTDvectors(genomesequence,mers)
                print str(vector)[1:-1]
                if first==0:
                    print str(name)+",",
                vector=[]
                genomesequence=''
        else:
            genomesequence=genomesequence+line.upper()
f.close()
if genomesequence!='':
    vector=computeRTDvectors(genomesequence,mers)
    print str(vector)[1:-1]
    vector=[]
    genomesequence=''


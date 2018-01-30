import PyInstEvodb
from PyInstEvodb import *

import numpy as np
import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import copy

def mlfun_ex(x, Obj):
    '''Take a language tree object, set the parameters, unpack them
       and then calculate the entire likelihood'''

    Obj.parameters = np.matrix(x)
    Obj.unpack()
    Obj.settimes()
    L1=-Obj.SplitLikelihood()
    L2=-Obj.jukescantorlikelihood()
    L3=-Obj.DeathLikelihood()
    L4=-np.max(Obj.OriginLikelihood())
    return -L1-L2-L3-L4

# Import all the routines that I have written.
    
Data   = pd.read_pickle(os.getcwd() + '//IEData//MasterData.pkl')
Splits = pd.read_pickle(os.getcwd() + '//IEData//Splits.pkl')
Depths = pd.read_pickle(os.getcwd() + '//IEData//Depths.pkl')

NaDeneRT    = PyInstEvodb.ResolvedTree(Data.loc[Data['ruhlen_1']=='NADENE'],'NTree1')
numbranches = NaDeneRT.interiorbranches
bInit       = np.matrix(-1-np.linspace(0,10,num=numbranches)/numbranches)
rInit=np.zeros((1,len(NaDeneRT.words)))
dparms=np.sum(NaDeneRT.deathmat[:,0]==0)
dInit=np.zeros((1,dparms))+1
eInit=np.matrix(5)
parmsInit=np.hstack((bInit,rInit,dInit,eInit))

NaDenePT=PyInstEvodb.ParameterizedTree(Data.loc[Data['ruhlen_1']=='NADENE'], 'NTree1', parmsInit)


min = np.array(Depths['min'].loc[Depths['phylum'] == 'NaDene'])
max = np.array(Depths['max'].loc[Depths['phylum'] == 'NaDene'])
NaDenePT.priordepth(min[0], max[0])
NaDenePT.splitinfo(Splits[Splits['phylum'] == 'NaDene'])

NaDenePT.settimes()

VInit=np.eye((np.shape(parmsInit)[1]))

x, y, z = PyInstEvodb.myMcMcSampler_mwg(PyInstEvodb.mlfun, np.matrix(parmsInit), VInit, 100, 10, .5, .28, NaDenePT)


# Values where the likelihood starts throwing an error.

x1 = x[39]
x2 = x[40]

# culprit

NaDenePT.parameters = np.matrix(x2)
NaDenePT.unpack()
NaDenePT.settimes()

NaDenePT.OriginLikelihood()

# Note that the large values are throwing the errors. 




### Basic walkthrough of the code

TM = np.copy(NaDenePT.resolvedtree[:,-1])   #Just to shorten the code a bit
bp = NaDenePT.branchpositions      #Likewise
LL = np.zeros((rows(TM), 1))
IndCount = np.zeros((rows(TM), 1))
Live = np.zeros((rows(TM), 1))        

# Get rid of zeros in the Distance matrix:
D = np.copy(NaDenePT.D)
D[D <= 0] = 10   
lnD = np.matrix(-np.log(D))
np.fill_diagonal(lnD,0) 

# Beginning of loop

for i in range(rows(NaDenePT.resolvedtree), rows(bp)):     
    print(i)
    id = NaDenePT.resolvedtree[bp[i,0], bp[i,1]]
    tu = np.where(NaDenePT.resolvedtree[:,bp[i,1]] == id)[0]
    Bhat = NaDenePT.depth*1000*np.copy(NaDenePT.filledtimeFractions[tu,bp[i,1]:])
    TreeHat = np.copy(NaDenePT.resolvedtree[tu,bp[i,1]:])
    IndCountHat = np.copy(IndCount[tu])
    LLHat = np.copy(LL[tu])
    
    z = 0
    while True:
        ids = uniq(TreeHat[:,z])
        nums = TreeHat[:,z]
        z += 1
        if len(ids) > 1:
            break

    TMHat = np.zeros((len(nums), 1))
    IMHat = np.copy(IndCountHat)   
    DHat = (lnD[tu,:])[:,tu] 
    for m in ids:
        posi = np.where(nums == m)[0]
        posni = np.where(nums != m)[0]
        toAdd = np.matrix(np.nansum(Bhat[posni,:], axis=1)).T
        for q in posi:
            maxFinder = DHat[q,posni].T
            for n in range(0, len(posni)):
                if IndCountHat[posni[n]] >= 1:
                    maxFinder[n]=maxFinder[n] + IndCountHat[posni[n]]*(np.log(IndCountHat[posni[n]]) - np.log(toAdd[n])) + LLHat[posni[n]]
            max = np.argmax(maxFinder)           # No need to do more - defaults to first entry if more than one
            TMHat[q] = toAdd[max]
            IMHat[q] =1 + IndCountHat[posni[max]]
            LLHat[q] = DHat[q,posni[max]]
        TMHat[posi] = toAdd[max]
        IMHat[posi] = IndCountHat[posni[max]] + 1
        
    for p in range(0, rows(tu)):
        if Live[tu[p]] == 1:
            LL[tu[p]] = LL[tu[p]] + IndCount[tu[p]]*(np.log(IndCount[tu[p]]) - np.log(TM[tu[p]])) + LLHat[p]
            TM[tu[p]] = np.copy(TMHat[p])
            IndCount[tu[p]] = np.copy(IMHat[p])
        else:
            LL[tu[p]] = LLHat[p]
            Live[tu[p]] = 1
            TM[tu[p]] = TMHat[p]
            IndCount[tu[p]] = IMHat[p] 

# Loop - step by step until the first error is thrown:
            
TM = np.copy(NaDenePT.resolvedtree[:,-1])   #Just to shorten the code a bit
bp = NaDenePT.branchpositions      #Likewise
LL = np.zeros((rows(TM), 1))
IndCount = np.zeros((rows(TM), 1))
Live = np.zeros((rows(TM), 1))        

# Get rid of zeros in the Distance matrix:
D = np.copy(NaDenePT.D)
D[D <= 0] = 10   
lnD = np.matrix(-np.log(D))
np.fill_diagonal(lnD,0)             
            
for i in range(21, 22):     
    print(i)
    id = NaDenePT.resolvedtree[bp[i,0], bp[i,1]]
    tu = np.where(NaDenePT.resolvedtree[:,bp[i,1]] == id)[0]
    Bhat = NaDenePT.depth*1000*np.copy(NaDenePT.filledtimeFractions[tu,bp[i,1]:])
    TreeHat = np.copy(NaDenePT.resolvedtree[tu,bp[i,1]:])
    IndCountHat = np.copy(IndCount[tu])
    LLHat = np.copy(LL[tu])
    
    z = 0
    while True:
        ids = uniq(TreeHat[:,z])
        nums = TreeHat[:,z]
        z += 1
        if len(ids) > 1:
            break

    TMHat = np.zeros((len(nums), 1))
    IMHat = np.copy(IndCountHat)   
    DHat = (lnD[tu,:])[:,tu] 
    for m in ids:
        posi = np.where(nums == m)[0]
        posni = np.where(nums != m)[0]
        toAdd = np.matrix(np.nansum(Bhat[posni,:], axis=1)).T
        for q in posi:
            maxFinder = DHat[q,posni].T
            for n in range(0, len(posni)):
                if IndCountHat[posni[n]] >= 1:
                    maxFinder[n]=maxFinder[n] + IndCountHat[posni[n]]*(np.log(IndCountHat[posni[n]]) - np.log(toAdd[n])) + LLHat[posni[n]]
            max = np.argmax(maxFinder)           # No need to do more - defaults to first entry if more than one
            TMHat[q] = toAdd[max]
            IMHat[q] =1 + IndCountHat[posni[max]]
            LLHat[q] = DHat[q,posni[max]]
        TMHat[posi] = toAdd[max]
        IMHat[posi] = IndCountHat[posni[max]] + 1
        
    for p in range(0, rows(tu)):
        if Live[tu[p]] == 1:
            LL[tu[p]] = LL[tu[p]] + IndCount[tu[p]]*(np.log(IndCount[tu[p]]) - np.log(TM[tu[p]])) + LLHat[p]
            TM[tu[p]] = np.copy(TMHat[p])
            IndCount[tu[p]] = np.copy(IMHat[p])
        else:
            LL[tu[p]] = LLHat[p]
            Live[tu[p]] = 1
            TM[tu[p]] = TMHat[p]
            IndCount[tu[p]] = IMHat[p] 
            
for i in range(22, 23):     
    print(i)
    id = NaDenePT.resolvedtree[bp[i,0], bp[i,1]]
    tu = np.where(NaDenePT.resolvedtree[:,bp[i,1]] == id)[0]
    Bhat = NaDenePT.depth*1000*np.copy(NaDenePT.filledtimeFractions[tu,bp[i,1]:])
    TreeHat = np.copy(NaDenePT.resolvedtree[tu,bp[i,1]:])
    IndCountHat = np.copy(IndCount[tu])
    LLHat = np.copy(LL[tu])
    
    z = 0
    while True:
        ids = uniq(TreeHat[:,z])
        nums = TreeHat[:,z]
        z += 1
        if len(ids) > 1:
            break

    TMHat = np.zeros((len(nums), 1))
    IMHat = np.copy(IndCountHat)   
    DHat = (lnD[tu,:])[:,tu] 
    for m in ids:
        posi = np.where(nums == m)[0]
        posni = np.where(nums != m)[0]
        toAdd = np.matrix(np.nansum(Bhat[posni,:], axis=1)).T
        for q in posi:
            maxFinder = DHat[q,posni].T
            for n in range(0, len(posni)):
                if IndCountHat[posni[n]] >= 1:
                    maxFinder[n]=maxFinder[n] + IndCountHat[posni[n]]*(np.log(IndCountHat[posni[n]]) - np.log(toAdd[n])) + LLHat[posni[n]]
            max = np.argmax(maxFinder)           # No need to do more - defaults to first entry if more than one
            TMHat[q] = toAdd[max]
            IMHat[q] =1 + IndCountHat[posni[max]]
            LLHat[q] = DHat[q,posni[max]]
        TMHat[posi] = toAdd[max]
        IMHat[posi] = IndCountHat[posni[max]] + 1
        
    for p in range(0, rows(tu)):
        if Live[tu[p]] == 1:
            LL[tu[p]] = LL[tu[p]] + IndCount[tu[p]]*(np.log(IndCount[tu[p]]) - np.log(TM[tu[p]])) + LLHat[p]
            TM[tu[p]] = np.copy(TMHat[p])
            IndCount[tu[p]] = np.copy(IMHat[p])
        else:
            LL[tu[p]] = LLHat[p]
            Live[tu[p]] = 1
            TM[tu[p]] = TMHat[p]
            IndCount[tu[p]] = IMHat[p] 

# So we have found the culprit - in TM we are trying to take the log of zero. We even pointed this out
            # as a potential problem...
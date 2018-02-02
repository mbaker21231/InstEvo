


import numpy as np
import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import copy

d = os.getcwd()

os.chdir(d + '\\Documents\\GitHub\\InstEvo\\')


import PyInstEvo
from PyInstEvo import *


# Import all the routines that I have written.
    
Data   = pd.read_pickle(os.getcwd() + '//IEData//MasterData.pkl')
Splits = pd.read_pickle(os.getcwd() + '//IEData//Splits.pkl')
Depths = pd.read_pickle(os.getcwd() + '//IEData//Depths.pkl')

EskimoALRT    = PyInstEvo.ResolvedTree(Data.loc[Data['ruhlen_1']=='ESKIMOAL'],'NTree1')
numbranches = EskimoALRT.interiorbranches
bInit       = np.matrix(-1-np.linspace(0,10,num=numbranches)/numbranches)
rInit=np.zeros((1,len(EskimoALRT.words)))
dparms=np.sum(EskimoALRT.deathmat[:,0]==0)
dInit=np.zeros((1,dparms))+1
eInit=np.matrix(5)
parmsInit=np.hstack((bInit,rInit,dInit,eInit))

EskimoPT=PyInstEvo.ParameterizedTree(Data.loc[Data['ruhlen_1']=='ESKIMOAL'], 'NTree1', parmsInit)


min = np.array(Depths['min'].loc[Depths['phylum'] == 'EskimoAl'])
max = np.array(Depths['max'].loc[Depths['phylum'] == 'EskimoAl'])
EskimoPT.priordepth(min[0], max[0])
EskimoPT.splitinfo(Splits[Splits['phylum'] == 'EkimoAl'])

EskimoPT.settimes()
EskimoPT.showtree()

VInit=np.eye((np.shape(parmsInit)[1]))

x, y, z = PyInstEvo.myMcMcSampler_mwg(PyInstEvo.mlfun, np.matrix(parmsInit), VInit, 100, 10, .5, .28, EskimoPT)


# Values where the likelihood starts throwing an error.

EskimoPT.OriginLikelihood()

# Note that the large values are throwing the errors. 


   # Initialize placeholders:
TM = np.copy(EskimoPT.resolvedtree[:,-1])   #Just to shorten the code a bit
bp = EskimoPT.branchpositions      #Likewise
LL = np.zeros((rows(TM), 1))
IndCount = np.zeros((rows(TM), 1))
Live = np.zeros((rows(TM), 1))        

        # Get rid of zeros in the Distance matrix:
D = np.copy(EskimoPT.D)
D[D <= 0] = 10   
lnD = np.matrix(-np.log(D))
np.fill_diagonal(lnD,0) 

    # The main recursive loop backwards through the branches:
#for i in range(rows(self.resolvedtree), rows(bp)):   
    # i runs from 9 to 15...let's pop first branch and go with it
    
    i = 9
    id = EskimoPT.resolvedtree[bp[i,0], bp[i,1]]
    tu = np.where(EskimoPT.resolvedtree[:,bp[i,1]] == id)[0]
    Bhat = EskimoPT.depth*1000*np.copy(EskimoPT.filledtimeFractions[tu,bp[i,1]:])
    TreeHat = np.copy(EskimoPT.resolvedtree[tu,bp[i,1]:])
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
           LL[tu[p]] = LL[tu[p]] + IndCount[tu[p]]*(np.log(IndCount[tu[p]]) - np.log(TMHat[p])) + LLHat[p]
           TM[tu[p]] = np.copy(TMHat[p])
           IndCount[tu[p]] = np.copy(IMHat[p])
        else:
            LL[tu[p]] = LLHat[p]
            Live[tu[p]] = 1
            TM[tu[p]] = TMHat[p]
            IndCount[tu[p]] = IMHat[p]   
    

##### New function that picks the maximum-likelihood path recursively. 
##### Exponential version of the model...            
##### First, initialize everything we need:
            
bp       = EskimoPT.branchpositions                          #branch positions
bp       = bp[rows(TM):]                                     #interior branch only
branches = EskimoPT.filledtimeFractions*EskimoPT.depth*1000  #branch lengths
TM       = EskimoPT.resolvedtree

LL    = np.zeros((rows(TM), 1))                       #initial Log likelihood
NN    = np.zeros((rows(TM), 1))                       #Counter for number of jumps
TT    = np.zeros((rows(TM), 1))                       #Holder for branch lengths

Live = np.zeros((rows(TM), 1))                     #Location is live indicator        

        # Get rid of zeros in the Distance matrix:
D = np.copy(EskimoPT.D)
D[D <= 0] = 10   
lnD = np.matrix(-np.log(D))
np.fill_diagonal(lnD,0)           

#for r,c in bp:
    r,c   = bp[0,0], bp[0,1]
    id    = TM[r,c]
    tu    = np.where(TM[:, c] == id)[0]
    bhat  = branches[tu, c:]
    THat  = TM[tu, c:]
    NNHat = NN[tu]
    LLHat = LL[tu]
    TTHat = np.zeros((len(nums), 1))
    DHat  = (lnD[tu, :])[:, tu]
    
    z = 0
    while True:
        ids  = uniq(THat[:, z])
        nums = THat[:, z]
        z    = z + 1
        if len(ids) > 1:
            break
    
    for m in ids:
        posi  = np.where(nums == m)[0]
        posni = np.where(nums != m)[0]
        toAdd = np.matrix(np.nansum(bhat[posni, :], axis=1)).T
        for q in posi:
            maxFinder = DHat[q, posni].T
            for n in range(0, len(posni)):
                if NNHat[posni[n]] >= 1:
                    maxFinder[n] = maxFinder[n] + NNHat[posni[n]]*(np.log(NNHat[posni[n]]) - np.log(toAdd[n])) +  LLHat[posni[n]]
            max = np.argmax(maxFinder)                
            TTHat[q] = toAdd[max]
            NNHat[q] = NNHat[q] + 1
            LLHat[q] = DHat[q, posni[max]]
            
    for p in range(0, rows(tu)):
        if Live[tu[p]] == 1:
            LL[tu[p]] = LL[tu[p]] + NN[tu[p]]*(np.log(NN[tu[p]]) - np.log(TTHat[p])) + LLHat[p]
            TT[tu[p]] = TTHat[p]
            NN[tu[p]] = NNHat[p]
        else:
            LL[tu[p]]   = LLHat[p]
            Live[tu[p]] = 1
            TT[tu[p]]   = TTHat[p]
            NN[tu[p]]   = NNHat[p]
        
# Iteration #2

    r,c   = bp[1,0], bp[1,1]
    id    = TM[r,c]
    tu    = np.where(TM[:, c] == id)[0]
    bhat  = branches[tu, c:]
    THat  = TM[tu, c:]
    NNHat = NN[tu]
    LLHat = LL[tu]
    TTHat = np.zeros((len(nums), 1))
    DHat  = (lnD[tu, :])[:, tu]
    
    z = 0
    while True:
        ids  = uniq(THat[:, z])
        nums = THat[:, z]
        z    = z + 1
        if len(ids) > 1:
            break
    
    for m in ids:
        posi  = np.where(nums == m)[0]
        posni = np.where(nums != m)[0]
        toAdd = np.matrix(np.nansum(bhat[posni, :], axis=1)).T
        for q in posi:
            maxFinder = DHat[q, posni].T
            for n in range(0, len(posni)):
                if NNHat[posni[n]] >= 1:
                    maxFinder[n] = maxFinder[n] + NNHat[posni[n]]*(np.log(NNHat[posni[n]]) - np.log(toAdd[n])) +  LLHat[posni[n]]
            max = np.argmax(maxFinder)                
            TTHat[q] = toAdd[max]
            NNHat[q] = NNHat[q] + 1
            LLHat[q] = DHat[q, posni[max]]
            
    for p in range(0, rows(tu)):
        if Live[tu[p]] == 1:
            LL[tu[p]] = LL[tu[p]] + NN[tu[p]]*(np.log(NN[tu[p]]) - np.log(TTHat[p])) + LLHat[p]
            TT[tu[p]] = TTHat[p]
            NN[tu[p]] = NNHat[p]
        else:
            LL[tu[p]]   = LLHat[p]
            Live[tu[p]] = 1
            TT[tu[p]]   = TTHat[p]
            NN[tu[p]]   = NNHat[p]
        
# Third iteration
            
    r,c   = bp[2,0], bp[2,1]
    id    = TM[r,c]
    tu    = np.where(TM[:, c] == id)[0]
    bhat  = branches[tu, c:]
    THat  = TM[tu, c:]
    NNHat = NN[tu]
    LLHat = LL[tu]
    TTHat = np.zeros((len(nums), 1))
    DHat  = (lnD[tu, :])[:, tu]
    
    z = 0
    while True:
        ids  = uniq(THat[:, z])
        nums = THat[:, z]
        z    = z + 1
        if len(ids) > 1:
            break
    
    for m in ids:
        posi  = np.where(nums == m)[0]
        posni = np.where(nums != m)[0]
        toAdd = np.matrix(np.nansum(bhat[posni, :], axis=1)).T
        for q in posi:
            maxFinder = DHat[q, posni].T
            for n in range(0, len(posni)):
                if NNHat[posni[n]] >= 1:
                    maxFinder[n] = maxFinder[n] + NNHat[posni[n]]*(np.log(NNHat[posni[n]]) - np.log(toAdd[n])) +  LLHat[posni[n]]
            max = np.argmax(maxFinder)                
            TTHat[q] = toAdd[max]
            NNHat[q] = NNHat[posni[max]] + 1
            LLHat[q] = DHat[q, posni[max]]
            
    for p in range(0, rows(tu)):
        if Live[tu[p]] == 1:
            LL[tu[p]] = LL[tu[p]] + NN[tu[p]]*(np.log(NN[tu[p]]) - np.log(TTHat[p])) + LLHat[p]
            TT[tu[p]] = TTHat[p]
            NN[tu[p]] = NNHat[p]
        else:
            LL[tu[p]]   = LLHat[p]
            Live[tu[p]] = 1
            TT[tu[p]]   = TTHat[p]
            NN[tu[p]]   = NNHat[p]            
                
             
# Use locs flag

use_distance = False
poisson      = True
            
bp       = EskimoPT.branchpositions                          #branch positions
bp       = bp[rows(TM):]                                     #interior branch only
branches = EskimoPT.filledtimeFractions*EskimoPT.depth*1000  #branch lengths
TM       = EskimoPT.resolvedtree

LL    = np.zeros((rows(TM), 1))                       #initial Log likelihood
NN    = np.zeros((rows(TM), 1))                       #Counter for number of jumps
TT    = np.zeros((rows(TM), 1))                       #Holder for branch lengths

Live = np.zeros((rows(TM), 1))                     #Location is live indicator        

        # Get rid of zeros in the Distance matrix:

if use_distance:
    D = np.copy(EskimoPT.D)
    D[D <= 0] = 10   
    lnD = np.matrix(-np.log(D))
    np.fill_diagonal(lnD,0) 
else:
    lnD = np.zeros(np.shape(EskimoPT.D))          

for bb in bp:
    r,c   = bb[0], bb[1]
    id    = TM[r,c]
    tu    = np.where(TM[:, c] == id)[0]
    bhat  = branches[tu, c:]
    THat  = TM[tu, c:]
    NNHat = NN[tu]
    LLHat = LL[tu]

    DHat  = (lnD[tu, :])[:, tu]

    z = 0
    while True:
        ids  = uniq(THat[:, z])
        nums = THat[:, z]
        z    = z + 1
        if len(ids) > 1:
            break

    TTHat = np.zeros((len(nums), 1))
    
    for m in ids:
        posi  = np.where(nums == m)[0]
        posni = np.where(nums != m)[0]
        toAdd = np.matrix(np.nansum(bhat[posni, :], axis=1)).T
        if poisson:
            toAdd[:, :] = 1
        for q in posi:
            maxFinder = DHat[q, posni].T
            for n in range(0, len(posni)):
                if NNHat[posni[n]] >= 1:
                    maxFinder[n] = maxFinder[n] + NNHat[posni[n]]*(np.log(NNHat[posni[n]]) - np.log(toAdd[n])) +  LLHat[posni[n]]
            max = np.argmax(maxFinder)
            TTHat[q] = toAdd[max]
            NNHat[q] = NNHat[q] + 1
            LLHat[q] = DHat[q, posni[max]]
            
    for p in range(0, rows(tu)):
        if Live[tu[p]] == 1:
            LL[tu[p]] = LL[tu[p]] + NN[tu[p]]*(np.log(NN[tu[p]]) - np.log(TTHat[p])) + LLHat[p]
            TT[tu[p]] = TTHat[p]
            NN[tu[p]] = NNHat[p]
        else:
            LL[tu[p]]   = LLHat[p]
            Live[tu[p]] = 1
            TT[tu[p]]   = TTHat[p]
            NN[tu[p]]   = NNHat[p]

# Last thing to do is to go through and collect the garbage...
     

    























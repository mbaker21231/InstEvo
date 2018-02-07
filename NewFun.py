


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

   # Initialize placeholders:

TM       = EskimoPT.resolvedtree 
bp       = EskimoPT.branchpositions      
bp       = bp[rows(TM):]
branches = EskimoPT.filledtimeFractions*EskimoPT.depth*1000

LL       = np.zeros(rows(TM))
LL_D     = np.zeros(rows(TM))
NN       = np.zeros(rows(TM))        
TT       = np.zeros(rows(TM))
DD       = np.zeros(rows(TM))
Live     = np.zeros(rows(TM))
        
D         = np.copy(EskimoPT.D)
np.fill_diagonal(D, 1)
lnD       = - np.log(D)

    # The main recursive loop backwards through the branches:

#for bb in bp:
#   bb      = bp[0,:]
#   bb      = bp[1,:]
    bb      = bp[2,:]
    bb      = bp[3,:]
    bb      = bp[4,:]
    bb      = bp[5,:]
    bb      = bp[6,:]
    bb      = bp[7,:]
    r, c    = bb[0], bb[1]
    id      = TM[r, c]
    tu      = np.where(TM[:, c] == id)[0]
    bhat    = branches[tu, c:]
    THat    = TM[tu, c:]
    NNCount = (np.isnan(bhat) == False).sum(axis=1) - 1
    LLHat   = LL[tu]
    for p in tu:
        if Live[p] == 1:
            LL[p]   = LL[p] + NN[p] * (np.log(NN[p]) - np.log(TT[p]))
            LL_D[p] = LL_D[p] + DD[p]
     
    LL_DHat = LL_D[tu]
     
    DHat    = (lnD[tu, :])[:, tu]
     
    z = 0
    while True:
        ids  = uniq(THat[:, z])
        nums = THat[:, z]
        z    = z + 1
        if len(ids) > 1:
            break
         
    TTHat = np.zeros(len(nums))
    NNHat = np.zeros(len(nums))
    DDHat = np.zeros(len(nums))
     
    for m in ids:
        posi   = np.where(nums == m)[0]
        posni  = np.where(nums != m)[0]
        toAdd  =np.nansum(bhat[posni, :], axis=1)
        for q in posi:
            maxFinder = NNCount[posni]*(np.log(NNCount[posni]) \
                               - np.log(toAdd)) + LL_DHat[posni] + DHat[q, posni]
            maxm      = np.argmax(maxFinder)
            TTHat[q]  = toAdd[maxm]
            NNHat[q]  = NNCount[posni[maxm]]
            DDHat[q]  = DHat[q, posni[maxm]] + LL_DHat[posni[maxm]]
     
    for p in range(0, rows(tu)):
        TT[tu[p]] = TTHat[p]
        NN[tu[p]] = NNHat[p]
        DD[tu[p]] = DDHat[p]
        if Live[tu[p]] == 0:
            Live[tu[p]] = 1
    
    # Everyone should be live, so now we can just collect everything
    
for p in tu:
    if Live[p] == 1:
        LL[p]   = LL[p] + NN[p] * (np.log(NN[p]) - np.log(TT[p]))
        LL_D[p] = LL_D[p] + DD[p]    
    
    
     
def pr(thing):
    for t in thing:
        print(t)


def OriginLikelihood(Tree):
    
    TM       = Tree.resolvedtree 
    bp       = Tree.branchpositions      
    bp       = bp[rows(TM):]
    branches = Tree.filledtimeFractions*Tree.depth*1000
 
    LL       = np.zeros(rows(TM))
    LL_D     = np.zeros(rows(TM))
    NN       = np.zeros(rows(TM))        
    TT       = np.zeros(rows(TM))
    DD       = np.zeros(rows(TM))
    Live     = np.zeros(rows(TM))
        
    D         = np.copy(Tree.D)/1000
    np.fill_diagonal(D, 1)
    lnD       = - np.log(D)

    for bb in bp:
        r, c    = bb[0], bb[1]
        id      = TM[r, c]
        tu      = np.where(TM[:, c] == id)[0]
        bhat    = branches[tu, c:]
        THat    = TM[tu, c:]
        NNCount = (np.isnan(bhat) == False).sum(axis=1) - 1
        LLHat   = LL[tu]
        for p in tu:
            if Live[p] == 1:
                LL[p]   = LL[p] + NN[p] * (np.log(NN[p]) - np.log(TT[p]))
                LL_D[p] = LL_D[p] + DD[p]
     
        LL_DHat = LL_D[tu]
     
        DHat    = (lnD[tu, :])[:, tu]
     
        z = 0
        while True:
            ids  = uniq(THat[:, z])
            nums = THat[:, z]
            z    = z + 1
            if len(ids) > 1:
                break
         
        TTHat = np.zeros(len(nums))
        NNHat = np.zeros(len(nums))
        DDHat = np.zeros(len(nums))
     
        for m in ids:
            posi   = np.where(nums == m)[0]
            posni  = np.where(nums != m)[0]
            toAdd  =np.nansum(bhat[posni, :], axis=1)
            for q in posi:
                maxFinder = NNCount[posni]*(np.log(NNCount[posni]) \
                    - np.log(toAdd)) + LL_DHat[posni] + DHat[q, posni]
                maxm      = np.argmax(maxFinder)
                TTHat[q]  = toAdd[maxm]
                NNHat[q]  = NNCount[posni[maxm]]
                DDHat[q]  = DHat[q, posni[maxm]] + LL_DHat[posni[maxm]]
     
        for p in range(0, rows(tu)):
            TT[tu[p]] = TTHat[p]
            NN[tu[p]] = NNHat[p]
            DD[tu[p]] = DDHat[p]
            if Live[tu[p]] == 0:
                Live[tu[p]] = 1
    
    for p in tu:
        if Live[p] == 1:
            LL[p]   = LL[p] + NN[p] * (np.log(NN[p]) - np.log(TT[p]))
            LL_D[p] = LL_D[p] + DD[p]    
        
    return LL_D + LL
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
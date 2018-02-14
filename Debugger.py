import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import mplleaflet
import re

import PyInstEvo

from scipy.optimize import minimize

Data   = pd.read_pickle(os.getcwd() + '//IEData//MasterData.pkl')
Splits = pd.read_pickle(os.getcwd() + '//IEData//Splits.pkl')
Depths = pd.read_pickle(os.getcwd() + '//IEData//Depths.pkl')

Data['ex_date'].loc[Data['name'] == 'NAMA'] = 1500
Data['ex_date'].loc[Data['name'] == 'GWI']  = 100
Data['deadOne'].loc[Data['name'] == 'NAMA'] = 0
Data['deadOne'].loc[Data['name'] == 'GWI']  = 0

KhoisanRT   = PyInstEvo.ResolvedTree(Data.loc[Data['ruhlen_1'] == 'KHOISAN'], 'KTree1') #Create a resolved tree
numbranches = KhoisanRT.interiorbranches                                                  #Get no. of interior branches
bInit       = np.matrix(- 1 - np.linspace(0,10,num=numbranches)/numbranches)              #Make a conformable set of parameters
rInit       = np.zeros((1, len(KhoisanRT.words)))                                         #initial rate parameters
dparms      = np.sum(KhoisanRT.deathmat[:,0] == 0)                                        #Number of death parameters needed
dInit       = np.zeros((1, dparms)) + 1                                                   #Values for death parameters
eInit       = np.matrix(5)                                                                #Overall depth parameter
parmsInit   = np.hstack((bInit, rInit, dInit, eInit))   

KhoisanPT=PyInstEvo.ParameterizedTree(Data.loc[Data['ruhlen_1']=='KHOISAN'],'KTree1',parmsInit)

min = np.array(Depths['min'].loc[Depths['phylum'] == 'Khoisan'])
max = np.array(Depths['max'].loc[Depths['phylum'] == 'Khoisan'])
KhoisanPT.priordepth(min,max)
KhoisanPT.splitinfo(Splits[Splits['phylum'] == 'Khoisan'])

KhoisanPT.settimes()
KhoisanPT.showtree()

lnProbs = KhoisanPT.OriginLikelihood()

def DiscretePicker(lnProbs, toExclude=False):

    '''Pick out an index from a group of log-probabilities, when
       some entries are to be excluded as specified in toExclude'''
       
    lnProbs2u = np.copy(lnProbs)
    
    if toExclude != False:
        lnProbs2u[toExclude] = -10000
    Probs = np.exp(lnProbs2u)/np.sum(np.exp(lnProbs2u))
    rs = runningsum(Probs)
    pick = np.random.uniform(0, 1)
    
    return np.max(np.where(rs < pick)[0]) 


def runningsum(X):

    '''Running sum of the rows of an array'''

    rs=np.array([0])

    for i in range(1, len(X)):
        rs = np.vstack((rs, np.array(np.sum(X[0:i]))))
    
    return rs

  
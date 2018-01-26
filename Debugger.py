# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 10:23:21 2018

@author: mjbaker
"""

import PyInstEvo
# Other modules needed

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import mplleaflet
import re

from scipy.optimize import minimize

# Import all the routines that I have written.

#%matplotlib inline

# Reading in the data
# Read in the Pickle files

Data   = pd.read_pickle(os.getcwd() + '//IEData//MasterData.pkl')
Splits = pd.read_pickle(os.getcwd() + '//IEData//Splits.pkl')
Depths = pd.read_pickle(os.getcwd() + '//IEData//Depths.pkl')

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

KhoisanPT.RouteChooser()
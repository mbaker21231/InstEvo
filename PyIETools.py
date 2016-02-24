# Latest version! Working!

import pandas as pd
import numpy as np
import re
from scipy.stats import norm
from numba import jit
from numba import int32
from tqdm import tqdm
from scipy.stats import multivariate_normal 
from numpy.random import normal

def worddolg(word):
    
    '''This function takes a word and converts it into one of ten Dogolpolsky 
       classes. The classes are based on the first letter of each word falling 
       into one of ten groups listed below in the Dog list'''
       
    if len(re.sub(' ', '', word)) == 0:
        return np.nan
    
    W  = re.sub(' ', '', word)
    W  = re.sub('u', 'w', W[0])

    Dog=['pbf', 't8dT!', 'szSZ', 'cjgkqGX7Cx', 'm', 'Nn45',
         'lLr', 'wv', 'y', 'ieE3auoh']
    
    check = [word.find(W) for word in Dog]
    check = np.where(np.array(check) > -1)

    return(np.asscalar(check[0]))


def worddolgs(words):

    '''This function transforms a list or array of words into 
       Dogolpolsy classes'''
   
    return [worddolg(word) for word in words]


def charnum_to_dummies(W):
    
    '''Render numeric state variables into state vectors. I.e., State 3 becomes
       [0 0 1 0 0 0 0 0 0 0]. Just to mention something that should be 
       checked later - the type of the argument W is a little too specific
       to my application, I think.'''
       
    def rowmaker(W, j):
        '''Subfunction used in processing'''
        
        Row = np.asarray(pd.get_dummies(W[:,j]))
        for i in range(0, len(W[:,j])):
            if np.isnan(W[i,j]):
                Row[i,:] = np.ones((1, 10))
        return(Row)
    
    States = rowmaker(W, 0)
    
    for i in range(1, np.shape(W)[1]): 
        States = np.hstack((States, rowmaker(W, i)))
    
    return(States)   


def uniq(X):
    
    '''Returns unique values of a matrix. Saves us the trouble of
       converting everything into an array every time we want to find
       unique values. Really just saves me the trouble of typing a bunch
       of things every time'''
    
    return np.unique(np.asarray(X))


def reindex(Tree):

    '''First reindexing function I wrote. Note as general purpose as
       reindexer below, which can deal with generic grouped ids in any
       order. Still, I dont want to mess with it because it works on 
       newly imported trees. 
       
       Note that the renumbering starts at zero - more Pythonic!'''

    NewTree = np.matrix(np.zeros(np.shape(Tree)))
 
    for i in range(0, np.shape(Tree)[1]):
        vals = uniq(np.asarray(Tree[:,i]))
        counter = 0   
        for j in vals:
            matching = np.where(Tree[:,i] == j)
            NewTree[matching[0],i] = counter
            counter = counter+1
    
    LastRow = np.matrix(np.arange(0,np.shape(Tree)[0])).T
    
    NewTree = np.hstack((NewTree,LastRow))      
    
    return NewTree


def reindexer(Mat, start, end):
    
    '''More expansive reindexing function which is more elegant. 
       Takes a group of unique indexes and generates commensurate ones that run
       from 0 to N. This function works better on trees where the numbers arent
       necessarily in order to begin with.'''
    
    MatMod = np.matrix(np.copy(Mat))

    for i in range(start, end):
        u,ind = np.unique(np.array(Mat[:,i]), return_index = True)
        dictionary = dict(zip(u[np.argsort(ind)], np.arange(0, len(u))))
        Y = np.copy(Mat[:,i])
        for k, v in dictionary.items(): Y[MatMod[:,i] == k] = v
        MatMod[:,i]=np.copy(Y)
    
    return(MatMod)


def comptree(Tree):
    
    '''This function just eliminates redundant columns of a tree, in effect 
       reducing it to a smaller set of nested panels'''
    
    NT = Tree[:,0]
    
    for j in range(1, np.shape(Tree)[1]):
        if (np.array_equal(NT[:,-1], Tree[:,j]) == False) and \
            np.array_equal(uniq(NT[:,-1]), uniq(Tree[:,j])) == False:
            NT = np.hstack((NT, Tree[:,j]))
    
    return(NT)     

@jit(int32[:](int32[:], int32))
def ps(X, col):

    '''Panel setup ala Stata\Mata - creates rows with beginning and ending 
       indices of panel indicator variable. Note it does not renumber or 
       anything like that. Also, note that this returns an array, not a matrix.
       The return type is integer because in the future Python will have 
       problems with non-integer indices in matrices.'''

    vals = np.unique(np.asarray(X[:,col]))
    A = np.zeros(shape=(len(vals), 2))

    for i in range(0, len(vals)):
        spots = np.where(X[:,col] == vals[i])
        A[i,0],A[i,1] = spots[0][0], spots[0][-1]

    return(A.astype(int))

@jit
def psm(A, X, i):

    '''Panel submatrix returns the ith panel submatrix of X, i indexed as 
       it is in Python'''
    
    vals = np.arange(int(X[i,0]), int(X[i,1]) + 1)
    
    return(A[vals,])

	# A few more functions that basically simplify common Mata operations

@jit
def J(r, c, val):
    
    '''Empty matrix of arbitrary dimension, filled with val, which 
       can be whatever. This also mimics the Stata/Mata J() function.'''
    
    Mat = np.zeros((r, c))
    Mat[:,:]=val
    
    return(Mat)

@jit
def rows(Mat):
    
    '''Get the number of rows of a matrix'''
    
    return np.shape(Mat)[0]

@jit
def cols(Mat):
    
    '''Get the number of columns in a matrix'''
    
    return np.shape(Mat)[1]    

    
def resolvedTree(TreeBase):
    
    '''Randomly resolves trees written in matrix form as
       nested panels into bifurcating trees. Note that this has one 
       glaring weakness at this point. It peels off single observations
       at a time, rather than grouping and splitting. This needs to be 
       fixed.'''
    
    TreeMod = np.matrix(TreeBase, copy=True)
    
    while True:
        
        change = False
        
        for i in reversed(range(0, np.shape(TreeMod)[1] - 1)):
            m = ps(TreeMod, i)
        
            for j in range(0, np.shape(m)[0]):
                mAbove = psm(TreeMod[:,i+1], m, j)
                uniqvals = uniq(mAbove)
                countAbove = np.shape(uniqvals)[0]
            
                if countAbove > 2:
                    splitout = np.random.choice(uniqvals)
                    index = np.where(splitout == uniqvals)[0]
                    maxIndex = np.where(uniqvals == TreeMod[:,i+1])[0].max()
                    leftovers = np.delete(uniqvals, index)
                    newcol = np.copy(TreeMod[:,i])
                    leftoverIndex = np.array(np.where(TreeMod[:,i+1] == leftovers)[0]).astype(int)
                    newcol[leftoverIndex,0] = newcol[leftoverIndex,0] + 1
                    newcol[maxIndex+1:,0] = newcol[maxIndex+1:,0] + 1
                    TreeMod = np.hstack((TreeMod[:,0:i], TreeMod[:,i],
                                       newcol, TreeMod[:,i+1:]))
                    change = True
        
        TreeMod = TreeMod[np.lexsort(TreeMod.T[::-1])][0]
        TreeMod = reindexer(TreeMod, 0, np.shape(TreeMod)[1] - 1)
        
        if change == False:
            break

    return TreeMod    


@jit
def marknodes(Tree):
    
    '''Marks where the nodes are and numbers the subsequent groups following
       the node (actually, doesn't actually mark nodes but the idea is clear)
       For a bifurcating tree represented in matrix form'''
       
    Nodes = J(rows(Tree), cols(Tree), 1)
    Dum = J(rows(Tree), 1, 1)
    for j in range(0, cols(Tree)):
        m = ps(Tree,j)
        for i in range(0, rows(m)):
            Dump=psm(Dum, m, i)
            Nodes[m[i,0]:m[i,1]+1,j] = sum(Dump)
    return Nodes

@jit
def bpa(Tree,b):
    
    '''Finds locations where branching parameters should go, 
       using nodemarkers. It then puts the parameters where they need
       to be. The parameters are actual, not hyperparameters.'''
    
    M = marknodes(Tree)
    M.astype(int)

    BR = J(rows(Tree), cols(Tree), None)
    BR[:,-1] = 0
    B = np.hstack((J(1, rows(Tree), 0), b)) 

    count = rows(Tree)

    for j in reversed(range(1, cols(Tree))):
        m = ps(Tree, j)
        for i in range(0, rows(m)):
            if M[m[i,1],j] < M[m[i,1],j - 1] and M[m[i,1],j] > 1:
                BR[m[i,0]:m[i,1] + 1,j] = B[0,count]
                count += 1
    
    return BR    

@jit
def timeFractions(Tree, b, fill=False, DM=np.ones((1, 1)), r=1, a=0):
    
    '''Takes branch parameters and a death matrix and computes fractional times
       that make up the branches as we move along. Note that if we do not submit 
       a death matrix, every row should add to one. the parameters r and a help 
       in ensuring that branches do not get implausibly short, at least in
       theory, but at the moment dont seem to be behaving correctly. This
       needs to be checked.'''
       
    if cols(b) == 1: 
        Z = bpa(Tree, 1)      #Only one parameter makes for a funny tree!
    else: 
        Z = bpa(Tree, b[0,0:-1]) 

    Z[:,0] = b[0,-1]

    V = J(rows(Z), cols(Z), np.nan)
    F = J(rows(Z), cols(Z), np.nan)

    F[:,0:1] = (1 + (1 - r)*(np.exp(Z[:,0:1]) + a))/(1 + np.exp(Z[:,0:1] + a))
    V[:,0:1] = (np.exp(Z[:,0:1]) + a)*r/(1 + np.exp(Z[:,0:1]) + a)
    
    # Code block to progressively keep track of time left and time elapsed. 
    i = 1
    while True:
        V[:,i:i+1] = F[:,i-1:i]*(np.exp(Z[:,i:i+1]) + a)*r/(1 + np.exp(Z[:,i:i+1]) + a)
        F[:,i:i+1] = F[:,i-1:i]*(1 + (1 - r)*(np.exp(Z[:,i:i+1])+a))/(1 + np.exp(Z[:,i:i+1] + a))
        ind=np.where((F[:,i:i+1] <= 1) == False)[0]   # This line and the next one edit nans to previous value. 
        F[ind,i:i+1] = F[ind,i-1:i]
        i += 1
        if i > cols(V):
            break  
        
    V[:,-1] = F[:,-2]
    VV=J(rows(V), cols(V), np.nan)
    VV2=J(rows(V), cols(V), np.nan)
    VV[:,-1] = V[:,-1]
    VV2[:,-1] = V[:,-1]
    
    for i in reversed(range(0, cols(Tree))):
        m = ps(Tree, i)                     
        for j in range(0, rows(m)):
            VV[m[j,0]:m[j,1]+1,i] = V[m[j,1],i]   
            VV2[m[j,1],i] = V[m[j,1],i]           
    
    # Roll back the last branches if we have an expiration date for anyone...
    # Note that if fill is true, ONLY the filled value is returned. 
    # Otherwise, the first thing returns is filled, and the second is not!
    
    if np.any(DM[:,0] == 0):            
        ind = np.where(DM[:,0] == 0)[0]    
        VV[ind,-1] = 1/(1 + np.exp(DM[ind[0],-1]))
        VV2[ind,-1] = 1/(1 + np.exp(DM[ind[0],-1]))

    if fill == True:
        return VV     
    else:
        return VV, VV2

def branchcount(Tree):
    
    '''Count the number of branches on a tree'''
    
    M = marknodes(Tree)
    count = rows(Tree)
    
    for j in reversed(range(1, cols(Tree))):
        m = ps(Tree, j)
        for i in range(0,rows(m)):
            if (M[m[i,1],j] < M[m[i,1],j-1]) and (M[m[i,1],j] > 1): count += 1 

    return count


def branchpos(Tree):

    '''Catalog the positions of branches in a tree matrix in a nX2 matrix.'''
    '''Should operate under assumption that the unfilled tree is to be 
       returned. '''

    TFS = timeFractions(Tree, np.matrix(np.zeros((1, branchcount(Tree)))))
    TreeBP = TFS[1]
    
    BPs=J(0, 2, np.nan)
    
    for j in reversed(range(0, cols(Tree))):
        for i in range(0, rows(Tree)):
            if TreeBP[i,j] >= 0:
                BPs=np.vstack((BPs, (i, j)))
    
    return BPs    

    
def BuiltBranchNodes(Treeo,Tree,ibs):

    '''Arguments here are the original tree in matrix form, the resolved tree, 
       and the row,column positions of the branches of the tree. The output 
       will be the set of branches in the tree that were part of the resolution
       of the tree, not the original, base Tree'''

    builtCols=J(0,1,np.nan)

    for i in range(0, cols(Tree) - 1):
        foundit = 0
        j = 0
        while True:
            if np.all(Tree[:,i] == Treeo[:,j]):
                foundit = 1
            else:
                j += 1
            if (foundit == 1) or (j > cols(Treeo) - 1):
                break
        if foundit == 0:
            builtCols=np.vstack((builtCols, i))

    # Second part - mark the nodes that are built according to the above code
    marker=np.zeros((rows(ibs), 1))

    for i in builtCols:
        pos = np.where(ibs[:,1] == i)
        marker[pos] = marker[pos]+1

    posOfBuilt = np.where(marker != 0)[0]
    builtBranches = ibs[posOfBuilt,:]       
        
    # Third part - double check the branches on the list just to be sure!
    builtBranchesFinal = J(0, 2, np.nan)
    for i in range(0, rows(builtBranches)):
        id = Tree[builtBranches[i,0],builtBranches[i,1]]    
        p1 = len(np.where(Tree[:,builtBranches[i,1]] == id)[0].T)
        flag = 0
        for j in range(1, cols(Treeo)):
            ido = Treeo[builtBranches[i,0],j]
            p2 = len(np.where(Treeo[:,j] == ido)[0].T)
            if p1 == p2:
                flag = 1
        if flag == 0:
            builtBranchesFinal = np.vstack((builtBranchesFinal, builtBranches[i,:]))    
    
    return(builtBranchesFinal)


def jcallQ(rates, branch):

    '''Function to make rates and branches into one large matrix of stacked 
       Jukes-Cantor transition  matrices.'''

    # Diagonals of vertically-stacked transition matrices

    eR = np.exp(-rates.T*branch)
    dR = np.kron((1/10 + 9/10*eR), np.eye(10))
    np.shape(dR)
    
    #Off diagonals

    odR1 = np.matrix(1/10*(1 - np.kron(eR,np.ones((10, 10)))))
    odR2 = (J(rows(eR)*10, 10, 1) - np.kron(J(rows(eR), 1, 1),np.eye(10)))
    odR = np.multiply(odR1, odR2)
    
    return dR + odR  


def vecvecmult(A, B):
    
    '''Function to multiply sequences of vectors - represented as rows of a 
       matrix - in parallel'''
       
    return np.matrix(np.sum(np.multiply(A, B), axis=1))


def matvecmult(A, B):

    '''Function to multiply sequences of matrices - represented as rows of a 
       meta-matrix with sequences of vectors. Matrices must all be of the same 
       dimension and square'''

    C = J(rows(A), 0, 0)
    dim = np.sqrt(cols(A))
    dim = dim.astype(int)
    
    for i in range(0, cols(A), dim):
        C = np.hstack((C, vecvecmult(A[:,i:i+dim], B)))
    
    return(C)  


def JCallFast(Tree,branchPos,rates,branchParms,States,depth,DM):

    '''Computes the likelihood of a tree given rates, branch parameters,
       depth, etc.'''
    TFS = timeFractions(Tree, branchParms, False, DM)  #Do we need to call this?
    TreeBL = TFS[1]
    
    S = np.copy(States)
    dim = cols(S)/10
    
    for i in range(0, rows(TreeBL)):
        bl = TreeBL[i,-1]
        Q = jcallQ(rates,bl)
        Qhat = np.reshape(Q, (dim, 10*10))   
        Shat = np.reshape(S[i,:],(dim, 10))
        Snew = matvecmult(Qhat,Shat)
        S[i,:] = np.reshape(Snew,(1,dim*10)) 

    S=np.log(S)    

    # Actual Tree Loop:
    for i in range(rows(Tree), rows(branchPos)):
        r,c = branchPos[i,:]
        id = Tree[r,c]
        p = (np.matrix(S[np.where(Tree[:,c]==id)[0],:]))
        bl = TreeBL[r,c]         # Retrieve branch to use
        Sn = np.nansum(p,axis=0) 
        Qn = jcallQ(rates,bl)
        Sn = np.reshape(Sn, (dim,cols(Sn)/dim))
        Qn = np.reshape(Qn, (rows(Qn)/10, cols(Qn)*10)) 
        ind = np.where(Tree[:,c] == id)[0]      
        S[ind,:] = np.nan
        S[r,:] = np.reshape(Sn, (1, rows(Sn)*cols(Sn)))     
        
    # Collapsation and report:
    S = np.matrix(S[-1,:])
    S = np.reshape(S, (dim, cols(S)/dim))    
    rowmax = np.amax(S,axis=1)
    S = rowmax + np.log(np.dot(np.exp(S - rowmax), np.ones((10, 1))))
    
    return(sum(S)) 
  
    
def FractionCommon(Tree,TreeBL,Splits,Names):
    
    '''Takes in the Tree - with the last column reflecting the ordering of the 
       Names.'''
    
    '''Computes the fraction of the Tree time that each of the groups in Splits 
       data were _the_same_ group'''
       
    Order = np.array(Tree[:,-1]).flatten().astype(int)
    SplitsArray = np.array(Splits)

    Names1 = np.array(SplitsArray[:,1]).flatten()
    Names2 = np.array(SplitsArray[:,2]).flatten()

    NameInd1 = np.where(Names[Order] == Names1)[0]
    NameInd2 = np.where(Names[Order] == Names2)[0]
    SplitIndices = np.vstack((NameInd1, NameInd2)).T
    JointTimes = J(rows(SplitIndices), 1, 0)

    i=0
    for Split in SplitIndices:
        ToSum = np.where((Tree[Split[0],:] == Tree[Split[1],:]))[1]
        JointTimes[i,0] = np.nansum(TreeBL[Split[1],ToSum])
        i += 1
    
    return np.hstack((SplitIndices,JointTimes))


def SplitLikelihood(Tree,TreeBL,Splits,Names,depth):

    ''' Uses FractionCommon and returns the log likelihood of the splits, 
        given the prior mean and standard deviation in the file.'''

    TimeSplit = np.matrix(1000*depth*(1 - FractionCommon(Tree, TreeBL,
                                                         Splits, Names)[:,2])).T
    SplitMat = np.matrix(Splits)
    Vals = (TimeSplit - SplitMat[:,3])/SplitMat[:,4]
    Vals = Vals.astype(float)
    
    return np.sum(norm.logpdf(Vals))     
    

def gcircledist(lat, lon):

    '''Compute great circle distance between points in miles.'''    
    
    miles = 6372.8
    degs = np.pi/180
    lonConverted = np.copy(lon)
    lonConverted[lonConverted < 0] = 360 + lonConverted[lonConverted < 0]

    rlat = lat*degs
    rlon = lonConverted*degs
    D=J(len(lat), len(lat), 0)
    for i in range(0, len(lat)):
        for j in range(0, len(lon)):
            if i == j:
                D[i, j] = 0
            else:
                D[i,j] = miles*np.arccos(np.sin(rlat[i])*np.sin(rlat[j]) + 
                    np.cos(rlat[i])*np.cos(rlat[j])*np.cos(rlon[i] - rlon[j]))
    
    return (D + D.T)/2    


def LnLikeOrigins(Tr,D,Tree,bp):
    
    '''A function that takes in a time-filled in Tree (Tr), a pairwise distance
       matrix D, the actual Tree in Matrix form (Tree), and a list of the 
       row,column branch positions of the Tree. It puts
       out the likelihood of each row in Tree being the point 
       of origin for the tree. Do we ever use this function? '''
    
    '''Some things to watch out for - the np.where function can behave in 
       funny ways! Also, it is critical that maxFinder and toAdd have the 
       right shape - column vectors, effectively.'''
    
    '''In trouble shooting, a frequent culprit causing problems were 
       programming snafus that left TM empty (equal to zero) while at the 
       same time attempting to take a log.'''
    
    '''Another thing to watch out for - branches that are so small they 
       violate Poisson rareness assumption!'''
    
    # Initialize placeholders:
    
    TM = Tr[:,-1]
    LL = np.zeros((rows(TM), 1))
    IndCount = np.zeros((rows(TM), 1))
    Live = np.zeros((rows(TM), 1))        

    # Get rid of zeros in the Distance matrix:

    D[D <= 0] = 10 
    lnD = np.matrix(-np.log(D))
    np.fill_diagonal(lnD, 0) 

    # The main recursive loop backwards through the branches:

    for i in range(rows(Tree), rows(bp)):     
        id = Tree[bp[i,0], bp[i,1]]
        tu = np.where(Tree[:,bp[i,1]] == id)[0]
        Bhat = Tr[tu,bp[i,1]:]
        TreeHat = Tree[tu,bp[i,1]:]
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
            posi =  np.where(nums == m)[0]
            posni = np.where(nums != m)[0]
            toAdd = np.matrix(np.nansum(Bhat[posni,:], axis = 1)).T
            for q in posi:
                maxFinder = DHat[q,posni].T
                for n in range(0, len(posni)):
                    if IndCountHat[posni[n]] >= 1:
                        maxFinder[n] = maxFinder[n] + IndCountHat[posni[n]]*(np.log(IndCountHat[posni[n]]) - np.log(toAdd[n])) + LLHat[posni[n]]
                max = np.argmax(maxFinder)         
                TMHat[q] = toAdd[max]
                IMHat[q] = 1 + IndCountHat[posni[max]]
                LLHat[q] = DHat[q,posni[max]]
            TMHat[posi] = toAdd[max]
            IMHat[posi] = IndCountHat[posni[max]] + 1
        
        for p in range(0, rows(tu)):
            if Live[tu[p]] == 1:
                LL[tu[p]] = LL[tu[p]] + IndCount[tu[p]]*(np.log(IndCount[tu[p]]) - np.log(TM[tu[p]]))+LLHat[p]
                TM[tu[p]] = np.copy(TMHat[p])
                IndCount[tu[p]] = np.copy(IMHat[p])
            else:
                LL[tu[p]] = LLHat[p]
                Live[tu[p]] = 1
                TM[tu[p]] = TMHat[p]
                IndCount[tu[p]] = IMHat[p]   
    
    return LL


def runningsum(X):

    '''Running sum of the rows of an array'''

    rs=np.array([0])

    for i in range(1, len(X)):
        rs = np.vstack((rs, np.array(np.sum(X[0:i,:]))))
    
    return rs

def DiscretePicker(lnProbs, toExclude=False):

    '''Pick out an index from a group of log-probabilities, when
       some entries are to be excluded as specified in toExclude'''
       
    lnProbs2u = np.copy(lnProbs)
    
    if toExclude != False:
        lnProbs2u[toExclude] = -10000
    Probs = np.exp(lnProbs2u)/np.sum(np.exp(lnProbs2u))
    rs = runningsum(Probs)
    pick = np.random.uniform(0, 1)
    
    return max(np.where(rs < pick)[0])    


def OriginLikelihood_args(Tree,bp,Dhat,TFR):
    
    '''A function that takes in a time-filled in Tree (Tr), a pairwise distance 
       matrix D, the actual Tree in Matrix form (Tree), and a list of the row,
       column branch positions of the Tree. It puts
       out the likelihood of each row in Tree being the point of origin for the
       tree. '''

    '''Some things to watch out for - the np.where function can behave in funny
       ways! Also, it is critical that maxFinder and toAdd have the right shape
       - column vectors, effectively.'''

    '''In trouble shooting, a frequent culprit causing problems were 
       programming snafus that left TM empty (equal to zero) while at the same 
       time attempting to take a log.'''
    
    '''Another thing to watch out for - branches that are so small they violate
       Poisson rareness assumption!'''

    # Initialize placeholders:
    
    TM = np.copy(Tree[:,-1])   #Just to shorten the code a bit
    LL = np.zeros((rows(TM),1))
    IndCount = np.zeros((rows(TM),1))
    Live = np.zeros((rows(TM),1))        

    # Get rid of zeros in the Distance matrix:
    
    D = np.copy(Dhat)
    D [D <= 0] = 10                                      # Here again I love Python!
    lnD = np.matrix(-np.log(D))
    np.fill_diagonal(lnD, 0) 

    # The main recursive loop backwards through the branches:
    for i in range(rows(Tree), rows(bp)):     
        id = Tree[bp[i,0], bp[i,1]]
        tu = np.where(Tree[:,bp[i,1]]==id)[0]
        Bhat = np.copy(TFR[tu,bp[i,1]:])
        TreeHat = np.copy(Tree[tu,bp[i,1]:])      # Better safe than sorry!
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
                for n in range(0,len(posni)):
                    if IndCountHat[posni[n]]>=1:
                        maxFinder[n] = maxFinder[n] + IndCountHat[posni[n]]*(np.log(IndCountHat[posni[n]]) - np.log(toAdd[n])) + LLHat[posni[n]]
                max = np.argmax(maxFinder)           # No need to do more - defaults to first entry if more than one
                TMHat[q] = toAdd[max]
                IMHat[q] = 1 + IndCountHat[posni[max]]
                LLHat[q] = DHat[q,posni[max]]
            TMHat[posi] = toAdd[max]
            IMHat[posi] = IndCountHat[posni[max]]+1
     
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
    
    return LL


    
    
def RouteChooser(ParameterizedTreeClass):
    
    '''Recursively applies OriginLikelihood to pick the most
       likely point occupied at any time.'''
       
    Tr = 100*ParameterizedTreeClass.filledtimeFractions*ParameterizedTreeClass.depth
    Tree = ParameterizedTreeClass.resolvedtree
    bp = ParameterizedTreeClass.branchpositions
    D = ParameterizedTreeClass.D

    branchRoute = J(0, 2, np.nan) 
    TreeHat = np.hstack((Tree[:,0:-1], np.matrix(np.arange(0,rows(Tree))).T))        #Renumbered tree

    # Compute origin probability....
    lnProbs = ParameterizedTreeClass.OriginLikelihood()
    
    init = DiscretePicker(lnProbs)
    Path = J(rows(Tree), cols(Tree), np.nan)
    Path[:,0] = init   
    
    for i in reversed(range(rows(Tree), rows(bp) - 1)): 
 
        ind = bp[i,:]          # Retrieve an index
        z = ind[1] - 1
        found = 0
    
        groupId = TreeHat[ind[0],ind[1]]
        colns = np.where(TreeHat[:,ind[1]] == groupId)[0]
        while True:
            if Path[max(colns),z] >= 0:               # Really just looking for np.nans!
                found = 1
                origin = Path[max(colns),z].astype(int)
            else: z -= 1
            if found == 1:
                break
            
        if len(np.where(colns != origin)[0]) == 1:
            destination = colns[np.where(colns != origin)[0]]
        elif any(colns == origin):
            destination = origin
        else:
            subTree = TreeHat[colns,ind[1]:]
            subTr = Tr[colns,ind[1]:]
            subBr = branchpos(subTree)
            subD = D[colns,:][:,colns]
            contProb = OriginLikelihood(subTree, subBr, subD, subTr)
            toGet = np.matrix(-np.log(D))[origin,:][:,colns].T
            contProb = contProb + toGet
            destination = colns[DiscretePicker(contProb, np.where(colns == origin)[0])]
        branchRoute = np.vstack(([origin,destination], branchRoute))
        colns = colns[np.where(colns!=origin)[0]]
        Path[colns,ind[1]] = destination    
        
        # Cleaning up a few things from the route.. 
 
    for i in range(0,rows(Path)):
        if (np.any(branchRoute[:,1]==i)):
            pass
        else:
            s = cols(Path)- 2    #Not the last column, but the one before that!
            originFound = 0
            while True:
                if Path[i,s] >= 0:               # Really just looking for np.nans!
                    originFound = 1
                else: s -= 1
                if originFound == 1:
                    if Path[i,s] != i:
                        branchRoute=np.vstack((branchRoute, (Path[i,s] ,i)))
                    break
    
    return Path, branchRoute    
	
@jit
def OriginLikelihood(PhyTree):
        
    '''A function that takes in a time-filled in Tree (Tr), a pairwise distance matrix D, the actual
       Tree in Matrix form (Tree), and a list of the row,column branch positions of the Tree. It puts
       out the likelihood of each row in Tree being the point of origin for the tree. '''
    
    '''Some things to watch out for - the np.where function can behave in funny ways! Also, it is critical
       that maxFinder and toAdd have the right shape - column vectors, effectively.'''
    
    '''In trouble shooting, a frequent culprit causing problems were programming snafus that
       left TM empty (equal to zero) while at the same time attempting to take a log.'''
    
    '''Another thing to watch out for - branches that are so small they violate Poisson rareness assumption!'''
    
    # Initialize placeholders:
    
    TM = np.copy(PhyTree.resolvedtree[:,-1])   #Just to shorten the code a bit
    bp = PhyTree.branchpositions      #Likewise
    LL = np.zeros((rows(TM), 1))
    IndCount = np.zeros((rows(TM), 1))
    Live = np.zeros((rows(TM), 1))        

    # Get rid of zeros in the Distance matrix:
    D = np.copy(PhyTree.D)
    D[D <= 0] = 10   
    lnD = np.matrix(-np.log(D))
    np.fill_diagonal(lnD,0) 

    # The main recursive loop backwards through the branches:
    for i in range(rows(PhyTree.resolvedtree), rows(bp)):     
        id = PhyTree.resolvedtree[bp[i,0], bp[i,1]]
        tu = np.where(PhyTree.resolvedtree[:,bp[i,1]] == id)[0]
        Bhat = PhyTree.depth*1000*np.copy(PhyTree.filledtimeFractions[tu,bp[i,1]:])
        TreeHat = np.copy(PhyTree.resolvedtree[tu,bp[i,1]:])
        IndCountHat = np.copy(IndCount[tu])
        LLHat = np.copy(LL[tu])
    
        z = 0
        while True:
            ids = np.unique(np.asarray(TreeHat[:,z]))             # Replaced again...
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
    
    return LL 	

@jit
def settimes(PhyTree):
        
    ''' Set times along the tree in both filled and unfilled form. We probably want a method
        to display and check that the time fractions are working okay.
    '''
    mind = PhyTree.depthmin
    maxd = PhyTree.depthmax
    midd = PhyTree.eparms
    PhyTree.depth = np.exp(midd)/(1 + np.exp(midd))*mind + 1/(1 + np.exp(midd))*maxd
    TFS = timeFractions(PhyTree.resolvedtree, PhyTree.bparms,False, PhyTree.deathmat)     
    PhyTree.filledtimeFractions = TFS[0]
    PhyTree.timeFractions = TFS[1]

def myMcMcSampler_mwg(lnlike, xinit, Vinit, draws, burn, damper, aopt, *args):
 
    nb     = np.shape(xinit)[1]
    xold   = xinit
    lam    = 2.38**2/nb*np.ones((1,nb))
    old    = lnlike(xold, *args)
    val    = []
    Accept = np.zeros((1,nb))
    alpha  = np.zeros((1,nb))
    xs     = []
    mu     = np.array(xold).flatten()

    for i in tqdm(range(0,draws), leave=True):
        accept = np.zeros((1,nb))
        for j in range(0,nb):
            xpro = np.copy(xold)
            xpro[0,j] = xold[0,j]+multivariate_normal.rvs(0)*np.sqrt(Vinit[j,j])*lam[0,j]
            xpro = np.reshape(np.matrix(xpro),(1,-1))
            pro = lnlike(xpro, *args)                        
            if np.isnan(xpro[0,j]):
                alpha[0,j] = 0
            elif pro > old:
                alpha[0,j] = 1
            else:
                alpha[0,j] = np.exp(pro - old)
            
            if np.random.uniform(0,1) < alpha[0,j]:
                old    = pro
                xold   = xpro
                accept[0,j] = 1
        
        lam = lam*np.exp(1/(i+1)**damper*(alpha-aopt)) 
        xs.append(xold)
        val.append(old)
        Accept = np.vstack((Accept,accept))
    
    return xs, val, Accept
	
def myMcMcSampler_global(lnlike, xinit, Vinit, draws, burn, damper, aopt, *args):
 
    nb   = np.shape(xinit)[1]
    xold = xinit
    lam  = 2.38**2/nb
    old  = lnlike(xold, *args)
    val  = []
    Accept = []
    xs     = []
    mu     = xold

    V      = np.matrix(Vinit)
    Vold   = np.eye(nb)

    for i in tqdm(range(0,draws), leave=True):
        accept = 0
        try:
            xpro = myNormalDrawer(mu,V)    # Not the cleanest way to do this
        except:
            xpro = myNormalDrawer(mu,Vold)
            V    = Vold
        xpro = np.reshape(np.matrix(xpro),(1,-1))
        pro = lnlike(xpro, *args)                        

        if np.any(np.isnan(xpro)):
            print('Having problems with positive definiteness...')
            alpha = 0
        elif pro > old:
            alpha = 1
        else:
            alpha = np.exp(pro - old)
            
        if np.random.uniform(0,1) < alpha:
            old    = pro
            xold   = xpro
            accept = 1
  
        lam = lam*np.exp(1/(i+1)**damper*(alpha-aopt)) 
        xs.append(xold)
        val.append(old)
        Accept.append(accept)
        mu = mu+1/(i+1+1)**damper*(xold-mu)
        Vold=np.copy(V)
        V=V+1/(i+1+1)**damper*((xold-mu).T*(xold-mu)-V)
        V=(V+V.T)/2

    return xs, val, Accept	
	
@jit
def myNormalDrawer(mean, Var):
    nb = np.shape(mean)[1]
    evec = np.matrix(normal(size = nb))
    cholVar=np.linalg.cholesky(Var)
    return mean + np.dot(evec, cholVar) 
	
def make_pos_def(X, eps = .001):
    newX=np.matrix(np.copy(X))
    counter=0
    while True:
        for i in range(np.shape(X)[0]):
            newX[i,i]=newX[i,i]+eps
        if is_pos_def(newX):
            break
        counter+=1
    print("rendered positive definite: ", counter, " iterations")
    return(newX)	

def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) > 0)

def mlfun(x, Obj):
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
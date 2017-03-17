# Collecting all the information for class building

import numpy as np
import matplotlib.pyplot as plt
import mplleaflet          

from PyIETools import *



class PhyloTree:
    
    '''A class for holding basic information and coincident data from a tree. 
       Note that it is assumed that Data is named accordingly (i.e., the 
       Swadesh words, names of groups, numeric nested panels, etc.) One can see
       below that words are called words* in the data, and the nested (possibly
       not completely resolved) tree comes in the form of variables named TR*. 
       Moreover, the class assumes that we have an instance. Also, location
       information is assumed to be labelled according to lat and lon in the 
       Data.'''
    
    # Effectively, global variables that are used by the methods below.
    
    words = ['word1', 'word2', 'word3', 'word11', 'word12', 'word18', 'word19',
           'word21', 'word22', 'word23', 'word25', 'word28', 'word30', 
           'word31', 'word34', 'word39', 'word40', 'word41', 'word43', 
           'word44', 'word47', 'word48', 'word51', 'word53', 'word54', 
           'word57', 'word58', 'word61', 'word66', 'word72', 'word74',
           'word75', 'word77', 'word82', 'word85', 'word86', 'word92', 
           'word95', 'word96', 'word100']
        
    namelist = [ 'TR' + str(x) for x in range(1, 17) ]  

    dim1 = len(words)
    DogList = ['p', 't', 's', 'c', 'm', 'N', 'l', 'w', 'y', 'i']

    # Creates a list of all the names of states variables. 
    
    statenames=[]
    for i in range(0, dim1):
        for j in range(0, 10):
            statenames.append(words[i] + str(i) + '_' + DogList[j])             
    
    # We don't want to carry around counters as part of the class!
    del i
    del j
    
    def __init__(self, Data, Title):

        '''Initiates a basic tree and its information. The Data should be the 
           data set, which the class then parses and organizes in the correct 
           format - i.e., picks out the words to use, translates these into
           Dolgopolsky classes, and then into states, and organizes the nested
           panel structure so that it is numbered nicely.'''

        # First, just save all the data associated with the Tree
        self.title = Title
        self.Data = Data
        
        # Squash the Basic tree to its essence, renumber it. 
        
        self.BaseTree = comptree(reindex(np.matrix(Data[self.namelist])))
        
        # Add a name vector to the base tree
        
        self.name = np.matrix(Data['name']).T # Transpose to column vector

        # Pull out the names part of the data and make form is right. 

        self.states = np.matrix(Data[self.statenames])  
        
        '''Words interpreted as horizontally stacked matrices of dummies 
           10 categories in total
           '''
        # Compute Great-Circle Distances between groups in the Tree:
 
        (self.lat, self.lon) = np.array(Data['lat']), np.array(Data['lon'])       
        
        '''Latitudes and Longitudes of each name on the Tree, which can then be used to create
           great circle distance between points. Our distance matrix is called D'''
        
        self.D=gcircledist(self.lat, self.lon)
        
        # Set the deathmatrix with parameter spaces (marked with 10 in the second column)
        # and whether or not the language is dead, marked with a zero in the first column
        # Note: we are going to have to think about how we want to calculate deaths, as we
        # have yet to put this in our calculations. 
        
        self.deathmat = np.hstack((np.matrix(Data.deadOne).T, 
                                   J(rows(Data), 1, 10)))
        
        self.deathdata = np.hstack((np.matrix(Data.ex_date).T,
                                  np.matrix(Data.ex_date_sd).T))
        
    def display(self):
        
        '''Print out the structure of the Tree...with names affixed on the end'''

        print(np.hstack((self.BaseTree.astype(int),self.name)))
        
    def matplot(self):
        
        '''The not-so-useful matrix plot...but at least fun to have as a method.'''
        
        valsInRow=J(1, cols(self.BaseTree), 0)
        for i in range(0, cols(self.BaseTree)):
                valsInRow[0,i] = len(uniq(self.BaseTree[:,i]))
        
        f = plt.figure(1)
        sp = f.add_subplot(111)
        mattoplot = np.hstack((self.BaseTree[:,:-1],
                               np.matrix(np.arange(0,rows(self.BaseTree))).T))
        sp.matshow(mattoplot/valsInRow, cmap='viridis')

    # Some methods to assign prior information:
    
    def priordepth(self, depthmin, depthmax):
        
        '''Define upper and lower bound for overall time depth of the tree. 
           Assumption will eventually be that the true depth is uniformly 
           distributed between these two points'''
        
        self.depthmin = depthmin
        self.depthmax = depthmax
        
    def splitinfo(self, Splits):
        
        '''Affix a set of split times to the tree. note that these come in in 
           the form of data from a dataframe, which we meld into an array. 
           The array has the format name1,name2, years_apart,stdev_years_apart.
           We drop the first column because it contains the (redundant)
           name of the Phylogeny.'''
        
        self.splittimes = np.array(Splits)[:,1:]



class ResolvedTree(PhyloTree):
    
    '''This class extends the previous class to include a fully resolved tree. Once the tree has been resolved, one
       can also add a count of the number of branches on the tree, and one can also add in dimension info, which 
       tells us where the parameter vector should be broken apart.'''
    
    def __init__(self, Data, Title):
        
        PhyloTree.__init__(self,Data,Title)

        # Adding in the key points of differentiation for the resolved tree:
        
        self.resolvedtree = resolvedTree(self.BaseTree)
        self.branchpositions = branchpos(self.resolvedtree).astype(int)
        self.numberbranches = branchcount(self.resolvedtree)
        self.interiorbranches = self.numberbranches + 1 - len(self.name)
        self.builtcols = BuiltBranchNodes(self.BaseTree, self.resolvedtree,
                                          self.branchpositions)
        
        # Add up the dimension information for the resolved tree:
        a = self.numberbranches + 1 - rows(self.resolvedtree)
        b = len(self.words) + a
        c = b + len(self.deathmat[:,0][self.deathmat[:,0] == 0].T)    
        self.dimInfo=(a, b, c, -1)
        
        # We need to be careful that we are getting the right order of everything...
        
        self.order=np.array(self.resolvedtree[:,-1].astype(int)).flatten()
        
        # Reorder the death matrix so it is conforming, along with the deathdata
        
        self.deathmat = self.deathmat[self.order,:]
        self.deathdata = self.deathdata[self.order,:]
        self.states = self.states[self.order,:]
        self.name = self.name[self.order,:]
        
    def matplot(self):
        
        '''The not-so-useful matrix plot...applied to the subclass, not the 
           original with a little method overloading. There is probably a much
           cleaner way of doing this!'''
        
        valsInRow = J(1, cols(self.resolvedtree), 0)
        for i in range(0, cols(self.resolvedtree)):
                valsInRow[0,i] = len(uniq(self.resolvedtree[:,i]))
        
        mattoplot = np.hstack((self.resolvedtree[:,:-1],
                               np.matrix(np.arange(0, rows(self.resolvedtree))).T))       
        f = plt.figure(1)
        sp = f.add_subplot(111)
        sp.matshow(mattoplot/valsInRow, cmap='viridis')


class ParameterizedTree(ResolvedTree):
    
    '''Here we add a parameter vector to the Tree. The dimension information 
       and all that have previously been set up.  '''
    
    def __init__(self,Data,Title,Parameters):
        ResolvedTree.__init__(self, Data, Title)
        self.parameters = Parameters
        self.unpack()

    def unpack(self):
        self.bparms = self.parameters[0,0:self.dimInfo[0]]
        self.rparms = self.parameters[0,self.dimInfo[0]:self.dimInfo[1]]
        self.dparms = self.parameters[0,self.dimInfo[1]:self.dimInfo[2]]
        self.eparms = self.parameters[0,self.dimInfo[3]]
        self.rates = np.exp(self.rparms)
    
    def settimes(self):
        
        ''' Set times along the tree in both filled and unfilled form. We probably want a method
            to display and check that the time fractions are working okay.
        '''
        mind = self.depthmin
        maxd = self.depthmax

        midd = self.eparms
        self.depth = np.exp(midd)/(1 + np.exp(midd))*mind + 1/(1 + np.exp(midd))*maxd
        TFS = timeFractions(self.resolvedtree, self.bparms,False, self.deathmat)     
        self.filledtimeFractions = TFS[0]
        self.timeFractions = TFS[1]
    
    def showparmcount(self):
        a = self.numberbranches+1-len(self.name)
        print("Number of branch parameters: ", a)
        b = len(self.words)
        print("Number of rate parameters:   ", b)
        DM = self.deathmat
        c = cols(DM[:,0][DM[:,0] == 0])
        print("Number of dead branches:     ", c)
        print("Overall depth parameters:    ", 1)
        print("Total parameters:            ", a + b + c + 1)
        print('')
        
    def showparameters(self):
        print('Branch parameters:    ', self.bparms)
        print('Rate parameters (ln): ', self.rparms)
        print('Death parameters:     ', self.dparms)
        print('Overall depth parms:  ', self.eparms)
        print('')
        
    '''We are in fact going to build a few methods right into the class...
       even though we have them as stand-alone functions. This will allow us to
       save time by computing things only once. These methods will allow us to 
       compute our entire likelihood without doing things
       over and over again (i.e., placing branches, computing time fractions, etc.)'''
    
    def jukescantorlikelihood(self):

        '''Computes Jukes Cantor likelihood of the tree, given the parameters 
           and the resolution.'''
        '''We make the sub-function a part of the actual function here. '''
        
        def jcallQ1(rates,branch):
            '''Function to make rates and branches into one large matrix of 
               stacked Jukes-Cantor transition
               matrices'''
            # Diagonals of vertically-stacked transition matrices
            eR = np.exp(-rates.T*branch)
            dR = np.kron((1/10 + 9/10*eR), np.eye(10))
            #Off diagonals
            odR1 = np.matrix(1/10*(1-np.kron(eR,np.ones((10,10)))))
            odR2 = (J(rows(eR)*10, 10, 1) - np.kron(J(rows(eR),1,1), np.eye(10)))
            odR = np.multiply(odR1, odR2)
            return dR+odR     
   
        TreeBL = self.timeFractions
    
        S = np.copy(self.states)
        dim = cols(S)/10
    
        # Initialization of the Tree
        for i in range(0, rows(TreeBL)):
            bl = TreeBL[i,-1]
            Q = jcallQ1(self.rates, bl)
            Qhat = np.reshape(Q,(dim, 10*10))   # Essentially following our old method...maybe not the most efficient...
            Shat = np.reshape(S[i,:],(dim, 10))
            Snew = matvecmult(Qhat, Shat)
            S[i,:] = np.reshape(Snew, (1, dim*10)) 

        S = np.log(S)    

        # Actual Tree Loop:
        for i in range(rows(self.BaseTree), rows(self.branchpositions)):
            r, c = self.branchpositions[i,:]
            id = self.resolvedtree[r,c]
            p = (np.matrix(S[np.where(self.resolvedtree[:,c] == id)[0],:]))
            bl = TreeBL[r,c]                # Retrieve the branch we are going to use
            Sn = np.nansum(p,axis=0) # Sum of logs of columns corresponds to multiplication
            Qn = jcallQ1(self.rates,bl)
            Sn = np.reshape(Sn,(dim,cols(Sn)/dim))
            Qn = np.reshape(Qn,(rows(Qn)/10,cols(Qn)*10)) 
            ind = np.where(self.resolvedtree[:,c]==id)[0]      # Be careful wiht this function! The tailing zero is needed I think!
            S[ind,:] = np.nan
            S[r,:] = np.reshape(Sn,(1,rows(Sn)*cols(Sn)))     # Needed or not?
        
    # Collapsation and report:
        S = np.matrix(S[-1,:])
        S = np.reshape(S, (dim, cols(S)/dim))    
        rowmax = np.amax(S,axis=1)
        S = rowmax + np.log(np.dot(np.exp(S - rowmax), np.ones((10, 1))))
    
        return(sum(S)) 
   
    def FractionCommon(self):
        
        '''Takes in the Tree - with the last column reflecting the ordering of 
           the Names.'''
        '''Computes the fraction of the Tree time that each of the groups in 
           Splits data 
           were _the_same_ group'''

        SplitsArray=np.array(self.splittimes)
        Names1=np.array(SplitsArray[:,0]).flatten()
        Names2=np.array(SplitsArray[:,1]).flatten()

        NameInd1=np.where(self.name==Names1)[0]
        NameInd2=np.where(self.name==Names2)[0]
        SplitIndices=np.vstack((NameInd1,NameInd2)).T
        JointTimes=J(rows(SplitIndices),1,0)
        i=0
        for Split in SplitIndices:
            ToSum=np.where((self.resolvedtree[Split[0],:]==self.resolvedtree[Split[1],:]))[1]
            JointTimes[i,0]=np.nansum(self.filledtimeFractions[Split[1],ToSum])
            i+=1
        return np.hstack((SplitIndices,JointTimes))  
    
    def SplitLikelihood(self):
        
        ''' Uses FractionCommon and returns the log likelihood of the splits, given the
            prior mean and standard deviation in the file.'''
        
        TimeSplit = np.matrix(self.depth*1000*(1 - self.FractionCommon()[:,2]))
        Vals = (TimeSplit - self.splittimes[:,2])/self.splittimes[:,3]
        Vals = Vals.astype(float)
        
        return np.sum(norm.logpdf(Vals))      
    
    def DeathLikelihood(self, refpoint=2000):
        
        '''We need a method to try to match the time of expiry for those languages that went extinct.
           We dont have one yet but have the mata code somewhere.'''
        
        timeFracs = np.matrix(np.nansum(self.filledtimeFractions, axis=1)).T
        
        Place = np.where(self.deathmat[:,0]==0)[0]
        
        EDates = np.matrix(refpoint - self.depth*(1 - timeFracs)*1000)
        self.expirydates = EDates
        DateErr = (EDates[Place] - self.deathdata[Place,0])/self.deathdata[Place,1]
        DateErr = DateErr.astype(float)
        return np.sum(norm.logpdf(DateErr))
    
    def OriginLikelihood(self):
        
        '''A function that takes in a time-filled in Tree (Tr), a pairwise distance matrix D, the actual
           Tree in Matrix form (Tree), and a list of the row,column branch positions of the Tree. It puts
           out the likelihood of each row in Tree being the point of origin for the tree. '''
    
        '''Some things to watch out for - the np.where function can behave in funny ways! Also, it is critical
           that maxFinder and toAdd have the right shape - column vectors, effectively.'''
    
        '''In trouble shooting, a frequent culprit causing problems were programming snafus that
           left TM empty (equal to zero) while at the same time attempting to take a log.'''
    
        '''Another thing to watch out for - branches that are so small they violate Poisson rareness assumption!'''
    
    # Initialize placeholders:
        TM = np.copy(self.resolvedtree[:,-1])   #Just to shorten the code a bit
        bp = self.branchpositions      #Likewise
        LL = np.zeros((rows(TM), 1))
        IndCount = np.zeros((rows(TM), 1))
        Live = np.zeros((rows(TM), 1))        

        # Get rid of zeros in the Distance matrix:
        D = np.copy(self.D)
        D[D <= 0] = 10   
        lnD = np.matrix(-np.log(D))
        np.fill_diagonal(lnD,0) 

    # The main recursive loop backwards through the branches:
        for i in range(rows(self.resolvedtree), rows(bp)):     
            id = self.resolvedtree[bp[i,0], bp[i,1]]
            tu = np.where(self.resolvedtree[:,bp[i,1]] == id)[0]
            Bhat = self.depth*1000*np.copy(self.filledtimeFractions[tu,bp[i,1]:])
            TreeHat = np.copy(self.resolvedtree[tu,bp[i,1]:])
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
    
        return LL    
    
    def TotalLikelihood(self):
        '''Compute the combined likelihood of the tree, given death times
           parameters, word transition rates, the location model, and the
           prior information on splits.'''
        LL = self.DeathLikelihood() + self.jukescantorlikelihood() + \
                np.max(self.OriginLikelihood()) + self.SplitLikelihood()
        return(LL)
    
    def showtree(self):

        '''Prints out a visual display of a parameterized tree. Eventually
           should have an option to save as data.'''        
       
        TFR = self.depth*self.filledtimeFractions
        Tree = self.resolvedtree.astype(int)
        Num = np.matrix(np.arange(0, rows(Tree))).T
        names = np.array(self.name)

        # Some place holders for translating the tree into data
            
        minPos = J(rows(Tree), cols(Tree), np.nan)
        maxPos = J(rows(Tree), cols(Tree), np.nan)
        avePos = J(rows(Tree), cols(Tree), np.nan)
            
        # Roll through the tree and positions point and lines accordingly.
            
        for i in range(0, cols(Tree)):
            m=ps(Tree, i)
            for j in range(0, rows(m)):
                bracketed = psm(Num, m, j)
                lowest = min(bracketed)
                highest = max(bracketed)
                minPos[m[j,0]:m[j,1]+1,i:i+1] = J(1 + m[j,1] - m[j,0], 1, lowest)
                maxPos[m[j,0]:m[j,1]+1,i:i+1] = J(1 + m[j,1] - m[j,0], 1, highest)
                average = np.mean(bracketed)
                avePos[m[j,0]:m[j,1]+1,i:i+1] = J(1 + m[j,1] - m[j,0], 1, average)

        # Block of code to to compile tree running sums. 
            
        runningsum = np.copy(TFR[:,0:1])
        for i in range(1,cols(TFR)):
            runningsum=np.hstack((runningsum, np.copy(np.matrix(np.nansum(TFR[:,:i+1],axis=1)).T)))
            
        # Reorganize the numbers so that the x-axis is the depth of the tree
            
        depth = np.max(runningsum)
        runningsum = runningsum - depth

        # Plot the root of the tree...

        plt.plot([-depth,-depth + TFR[0,0]], [(rows(TFR) - 1)/2,(rows(TFR) - 1)/2])

        for c in range(0, cols(TFR) - 1):
            for r in range(0, rows(Tree)):
                y = [avePos[r,c], avePos[r,c+1]]
                x = [runningsum[r,c], runningsum[r,c+1]]
                plt.plot(x, y)
                plt.axis([-depth - 1,np.max(runningsum) + 1.5,-1,rows(Tree) + 1],yticks='none')
        
        names = list(names)        
        for n in range(0, len(names)):
            plt.annotate(names[n][0], xy = (runningsum[n,-1] + .2, n))      
            
    def RouteChooser(self):
        Tr = self.filledtimeFractions
        Tree = self.resolvedtree
        bp = self.branchpositions
        D = self.D

        branchRoute = J(0, 2, np.nan)   # Arrows
        TreeHat = np.hstack((Tree[:,0:-1], np.matrix(np.arange(0, rows(Tree))).T))        #Renumbered tree

        # Compute origin probability....
    
        lnProbs = self.OriginLikelihood()
    
        init = DiscretePicker(lnProbs)
        Path = J(rows(Tree), cols(Tree), np.nan)
        Path[:,0] = init   
    
        for i in reversed(range(rows(Tree), rows(bp) - 1)): 
            ind = bp[i,:]          # Retrieve an index
            z = ind[1] - 1
            found = 0
    
            groupId = TreeHat[ind[0], ind[1]]
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
                contProb = OriginLikelihood_args(subTree, subBr, subD, subTr)
                toGet = np.matrix(-np.log(D))[origin,:][:,colns].T
                contProb = contProb + toGet
                destination = colns[DiscretePicker(contProb,np.where(colns == origin)[0])]
            branchRoute = np.vstack(([origin,destination], branchRoute))
            colns = colns[np.where(colns != origin)[0]]
            Path[colns,ind[1]] = destination    
        
            # Cleaning up a few things from the route.. 
    
        for i in range(0, rows(Path)):
            if (np.any(branchRoute[:,1] == i)):
                pass
            else:
                s = cols(Path) - 2    #Not the last column, but the one before that!
                originFound = 0
                while True:
                    if Path[i,s] >= 0:               # Really just looking for np.nans!
                        originFound = 1
                    else: s -= 1
                    if originFound == 1:
                        if Path[i,s] != i:
                            branchRoute = np.vstack((branchRoute, (Path[i,s], i)))
                        break
        self.Path = Path
        self.branchRoute = branchRoute

    def TimeInPlace(self):

        '''Create estimates of time in place'''

        TimePlace = []
        for i in range(0, rows(self.Path)):
            firstAppearance = np.where(self.Path[i,:] == i)[0]
            if (len(firstAppearance) > 0):
                TimePlace.append(firstAppearance[0])
            else:
                TimePlace.append(-1)

        TimeInPlace = []
        for i in range(0, rows(self.Path)):
            toAdd = np.nansum(self.filledtimeFractions[i,TimePlace[i]:])
            TimeInPlace.append(toAdd)

        self.timeinplace = TimeInPlace

    def latlonplot(self, htmlfile, Cmap = "Reds"):

        '''
        In the future, we should add in a color choice method...or a save the data method.
        We can always export the previous data, however. This isn't working so 
        we might want to find a different way of doing this...Is there another
        way to export to a leaflet map?
        '''
        
        y = np.asarray(self.lat.astype(float)).flatten().tolist()
        x = np.asarray(self.lon.astype(float)).flatten().tolist()
        
        Order=np.asarray(self.resolvedtree[:,-1]).astype(int).flatten()

        # It might be better to reorder this elsewhere!
        
        y = [y[i] for i in Order]
        x = [x[i] for i in Order]
        
        name = np.asarray(self.name).flatten().tolist()

        s = [3000*n for n in self.timeinplace]

        myfig = plt.figure()
        mypic = myfig.add_subplot(111)
        mypic.scatter(x, y, c=s, s=s, cmap=Cmap)

        for i, txt in enumerate(name):
            mypic.annotate(name[i], (x[i], y[i]))
        
        mplleaflet.show(fig = myfig, path = htmlfile)


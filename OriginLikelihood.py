# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 09:09:11 2018

@author: Matthew
"""
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
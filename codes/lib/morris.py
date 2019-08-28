#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
"""
    Morris Method / Elementary Effects
    It includes optimised sampling of trajectories including optional groups
    as well as calculation of the Morris measures mu, stddev and mu*


    Definitions
    -----------
    def morris_sampling(NumFact, LB, UB, N=500, p=4, r=10, GroupMat=np.array([]), Diagnostic=0):
    def elementary_effects(NumFact, Sample, OutFact, Output, p=4, Group=[], Diagnostic=False):


    Input
    -----
    morris_sampling
        NumFact           number of factors examined
        LB                [NumFact] Lower Bound for each factor in list or array
        UB                [NumFact] Upper Bound for each factor in list or array

    elementary_effects
        NumFact           Number of factors
        Sample            Matrix of the Morris sampled trajectories
        OutFact           Matrix with the factor changings as specified in Morris sampling
        Output            Matrix of the output(s) values in correspondence of each point of each trajectory


    Optional Input
    --------------
    morris_sampling
        N                 Total number of trajectories (default: 500)
        p                 Number of levels (default: 4)
        r                 Final number of optimal trajectories (default: 10)
        GroupMat          [NumFact,NumGroups] Matrix describing the groups. (default: np.array([]))
                          Each column represents a group and its elements are set to 1 in correspondence
                          of the factors that belong to the fixed group. All the other elements are zero.
        Diagnostic        1=plot the histograms and compute the efficiency of the samplign or not,
                          0 otherwise (default)

    elementary_effects
        p                 Number of levels
        Group             [NumFactor, NumGroups] Matrix describing the groups.
                          Each column represents one group. The element of each column are zero
                          if the factor is not in the group. Otherwise it is 1.
        Diagnostic        True:  print out diagnostics
                          False: otherwise (default)

    Output
    ------
    morris_sampling
        [OptMatrix, OptOutVec]

    elementary_effects
        OutMatrix(NumFact*Output.shape[1], 3) = [Mu*, Mu, StDev]
            for each output it gives the three measures of each factor


    Notes
    -----
    The functions morris_sampling and elementary_effects are shortcuts for the functions
    Optimized_Groups and Morris_Measure_Groups of F. Campolongo and J. Cariboni ported to Python by Stijn Van Hoey.


    References
    ----------
    Saltelli, A., Chan, K., & Scott, E. M. (2000). Sensitivity Analysis.
        Wiley Series in Probability and Statistics, John Wiley & Sons, New York, 1-504. - on page 68ff


    Examples
    --------
    None


    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT License.

    Copyright (c) 2012-2017 S Van Hoey, Matthias Cuntz - mc (at) macu (dot) de

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.


    History
    -------
    Written original Matlab code by F. Campolongo, J. Cariboni, JRC - IPSC Ispra, Varese, IT
        Last Update: 15 November 2005 by J.Cariboni
        http://sensitivity-analysis.jrc.ec.europa.eu/software/index.htm
        now at: https://ec.europa.eu/jrc/en/samo/simlab
    Modified, S. Van Hoey, May 2012 - ported to Python
              MC, Oct 2013 - adapted to JAMS Python package and ported to Python 3
              MC, Dec 2017 - from exponential time increase with number of trajectories to linear increase
                           - by using in Optimised_Groups one call to cdist from scipy.spatial.distance
                             and removed one loop in a loop over total number of trajectories
                           - several little improvements on speed
"""
import numpy as np


def  Sampling_Function_2(p, k, r, LB, UB, GroupMat=np.array([])):
    """
        Morris sampling function


        Definition
        ----------
        def  Sampling_Function_2(p, k, r, LB, UB, GroupMat=np.array([])):


        Input
        -----
        p                               number of intervals considered in [0,1]
        k                               number of factors examined (sizea=k)
                                        In case groups are chosen the number of factors is stored in NumFact and sizea
                                        becomes the number of created groups (sizea=GroupMat.shape[1]).
        r                               sample size
        LB(sizea)                       Lower Bound for each factor in list or array
        UB(sizea)                       Upper Bound for each factor in list or array


        Optional Input
        --------------
        GroupMat(NumFact,GroupNumber)   Array which describes the chosen groups. (default: np.array([]))
                                        Each column represents a group and its elements are set to 1
                                        in correspondence of the factors that belong to the fixed group. All
                                        the other elements are zero.


        Output
        ------
        [Outmatrix, OutFact]
        Outmatrix(sizeb*r, sizea)       for the entire sample size computed In(i,j) matrices
        OutFact(sizea*r)                for the entire sample size computed Fact(i,1) vectors


        Notes
        -----
        Local Variables
            NumFact                     number of factors examined in the case when groups are chosen
            GroupNumber                 Number of groups (eventually 0)
            sizeb                       sizea+1
            randmult(sizea)             vector of random +1 and -1
            perm_e(sizea)               vector of sizea random permutated indeces
            fact(sizea)                 vector containing the factor varied within each traj
            DDo(sizea,sizea)            D*       in Morris, 1991
            A(sizeb,sizea)              Jk+1,k   in Morris, 1991
            B(sizeb,sizea)              B        in Morris, 1991
            Po(sizea,sizea)             P*       in Morris, 1991
            Bo(sizeb,sizea)             B*       in Morris, 1991
            Ao(sizeb)                   Jk+1     in Morris, 1991
            xo(sizea)                   x*       in Morris, 1991 (starting point for the trajectory)
            In(sizeb,sizea)             for each loop orientation matrix. It corresponds to a trajectory
                                        of k step in the parameter space and it provides a single elementary
                                        effect per factor
            Fact(sizea)                 for each loop vector indicating which factor or group of factors
                                        has been changed in each step of the trajectory
            AuxMat(sizeb,sizea)         Delta*0.5*((2*B - A) * DD0 + A) in Morris, 1991.
                                        The AuxMat is used as in Morris design for single factor analysis, while
                                        it constitutes an intermediate step for the group analysis.

        Note: B0 is constructed as in Morris design when groups are not considered.
              When groups are considered the routine follows the following steps:
                  1. Creation of P0 and DD0 matrices defined in Morris for the groups.
                     This means that the dimensions of these 2 matrices are (GroupNumber,GroupNumber).
                  2. Creation of AuxMat matrix with (GroupNumber+1,GroupNumber) elements.
                  3. Definition of GroupB0 starting from AuxMat, GroupMat and P0.
                  4. The final B0 for groups is obtained as [ones(sizeb,1)*x0' + GroupB0].
                     The P0 permutation is present in GroupB0 and it's not necessary to permute the
                     matrix (ones(sizeb,1)*x0') because it's already randomly created.


        References
        ----------
        Saltelli, A., Chan, K., & Scott, E. M. (2000). Sensitivity Analysis.
            Wiley Series in Probability and Statistics, John Wiley & Sons, New York, 1-504. - on page 68ff


        Examples
        --------
        None

        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT License.

        Copyright (c) 2012-2013 S Van Hoey, Matthias Cuntz - mc (at) macu (dot) de

        Permission is hereby granted, free of charge, to any person obtaining a copy
        of this software and associated documentation files (the "Software"), to deal
        in the Software without restriction, including without limitation the rights
        to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
        copies of the Software, and to permit persons to whom the Software is
        furnished to do so, subject to the following conditions:

        The above copyright notice and this permission notice shall be included in all
        copies or substantial portions of the Software.

        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
        LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
        OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
        SOFTWARE.


        History
        -------
        Written original Matlab code by F. Campolongo, J. Cariboni, JRC - IPSC Ispra, Varese, IT
            Last Update: 15 November 2005 by J.Cariboni
            http://sensitivity-analysis.jrc.ec.europa.eu/software/index.htm
            now at: https://ec.europa.eu/jrc/en/samo/simlab
        Modified, S. Van Hoey, May 2012 - ported to Python
                  MC, Oct 2013 - adapted to JAMS Python package and ported to Python 3
    """
    # Parameters and initialisation of the output matrix
    sizea = k
    Delta = p/(2.*(p-1.))
    NumFact = sizea
    if GroupMat.shape[0]==GroupMat.size:
        Groupnumber=0
    else:
        Groupnumber = GroupMat.shape[1]    #size(GroupMat,2)
        sizea = GroupMat.shape[1]

    sizeb = sizea + 1

    Outmatrix = np.zeros(((sizea+1)*r,NumFact))
    OutFact = np.zeros(((sizea+1)*r,1))
    # For each i generate a trajectory
    for i in range(r):
        Fact=np.zeros(sizea+1)
        # Construct DD0
        DD0 = np.matrix(np.diagflat(np.sign(np.random.random(k)*2-1)))

        # Construct B (lower triangular)
        B = np.matrix(np.tri((sizeb), sizea,k=-1, dtype=int))

        # Construct A0, A
        A0 = np.ones((sizeb,1))
        A = np.ones((sizeb,NumFact))

        # Construct the permutation matrix P0. In each column of P0 one randomly chosen element equals 1
        # while all the others equal zero.
        # P0 tells the order in which order factors are changed in each
        # Note that P0 is then used reading it by rows.
        I = np.matrix(np.eye(sizea))
        P0 = I[:,np.random.permutation(sizea)]

        # When groups are present the random permutation is done only on B. The effect is the same since
        # the added part (A0*x0') is completely random.
        if Groupnumber != 0:
            B = B * (np.matrix(GroupMat)*P0.transpose()).transpose()

        # Compute AuxMat both for single factors and groups analysis. For Single factors analysis
        # AuxMat is added to (A0*X0) and then permutated through P0. When groups are active AuxMat is
        # used to build GroupB0. AuxMat is created considering DD0. If the element on DD0 diagonal
        # is 1 then AuxMat will start with zero and add Delta. If the element on DD0 diagonal is -1
        # then DD0 will start Delta and goes to zero.
        AuxMat = Delta* 0.5 *((2*B - A) * DD0 + A)

        # a --> Define the random vector x0 for the factors. Note that x0 takes value in the hypercube
        # [0,...,1-Delta]*[0,...,1-Delta]*[0,...,1-Delta]*[0,...,1-Delta]
        # Original in Stijn Van Hoey's version
        #   xset=np.arange(0.0,1.0-Delta,1.0/(p-1))
        # Jule's version from The Primer
        #   xset=np.arange(0.0,1.0-1.0/(p-1),1.0/(p-1))
        # Matthias thinks that the difference between Python and Matlab is that Python is not taking
        # the last element; therefore the following version
        xset=np.arange(0.0,1.00000001-Delta,1.0/(p-1))
        x0 = np.matrix(xset.take(list(np.ceil(np.random.random(k)*np.floor(p/2))-1)))  #.transpose()

        # b --> Compute the matrix B*, here indicated as B0. Each row in B0 is a
        # trajectory for Morris Calculations. The dimension  of B0 is (Numfactors+1,Numfactors)
        if Groupnumber != 0:
            B0 = (A0*x0 + AuxMat)
        else:
            B0 = (A0*x0 + AuxMat)*P0

        # c --> Compute values in the original intervals
        # B0 has values x(i,j) in [0, 1/(p -1), 2/(p -1), ... , 1].
        # To obtain values in the original intervals [LB, UB] we compute
        # LB(j) + x(i,j)*(UB(j)-LB(j))
        In=np.tile(LB, (sizeb,1)) + np.array(B0)*np.tile((UB-LB), (sizeb,1)) #array!! ????

        # Create the Factor vector. Each component of this vector indicate which factor or group of factor
        # has been changed in each step of the trajectory.
        for j in range(sizea):
            Fact[j] = np.where(P0[j,:])[1]
        Fact[sizea] = int(-1)  #Enkel om vorm logisch te houden. of Fact kleiner maken

        #append the create traject to the others
        Outmatrix[i*(sizea+1):i*(sizea+1)+(sizea+1),:]=np.array(In)
        OutFact[i*(sizea+1):i*(sizea+1)+(sizea+1)]=np.array(Fact).reshape((sizea+1,1))

    return Outmatrix, OutFact



def Optimized_Groups(NumFact, LB, UB, N=500, p=4, r=10, GroupMat=np.array([]), Diagnostic=0):
    """
        Optimisation in the choice of trajectories for Morris experiment,
        that means elementary effects


        Definition
        ----------
        def Optimized_Groups(NumFact, LB, UB, N=500, p=4, r=10, GroupMat=np.array([]), Diagnostic=0):


        Input
        -----
        NumFact           Number of factors
        LB                [NumFact] Lower bound of the uniform distribution for each factor
        UB                [NumFact] Upper bound of the uniform distribution for each factor


        Optional Input
        --------------
        N                 Total number of trajectories (default: 500)
        p                 Number of levels (default: 4)
        r                 Final number of optimal trjectories (default: 10)
        GroupMat          [NumFact,NumGroups] Matrix describing the groups.  (default: np.array([]))
                          Each column represents a group and its elements are set to 1 in correspondence
                          of the factors that belong to the fixed group. All the other elements are zero.
        Diagnostic        1=plot the histograms and compute the efficiency of the samplign or not,
                          0 otherwise (default)


        Output
        ------
        [OptMatrix, OptOutVec]


        References
        ----------
        Saltelli, A., Chan, K., & Scott, E. M. (2000). Sensitivity Analysis.
            Wiley Series in Probability and Statistics, John Wiley & Sons, New York, 1-504. - on page 68ff


        Examples
        --------
        None

        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT License.

        Copyright (c) 2012-2013 S Van Hoey, Matthias Cuntz - mc (at) macu (dot) de

        Permission is hereby granted, free of charge, to any person obtaining a copy
        of this software and associated documentation files (the "Software"), to deal
        in the Software without restriction, including without limitation the rights
        to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
        copies of the Software, and to permit persons to whom the Software is
        furnished to do so, subject to the following conditions:

        The above copyright notice and this permission notice shall be included in all
        copies or substantial portions of the Software.

        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
        LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
        OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
        SOFTWARE.


        History
        -------
        Written original Matlab code by F. Campolongo, J. Cariboni, JRC - IPSC Ispra, Varese, IT
            Last Update: 15 November 2005 by J.Cariboni
            http://sensitivity-analysis.jrc.ec.europa.eu/software/index.htm
            now at: https://ec.europa.eu/jrc/en/samo/simlab
        Modified, S. Van Hoey, May 2012 - ported to Python
                  MC, Oct 2013 - adapted to JAMS Python package and ported to Python 3
    """
    from scipy.spatial import distance

    LBt = np.zeros(NumFact)
    UBt = np.ones(NumFact)

    # np.random.seed(seed=1025)
    OutMatrix, OutFact = Sampling_Function_2(p, NumFact, N, LBt, UBt, GroupMat)     #Version with Groups

    try:
        Groupnumber = GroupMat.shape[1]
    except:
        Groupnumber = 0

    if Groupnumber != 0:
        sizeb = Groupnumber + 1
    else:
        sizeb = NumFact + 1

    # Compute the distance between all pair of trajectories (sum of the distances between points)
    # The distance matrix is a matrix N*N
    # The distance is defined as the sum of the distances between all pairs of points
    #   if the two trajectories differ, 0 otherwise
    Dist = np.zeros((N,N))
    Diff_Traj = np.arange(0.0,N,1.0)
    MyDistall = distance.cdist(OutMatrix,OutMatrix)
    for j in range(N):   # combine all trajectories: eg N=3: 0&1; 0&2; 1&2 (is not dependent from sequence)
        for z in range(j+1,N):
            MyDist = MyDistall[sizeb*j:sizeb*(j+1),sizeb*z:sizeb*(z+1)]
            if np.where(MyDist==0.)[0].size == sizeb:
                # Same trajectory. If the number of zeros in Dist matrix is equal to
                # (NumFact+1) then the trajectory is a replica. In fact (NumFact+1) is the maximum number of
                # points that two trajectories can have in common
                Dist[j,z] = 0.
                Dist[z,j] = 0.
                # Memorise the replicated trajectory
                Diff_Traj[z] = -1.  #the z value identifies the duplicate
            else:
                # Define the distance between two trajectories as
                # the minimum distance among their points
                dd = np.sum(MyDist)
                Dist[j,z] = dd
                Dist[z,j] = dd

    # prepare array with excluded duplicates (alternative would be deleting rows)
    iidup = np.where(Diff_Traj==-1.)[0]
    dupli = iidup.size
    iiind = np.where(Diff_Traj!=-1.)[0]
    New_N = iiind.size # N - iidup.size
    New_OutMatrix = np.zeros((sizeb*New_N,NumFact))
    New_OutFact   = np.zeros((sizeb*New_N,1))

    # Eliminate replicated trajectories in the sampled matrix
    ID=0
    for i in range(N):
        if Diff_Traj[i] != -1.:
            New_OutMatrix[ID*sizeb:(ID+1)*sizeb,:] = OutMatrix[i*sizeb:(i+1)*sizeb,:]
            New_OutFact[ID*sizeb:(ID+1)*sizeb,:] = OutFact[i*sizeb:(i+1)*sizeb,:]
            ID+=1

    # Select in the distance matrix only the rows and columns of different trajectories
    #   Dist_Diff = np.delete(Dist_Diff,np.where(Diff_Traj==-1.)[0])
    Dist_Diff = Dist[iiind,:] #moet 2D matrix zijn... wis rijen ipv hou bij
    Dist_Diff = Dist_Diff[:,iiind] #moet 2D matrix zijn... wis rijen ipv hou bij

    # Select the optimal set of trajectories
    Traj_Vec = np.zeros((New_N,r), dtype=np.int)
    OptDist  = np.zeros((New_N,r))
    for m in range(New_N):                # each row in Traj_Vec
        Traj_Vec[m,0]  = m
        for z in range(1,r):              # elements in columns after first
            New_Dist_Diff = np.sqrt(np.sum(Dist_Diff[Traj_Vec[m,:z],:]**2, axis=0))
            ii            = New_Dist_Diff.argmax()
            Traj_Vec[m,z] = ii
            OptDist[m,z]  = New_Dist_Diff[ii]

    # Construct optimal matrix
    SumOptDist = np.sum(OptDist, axis=1)
    # Find the maximum distance
    Pluto = np.where(SumOptDist == SumOptDist.max())[0]
    Opt_Traj_Vec = Traj_Vec[Pluto[0],:]

    OptMatrix = np.zeros((sizeb*r,NumFact))
    OptOutVec = np.zeros((sizeb*r,1))

    for k in range(r):
        OptMatrix[k*sizeb:(k+1)*sizeb,:] = New_OutMatrix[sizeb*Opt_Traj_Vec[k]:sizeb*(Opt_Traj_Vec[k]+1),:]
        OptOutVec[k*sizeb:(k+1)*sizeb,:] = New_OutFact[sizeb*Opt_Traj_Vec[k]:sizeb*(Opt_Traj_Vec[k]+1),:]

    #----------------------------------------------------------------------
    # Compute values in the original intervals
    # Optmatrix has values x(i,j) in [0, 1/(p -1), 2/(p -1), ... , 1].
    # To obtain values in the original intervals [LB, UB] we compute
    # LB(j) + x(i,j)*(UB(j)-LB(j))
    OptMatrix_b = OptMatrix.copy()
    OptMatrix = np.tile(LB, (sizeb*r,1)) + OptMatrix*np.tile((UB-LB), (sizeb*r,1))

    if Diagnostic==True:
        # Clean the trajectories from repetitions and plot the histograms
        hplot=np.zeros((2*r,NumFact))

        for i in range(NumFact):
            for j in range(r):
                # select the first value of the factor
                hplot[j*2,i] = OptMatrix_b[j*sizeb,i]

                # search the second value
                for ii in range(1,sizeb):
                    if OptMatrix_b[j*sizeb+ii,i] != OptMatrix_b[j*sizeb,i]:
                        kk = 1
                        hplot[j*2+kk,i] = OptMatrix_b[j*sizeb+ii,i]

        import matplotlib.pyplot as plt
        fig=plt.figure()
        fig.suptitle('New Strategy')
        DimPlots = NumFact//2 + (1 if NumFact%2 > 0 else 0)
        for i in range(NumFact):
            ax=fig.add_subplot(DimPlots,2,i+1)
            ax.hist(hplot[:,i],p)

        # Plot the histogram for the original samplng strategy
        # Select the matrix
        OrigSample = OutMatrix[:r*(sizeb),:]
        print(OrigSample)
        Orihplot = np.zeros((2*r,NumFact))
        print(Orihplot)

        for i in range(NumFact):
            for j in range(r):
                # select the first value of the factor
                Orihplot[j*2,i] = OrigSample[j*sizeb,i]

                # search the second value
                for ii in range(1,sizeb):
                    if OrigSample[j*sizeb+ii,i] != OrigSample[j*sizeb,i]:
                        kk = 1
                        Orihplot[j*2+kk,i] = OrigSample[j*sizeb+ii,i]

        fig=plt.figure()
        fig.suptitle('Old Strategy')
        DimPlots = NumFact//2 + (1 if NumFact%2 > 0 else 0)
        for i in range(NumFact):
            ax=fig.add_subplot(DimPlots,2,i+1)
            ax.hist(Orihplot[:,i],p)
            #        plt.title('Old Strategy')
        print('hplotten')
        print(hplot)

        # Measure the quality of the sampling strategy
        levels=np.arange(0.0,1.1,1.0/(p-1))
        NumSPoint=np.zeros((NumFact,p))
        NumSOrigPoint=np.zeros((NumFact,p))
        for i in range(NumFact):
            for j in range(p):
                # For each factor and each level count the number of times the factor is on the level
                #This for the new and original sampling
                NumSPoint[i,j] = np.where(np.abs(hplot[:,i]-np.tile(levels[j], hplot.shape[0]))<1e-5)[0].size
                NumSOrigPoint[i,j] = np.where(np.abs(Orihplot[:,i]-np.tile(levels[j], Orihplot.shape[0]))<1e-5)[0].size

        # The optimal sampling has values uniformly distributed across the levels
        OptSampl = 2.*r/p
        QualMeasure = 0.
        QualOriMeasure = 0.
        for i in range(NumFact):
            for j in range(p):
                QualMeasure = QualMeasure + np.abs(NumSPoint[i,j]-OptSampl)
                QualOriMeasure = QualOriMeasure + np.abs(NumSOrigPoint[i,j]-OptSampl)

        QualMeasure = 1. - QualMeasure/(OptSampl*p*NumFact)
        QualOriMeasure = 1. - QualOriMeasure/(OptSampl*p*NumFact)

        print('The quality of the sampling strategy changed from {:f} with the old strategy to {:f} '
               'for the optimized strategy'.format(QualOriMeasure,QualMeasure))

        plt.show()

    return OptMatrix, OptOutVec[:,0]


def Morris_Measure_Groups(NumFact, Sample, OutFact, Output, p=4, Group=[], Diagnostic=False):
    """
        Given the Morris sample matrix, the output values and the group matrix compute the Morris measures.


        Definition
        ----------
        def Morris_Measure_Groups(NumFact, Sample, OutFact, Output, p=4, Group=[], Diagnostic=False):


        Input
        -----
        NumFact           Number of factors
        Sample            Matrix of the Morris sampled trajectories
        OutFact           Matrix with the factor changings as specified in Morris sampling
        Output            Matrix of the output(s) values in correspondence of each point of each trajectory


        Optional Input
        --------------
        p                 Number of levels
        Group             [NumFact, NumGroups] Matrix describing the groups.
                          Each column represents one group. The element of each column are zero
                          if the factor is not in the group. Otherwise it is 1.
        Diagnostic        True:  print out diagnostics
                          False: otherwise (default)


        Output
        ------
        OutMatrix(NumFact*Output.shape[1], 3) = [Mu*, Mu, StDev]
            for each output it gives the three measures of each factor


        References
        ----------
        Saltelli, A., Chan, K., & Scott, E. M. (2000). Sensitivity Analysis.
            Wiley Series in Probability and Statistics, John Wiley & Sons, New York, 1-504. - on page 68ff


        Examples
        --------
        None

        License
        -------
        This file is part of the JAMS Python package, distributed under the MIT License.

        Copyright (c) 2012-2017 S Van Hoey, Matthias Cuntz - mc (at) macu (dot) de

        Permission is hereby granted, free of charge, to any person obtaining a copy
        of this software and associated documentation files (the "Software"), to deal
        in the Software without restriction, including without limitation the rights
        to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
        copies of the Software, and to permit persons to whom the Software is
        furnished to do so, subject to the following conditions:

        The above copyright notice and this permission notice shall be included in all
        copies or substantial portions of the Software.

        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
        LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
        OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
        SOFTWARE.


        History
        -------
        Written original Matlab code by F. Campolongo, J. Cariboni, JRC - IPSC Ispra, Varese, IT
            Last Update: 15 November 2005 by J.Cariboni
            http://sensitivity-analysis.jrc.ec.europa.eu/software/index.htm
            now at: https://ec.europa.eu/jrc/en/samo/simlab
        Modified, S. Van Hoey, May 2012 - ported to Python
                  MC, Oct 2013 - adapted to JAMS Python package and ported to Python 3
                  MC, Dec 2017 - deal with single trajectories
                  MC, Dec 2017 - set Delta=A and not A[np.where(A)]; deal with Delta==0.
    """

    try:
        NumGroups = Group.shape[1]
        if Diagnostic: print('{:d} Groups are used'.format(NumGroups))
    except:
        NumGroups = 0
        if Diagnostic: print('No Groups are used')

    Delt = p/(2.*(p-1.))

    if NumGroups != 0:
        sizea = NumGroups
        GroupMat = Group
        GroupMat = GroupMat.transpose()
        if Diagnostic: print(NumGroups)
    else:
        sizea = NumFact
    sizeb = sizea+1

    r = Sample.shape[0]/sizeb

    try:
        NumOutp = Output.shape[1]
    except:
        NumOutp = 1
        Output=Output.reshape((Output.size,1))


    # For each Output
    if NumGroups == 0:
        OutMatrix=np.zeros((NumOutp*NumFact,3)) # for every output: every factor is a line, columns are mu*, mu and std
    else:
        OutMatrix=np.zeros((NumOutp*NumFact,1)) # for every output: every factor is a line, column is mu*

    SAmeas_out=np.zeros((NumOutp*NumFact,r))

    for k in range(NumOutp):
        OutValues = Output[:,k]

        # For each trajectory
        SAmeas = np.zeros((NumFact,r)) #vorm afhankelijk maken van group of niet...
        for i in range(r):
            # For each step j in the trajectory
            # Read the orientation matrix fact for the r-th sampling
            # Read the corresponding output values
            # Read the line of changing factors
            Single_Sample    = Sample[i*sizeb:(i+1)*sizeb,:]
            Single_OutValues = OutValues[i*sizeb:(i+1)*sizeb]
            Single_Facts     = np.array(OutFact[i*sizeb:(i+1)*sizeb], dtype=np.int) # gives factor in change (or group)

            A = (Single_Sample[1:sizeb,:]-Single_Sample[:sizea,:]).transpose()
            Delta=A[np.where(A)] # AAN TE PASSEN?

            if Diagnostic:
                print('A: ', A)
                print('Delta: ', Delta)
                print('Single_Facts: ', Single_Facts)

            # For each point of the fixed trajectory, i.e. for each factor, compute the values of the Morris function.
            for j in range(sizea):
                if NumGroups != 0:  #work with groups
                    Auxfind=A[:,j]
                    Change_factor = np.where(np.abs(Auxfind)>1e-010)[0]
                    for gk in Change_factor:
                        SAmeas[gk,i] = np.abs((Single_OutValues[j] - Single_OutValues[j+1])/Delt)   #nog niet volledig goe
                else:
                    if Delta[j] > 0.0:
                        SAmeas[Single_Facts[j],i] = (Single_OutValues[j+1] - Single_OutValues[j])/Delt
                    else:
                        SAmeas[Single_Facts[j],i] = (Single_OutValues[j] - Single_OutValues[j+1])/Delt

        # Compute Mu AbsMu and StDev
        if np.isnan(SAmeas).any():
            AbsMu = np.zeros(NumFact)
            Mu    = np.zeros(NumFact)
            Stdev = np.zeros(NumFact)

            for j in range(NumFact):
                SAm = SAmeas[j,:]
                SAm = SAm[~np.isnan(SAm)]
                rr  = np.float(SAm.size)
                AbsMu[j] = np.sum(np.abs(SAm))/rr
                if NumGroups == 0:
                    Mu[j] = SAm.mean()
                    if SAm.size > 1:
                        Stdev[j] = np.std(SAm, dtype=np.float64, ddof=1) #ddof: /N-1 instead of /N
                    else:
                        Stdev[j] = 0.
        else:
            AbsMu = np.sum(np.abs(SAmeas), axis=1)/r
            if NumGroups == 0:
                Mu = SAmeas.mean(axis=1)
                if SAmeas.shape[1] > 1:
                    Stdev = np.std(SAmeas, dtype=np.float64, ddof=1, axis=1) #ddof: /N-1 instead of /N
                else:
                    Stdev = np.zeros(SAmeas.shape[0])
            else:
                Mu    = np.zeros(NumFact)
                Stdev = np.zeros(NumFact)

        OutMatrix[k*NumFact:k*NumFact+NumFact,0] = AbsMu
        if NumGroups == 0:
            OutMatrix[k*NumFact:k*NumFact+NumFact,1] = Mu
            OutMatrix[k*NumFact:k*NumFact+NumFact,2] = Stdev

        SAmeas_out[k*NumFact:k*NumFact+NumFact,:] = SAmeas


    return SAmeas_out, OutMatrix


def morris_sampling(NumFact, LB, UB, N=500, p=4, r=10, GroupMat=np.array([]), Diagnostic=0):
    """
        Wrapper function for Optimized_Groups.
        def Optimized_Groups(NumFact, LB, UB, N=500, p=4, r=10, GroupMat=np.array([]), Diagnostic=0):
    """
    return Optimized_Groups(NumFact, LB, UB, N, p, r, GroupMat, Diagnostic)


def elementary_effects(NumFact, Sample, OutFact, Output, p=4, Group=[], Diagnostic=False):
    """
        Wrapper function for Morris_Measure_Groups.
        def Morris_Measure_Groups(NumFact, Sample, OutFact, Output, p=4, Group=[], Diagnostic=False):
    """
    return Morris_Measure_Groups(NumFact, Sample, OutFact, Output, p, Group, Diagnostic)


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # NumFact = 15
    # LB = np.arange(NumFact)
    # UB = 2.*LB + 1.
    # N = 100
    # p = 6
    # r = 10
    # Diagnostic = 0
    # out = np.random.random(r*(NumFact+1))
    # mat, vec = morris_sampling(NumFact, LB, UB, N=N, p=p, r=r, Diagnostic=Diagnostic)
    # sa, res = elementary_effects(NumFact, mat, vec, out, p=p)
    # print(res) # (NumFact,3) = AbsMu, Mu, Stddev
    # print(sa) #  (NumFact,r) individual elementary effects for all parameters

    # NumFact = 15
    # LB = np.arange(NumFact)
    # UB = 2.*LB + 1.
    # N = 10
    # p = 6
    # r = 1
    # Diagnostic = 0
    # out = np.random.random(r*(NumFact+1))
    # mat, vec = morris_sampling(NumFact, LB, UB, N=N, p=p, r=r, Diagnostic=Diagnostic)
    # sa, res = elementary_effects(NumFact, mat, vec, out, p=p)
    # print(res) # (NumFact,3) = AbsMu, Mu, Stddev
    # print(sa) #  (NumFact,1)

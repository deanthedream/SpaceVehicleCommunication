"""
A 2d discrete distribution sampling method for a discrete array of bins
This function turns the 2d array into a 1d array, normalizes the bin counts,
calculates the CDF, and uses an inverse transform sampler to generate samples from this distribution

Attributes:
    xedges (numpy array) - a 1d numpy array of the bin xedges
    yedges (numpy array) - a 1d numpy array of the bin yedges
    counts (2d numpy array) - a 2d numpy array containing bin counts (or normalized bin counts)
        with size (len(xedges)-1) by (len(yedges) - 1)

Written by: Dean Keithly
Written On: 7/1/2019
"""
import numpy as np
#from EXOSIMS.util.InverseTransformSampler import InverseTransformSampler
import matplotlib.pyplot as plt

class sampleDistribution:
    def __init__(self, counts, xedges=None, yedges=None, xcent=None, ycent=None):
        # init constructs the sampler
        if not xcent is None: #prefer xcent to xedges
            self.xcent = xcent
            self.xmin = np.min(self.xcent) # minimum x value
            self.xmax = np.max(self.xcent) # maximum x value
        elif not xedges is None:
            self.xedges = xedges
            self.xmin = np.min(self.xedges) # minimum x value
            self.xmax = np.max(self.xedges) # maximum x value
            self.xcent = 0.5*(self.xedges[1:]+self.xedges[:-1])
        else:
            raise Exception

        if not ycent is None:
            self.ycent = ycent
            self.ymin = np.min(self.ycent) # minimum y value
            self.ymax = np.max(self.ycent) # maximum y value
        elif not yedges is None:
            self.yedges = yedges
            self.ymin = np.min(self.yedges) # minimum y value
            self.ymax = np.max(self.yedges) # maximum y value
            self.ycent = 0.5*(self.yedges[1:]+self.yedges[:-1])
        else:
            raise Exception
        
        self.counts = counts
        self.numCounts = np.sum(self.counts) # the total number of counts

        assert len(self.xcent) == len(counts[0,:]), 'counts has invalid x dimension, check xedges in'
        assert len(self.ycent) == len(counts[:,0]), 'counts has invalid y dimension, check xedges in'

        # concatenate the x and y into single vector
        self.scent = np.arange(len(self.xcent)*len(self.ycent))

        # probability mass function of value being in each x,y bin as indexed by s
        self.indexedPMF = [self.counts[int(np.floor(self.scent[i]/len(self.xcent))),np.mod(self.scent[i],len(self.xcent))] for i in np.arange(len(self.scent))]# for each scent, this is the probability of being that s
        assert len(self.xcent)*len(self.ycent), 'indexedPDF is not correct size'
        
        # cumulative mass function function of being below x,y bin as indexed by s
        self.indexedCMF = np.cumsum(self.indexedPMF)/self.numCounts # for each scent, this is the cumulative distribution of being below that s normalized to 1

    def plotPMF(self):
        plt.close(65489651)
        fig = plt.figure(65489651)
        plt.plot(self.scent,self.indexedPMF)
        plt.title('indexedPMF')
        plt.xlabel('s edges')
        plt.ylabel('Probability Mass Function of s')
        plt.show(block=False)

    def plotCMF(self):
        plt.close(983173)
        fig = plt.figure(983173)
        plt.plot(self.scent,self.indexedCMF)
        plt.title('indexedCMF')
        plt.xlabel('s edges')
        plt.ylabel('Cumulative Mass Function of s')
        plt.show(block=False)

    def plot2dPMF(self,x,y):
        """ Creates a 2d Scatter plot of sampled data
        x (numpy array) - x values
        y (numpy array) - y values
        """
        plt.close(65183)
        fig = plt.figure(65183)
        plt.scatter(x, y, alpha=0.1)
        plt.xlabel('x cent')
        plt.ylabel('y cent')
        plt.show(block=False)


    def discreteInverseTransformSampler(self,nsamples):
        """
        Args:
            nsamples (integer) - number of samples to generate
        Returns:
            xind (numpy array) - x index of each sample with len(nsamples)
            yind (numpy array) - y index of each sample with len(nsamples)
            x (numpy array) - x value of each sample with len(nsamples)
            y (numpy array) - y value of each sample with len(nsamples)
        """
        nsamples = int(nsamples)
        # Generate random numbers between 0 and 1
        Ps = np.random.uniform(low=0., high=1.,size=nsamples) #
        
        self.sInds = np.asarray([np.where(Ps[i]<=self.indexedCMF)[0][0] for i in np.arange(nsamples)])
        xind = np.asarray([int(np.floor(self.sInds[i]/len(self.xcent))) for i in np.arange(nsamples)])
        #yind = np.asarray([self.sInds - np.mod(self.sInds[i],len(self.xcent)) for i in np.arange(nsamples)], dtype=int)
        yind = np.asarray([np.mod(self.sInds[i],len(self.xcent)) for i in np.arange(nsamples)])
        x = np.asarray([self.xcent[i] for i in xind], dtype=float)
        y = np.asarray([self.ycent[i] for i in yind], dtype=float)

        return xind, yind, x, y

    def discreteInverseTransformSampler_given_1DCMF(self, cmf, xcent, nsamples):
        """ This function samples a 1D CMF for nsamples
        Args:
            cmf (numpy array) - cumulative probability of x being below index
            xcent (numpy array) - x values of each cmf value
            nsamples (integer) - number of samples to generate
        Returns:
            xvalues (numpy array) - the x values generated
        """
        Ps = np.random.uniform(low=0., high=1.,size=nsamples) # generates uniform random numbers
        cmfInds = np.asarray([np.where(Ps[i]<=cmf)[0][0] for i in np.arange(nsamples)]) # calculates the cmf inds
        xvalues = np.asarray([xcent[cmfInd] for cmfInd in cmfInds])
        return xvalues


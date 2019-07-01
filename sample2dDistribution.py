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

class sample2dDistribution:
    def __init__(self, xedges, yedges, counts):
        # init constructs the sampler
        self.xedges = xedges
        self.xmin = np.min(self.xedges) # minimum x value
        self.xmax = np.max(self.xedges) # maximum x value
        self.xcent = 0.5*(self.xedges[1:]+self.xedges[:-1])
        self.yedges = yedges
        self.ymin = np.min(self.yedges) # minimum y value
        self.ymax = np.max(self.yedges) # maximum y value
        self.ycent = 0.5*(self.yedges[1:]+self.yedges[:-1])
        self.counts = counts
        self.numCounts = np.sum(self.counts) # the total number of counts

        assert len(xedges) - 1 == len(counts[0,:]), 'counts has invalid x dimension'
        assert len(yedges) - 1 == len(counts[:,0]), 'counts has invalid y dimension'

        # concatenate the x and y into single vector
        self.scent = np.arange(len(xcent)*len(ycent))

        # probability mass function of value being in each x,y bin as indexed by s
        self.indexedPMF = [self.counts[np.floor(self.scent[i]/len(self.xcent)),self.scent[i] - np.mod(self.scent[i],len(self.xcent))] for i in np.arange(len(self.scent))]# for each scent, this is the probability of being that s
        assert len(self.xcent)*len(self.ycent), 'indexedPDF is not correct size'
        
        # cumulative mass function function of being below x,y bin as indexed by s
        self.indexedCMF = np.cumsum(indexPMF)# for each scent, this is the cumulative distribution of being below that s

    def plotPMF(self):
        plt.close(65489651)
        fig = plt.figure(65489651)
        plt.plot(self.scent,self.indexedPMF)
        plt.xlabel('s edges')
        plt.ylabel('Probability Mass Function of s')
        plt.show(block=False)

    def plotCMF(self):
        plt.close(983173)
        fig = plt.figure(983173)
        plt.plot(self.scent,self.indexedCMF)
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
        plt.plot(x, y, alpha=0.1)
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
        # Generate random numbers between 0 and 1
        Ps = np.random.uniform(low=0., high=1.,size=nsamples) #
        
        self.sInds = np.asarray([np.where(Ps[i]<=self.indexedCMF)[0][0] for i in np.arange(len(nsamples))])
        xind = np.asarray([np.floor(self.sInds[i]/len(self.xcent)) for i in np.arange(len(nsamples))])
        yind = np.asarray([self.sInds - np.mod(self.sInds[i],len(self.xcent)) for i in np.arange(len(nsamples))])
        x = self.xcent[xind]
        y = self.ycent[yind]

        return xind, yind, x, y




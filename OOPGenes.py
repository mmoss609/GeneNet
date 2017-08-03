#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

An incredibly inneficient script to create a gene regulatory net from 
RNA-seq data. We'' call this version 1 for now. Who knows if I'll choose to fix it. 
We'll see if it works/ if Camilla wants me to make it long lasting for future analyses

Made more efficient by adding parallel capabilities if those work

Created on Mon Jul 31 12:32:24 2017

@author: moss
"""

'''Import necessary libraries'''
'''These are mostly just needed for i/o at start of script'''
import pandas as pd #For excel. Going to use for gene names. Skipping for now for ease of programming
import os #For setting working directory
'''These are important for functional aspects of code'''
import numpy as np 
import scipy.stats as spt
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import multiprocessing

os.chdir('/Users/moss/Documents')
num_cores = multiprocessing.cpu_count()

'''Creates gene object that for now just has the different experiments/replicates as input parameters, but will also have the name soon'''
class Gene:
    def __init__(self,dataFile):
        '''Defines characteristics of each gene object'''
        #self.name = name
        self.data = dataFile
        self.reps = []
    
    def geneArrays(self,geneReps):
        for j in range(1,len(np.transpose(self.data))):
            Gene.reps.append(data[j])
        
def Correlate(Gene1,Gene2):
    '''Function for finding spearman correlation coefficient between 2 gene objects. Used in a parallel loop later'''

    cor = spt.spearmanr(gene1,gene2)

    return cor[1]

#data = pd.read_csv('FullGeneListwReplicates.csv',sep=',',header=0,usecols = [1,2,3,4,5,6])
        
'''Current main body of program. Imports data and then creates gene objects for each gene in the data set'''
data = np.genfromtxt('FullGeneListwReplicates.csv',delimiter=',')
genes = []
for i in range(1,len(data)):
    genes.append(Gene.geneArrays(data[i]))
    #genes.append(Gene(data[i,1],data[i,2],data[i,3],data[i,4],data[i,5],data[i,6]))

'''Runs correlation function for all gene pairs in dataset. This is the computationally intensive part, as it first allocates the space for the large correlation matrix,
    #then calculates the gene-gene correlation, and then puts those in the proper space in the correlation matrix. Later, I'll be able to associate this with gene names
    #for now, just indeces'''
corr = np.zeros((len(data),len(data)))

itsI = 0
itsJ = 0
'''for itsI in range(len(data)-1):
    if np.mod(itsI,50) == 0:
                print('I: {0}'.format(its))
    for itsJ in range(len(data)-1):
        if np.mod(itsJ,50) == 0:
                print('J: {0}'.format(itsJ))
        corr[itsI,itsJ]=Gene.Correlate(genes[itsI],genes[itsJ])'''

'''Performs the computations but runs it in parallel. Hopefully substantially faster than original function'''

corr = Parallel(n_jobs = num_cores,backend='threading')(delayed(Correlate)(genes[itsI],genes[itsJ]) for itsI in range(len(data)-1) for itsJ in range(len(data)-1))

'''Copy that correlation matrix and then set a threshold that removes all connections with coefficient -0.9 <= p <= 0.9
    (considers them as noise) to aid in viewing connectivity after. I can definitely build this into the for loop that calculates
    the corr coefficients if i want to save memory and don't care about the actual connection strength above or below the thresholds'''

for row in range(len(corr)):
    for col in range(len(corr)):
        if corr[row][col] >= 0.9:
            corr[row][col] = 1
        elif corr[row][col] <= -0.9:
            corr[row][col] = 1
        else:
            corr[row][col] = 0


'''Download biopython and use that for the clustering?'''


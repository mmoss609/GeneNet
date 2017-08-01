# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 14:59:50 2017

This is probably just gonna be a script dump where i comment 
out whatever I'm not using. We'll see though. For now
it has a correlation function that doesn't work and a 
script to remove duplicates from a list of strings

@author: matt
"""

import numpy as np
import csv
import pandas as pd


x = open('AllRNACopy.csv')
y = open('AllImmuneGenes.csv')
geneList = map(str.upper,y.readlines())
geneList = sorted(geneList)
x = map(str.upper,x.readlines())
#y = np.correlate(x[:][2],x[:][3],"full")
#y = []
#
#for i in xrange(len(geneList[1])):
#    for j in xrange(len(geneList[2])):
#        y.extend(np.correlate(geneList[i][1],geneList[j+1][2]))
newList=[]
for i in range(len(geneList)-1):
    if geneList[i] == geneList[i+1]:
        pass
    else:
        newList.append(geneList[i])
commonGenes=[]
for i in range(len(x)):
    for j in range(len(newList)):
        if x[i] == newList[j]:
            commonGenes.append(newList[j])

commonGenes = sorted(commonGenes)
df = pd.DataFrame(commonGenes)
df.to_csv('AllImmuneGenesSorted.csv',index=False)

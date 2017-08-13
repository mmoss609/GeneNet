# -*- coding: utf-8 -*-
"""
Created on Sat Aug 12 22:33:35 2017

@author: Moss
"""
import os
import pandas as pd
from numpy import inf

os.chdir('C:\\Users\\Moss\\Downloads')
x = pd.read_table('DMSO-hdac.tabular',header=0)
names = x['id'].astype('str')
folds = x['log2FoldChange']
'''for j in range(len(folds)):
    if folds[j] == '-inf':
        folds[j] = 0
    elif folds[j] == 'inf':
        folds[j] = 0'''

folds[folds==-inf] = 0
folds[folds==inf] = 0
for i in range(len(names)):
    names[i] = names[i].upper()
    
x.log2FoldChange = folds
x.id = names

x = x.sort_values('pval',ascending=True)
d = {'ID':x.id,'LOG2FOLDCHANGE':x.log2FoldChange}
df = pd.DataFrame(data=d)
filename = 'InputForGSEA.rnk'
df.to_csv(filename,sep='\t',index=False,encoding='utf-8')
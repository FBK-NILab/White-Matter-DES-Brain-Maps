#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 15:47:47 


@author: ludovicocoletta
"""

import pandas as pd
import numpy as np
from scipy.stats import ttest_rel
from scipy.stats import mannwhitneyu
from pingouin import compute_effsize # Just to double check
import os

def cohen_d_dep(a,b):
    
    # As per 
    #Lakens, D., 2013. Calculating and reporting effect sizes to facilitate cumulative science: a practical primer for t-tests and ANOVAs. 
    #Front. Psychol. 4, 863. https://doi.org/10.3389/fpsyg.2013.00863
    
    mean_a=a.mean()
    var_a=a.var(ddof=1)
    mean_b=b.mean()
    var_b=b.var(ddof=1)
    
    num=mean_a-mean_b
    den=np.sqrt(((var_a+var_b))/2.)
    
    return num/den
    
def main():
    
    func='PHONOLOGICAL'
    
    X=pd.read_csv(os.path.join('maps_2mm_cleaned',func+'.csv'),index_col=0,header='infer').to_numpy().astype('float') # change this path
    
    stat=ttest_rel(X[:,0],X[:,1],alternative='greater') # firts col is positive netw, second column anticorrelated netw
    print(func,X[:,0].mean(),X[:,1].mean(),cohen_d_dep(X[:,0], X[:,1]),stat)
    
if __name__ == "__main__":
        main() 

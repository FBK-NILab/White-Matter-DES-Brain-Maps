#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 18:17:29 2023

@author: ludovicocoletta
"""

import numpy as np
import pandas as pd
from sklearn.metrics import r2_score
from sklearn.model_selection import LeaveOneOut
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import ShuffleSplit
from sklearn.linear_model import QuantileRegressor
from sklearn.preprocessing import StandardScaler
import time

def main():
    
    df=pd.read_csv('intersection_with_hubs.csv')
    
    X=df.iloc[:,0].to_numpy().astype(float).reshape((-1,1)) # df.iloc[:,1] for phono and df.iloc[:,5] for Tot lesion volume

    y=df.iloc[:,6].to_numpy().astype(float)
    
    loo=list(LeaveOneOut().split(X,y))
    
    r2_values=np.zeros(100)
    
    pred_across_trials=[None]*100
    
    for it_index, it in enumerate(r2_values):
    
        print(it_index)
        predictions_by_subjects=np.zeros((len(df)))
        
        for split_index, split in enumerate(loo):
            
            
            pipeline_1 = Pipeline([
                ('norm', StandardScaler()),
                ("regr", QuantileRegressor(solver='highs-ds'))
                        ])   
        
            param_grid_1 = dict(
                regr__quantile = np.linspace(0.05,0.95,19),
                regr__alpha=np.logspace(-5,5,10)          
            )        
            
    
            grid_search = GridSearchCV(pipeline_1, param_grid=param_grid_1, scoring='r2',verbose=0,cv=ShuffleSplit(n_splits=10,test_size=0.5),n_jobs=4)
            grid_search.fit(X[split[0]],y[split[0]])
            predictions_by_subjects[split_index]=grid_search.predict(X[split[1]])

        r2_values[it_index]=r2_score(y,predictions_by_subjects)
        pred_across_trials[it_index]=predictions_by_subjects
        time.sleep(8)
        
    np.save('DES_SEMANTIC_MODEL_PREDS.npy', pred_across_trials)
        
if __name__ == "__main__":
    main() 

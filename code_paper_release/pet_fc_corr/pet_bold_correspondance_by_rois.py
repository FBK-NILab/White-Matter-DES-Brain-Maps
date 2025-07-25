#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 16:16:06 2023

@author: ludovicocoletta
"""

from brainsmash.mapgen.base import Base
from brainsmash.mapgen.eval import base_fit
from brainsmash.mapgen.stats import nonparp
import glob
import numpy as np
import time
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr,pearsonr
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import fdrcorrection

def ccc(x,y):
    
    ''' Concordance Correlation Coefficient.
    Taken from: https://nirpyresearch.com/concordance-correlation-coefficient/
    '''
    sxy = np.sum((x - x.mean())*(y - y.mean()))/x.shape[0]
    rhoc = 2*sxy / (np.var(x) + np.var(y) + (x.mean() - y.mean())**2)
    return rhoc

def main():

    n_shuffle=10000
    thr_perc=[0]
    
    rois=sorted(glob.glob('JHU-ICBM-labels-2mm_cog/roi_*txt'))
    
    cog=np.asarray([np.loadtxt(ii) for ii in rois])
    
    dist_matrix=squareform(pdist(cog, metric='euclidean'))

    
    sba_files=sorted(glob.glob('sba/*npy'))
    pet_files=sorted(glob.glob('pet/*npy'))
        
    conc_pet=np.zeros((dist_matrix.shape[1],dist_matrix.shape[1],len(pet_files)))
    conc_sba=np.zeros((dist_matrix.shape[1],dist_matrix.shape[1],len(sba_files)))
    
    for sub_index, sub_file in enumerate(sba_files):
        dummy_ts=np.load(sub_file)
        dummy_fc_matrix=np.corrcoef(dummy_ts.T)
        np.fill_diagonal(dummy_fc_matrix,0)
        conc_sba[:,:,sub_index]=np.arctanh(dummy_fc_matrix)
    mean_fc_sba=np.tanh(conc_sba.mean(axis=2))
    
    for sub_index, sub_file in enumerate(pet_files):
        dummy_ts=np.load(sub_file)
        dummy_fc_matrix=np.corrcoef(dummy_ts.T)
        np.fill_diagonal(dummy_fc_matrix,0)
        conc_pet[:,:,sub_index]=np.arctanh(dummy_fc_matrix)
    mean_mc_pet=np.tanh(conc_pet.mean(axis=2))
    
    whole_brain=1
    
    pvals_across_rois=[]
    sim_across_rois=[]
    null_model_fit_across_rois=[]
    
    for roi in range(0,dist_matrix.shape[0]):
        
        roi_of_int=roi
        bold_map=mean_fc_sba[roi_of_int,:]
        pet_map=mean_mc_pet[roi_of_int,:]
        dist_matrix_copy=dist_matrix.copy()
        
        if (roi_of_int%2==0) and (whole_brain)==0:
            rois=sorted(list(set(list(range(0,6))+list(range(0,dist_matrix.shape[0],2)))))
            ixgrid = np.ix_(rois,rois)
            dist_matrix_copy=dist_matrix_copy[ixgrid]
            bold_map=bold_map[ixgrid[0]]
            pet_map=pet_map[ixgrid[0]]
        elif (roi_of_int%2==1) and (whole_brain)==0:
            rois=sorted(list(set(list(range(0,6))+list(range(1,dist_matrix.shape[0],2)))))
            ixgrid = np.ix_(rois,rois)
            dist_matrix_copy=dist_matrix_copy[ixgrid]
            bold_map=bold_map[ixgrid[0]]
            pet_map=pet_map[ixgrid[0]]
        else:
            pass
            
        
        sim_across_thr=np.zeros((len(thr_perc),2))
        pvals_across_thr=np.zeros((len(thr_perc),2))
        null_model_fit=np.zeros((len(thr_perc)))
        
        for thr_index,thr in enumerate(thr_perc):
            
            print(thr_index)
            pet_map_copy=pet_map.copy().flatten()
            pet_map_copy[pet_map_copy<np.percentile(pet_map_copy,thr)]=0
            
            bold_map_copy=bold_map.copy().flatten()
            bold_map_copy[bold_map_copy<np.percentile(bold_map_copy,thr)]=0
            
            keywords={'resample': True}
            emp_var,u0,surr_var=base_fit(bold_map_copy, dist_matrix_copy, nsurr=1000,return_data=True,**keywords)
            
            #null_model_fit[thr_index]=np.percentile(([(spearmanr(emp_var,surr_var[iii,:])[0]) for iii in range(0,surr_var.shape[0])]),[50])[0]
            null_model_fit[thr_index]=pearsonr(emp_var,surr_var.mean(axis=0))[0]
        
            base=Base(x=bold_map_copy,D=dist_matrix_copy,resample=True,n_jobs=-1)                       
            surrogates=base(n=n_shuffle)
            
            dist=[pearsonr(surrogates[iii,:], pet_map_copy)[0] for iii in range(0,n_shuffle)]
            stat=pearsonr(bold_map_copy,pet_map_copy)[0]
            sim_across_thr[thr_index,0]=stat
            pvals_across_thr[thr_index,0]=1-(np.where(dist<stat)[0].shape[0]/n_shuffle)
        
            dist=[spearmanr(surrogates[iii,:], pet_map_copy)[0] for iii in range(0,n_shuffle)]
            stat=spearmanr(bold_map_copy,pet_map_copy)[0]
            sim_across_thr[thr_index,1]=stat
            pvals_across_thr[thr_index,1]=1-(np.where(dist<stat)[0].shape[0]/n_shuffle)
            
            del(base,surrogates,dist)
        
        pvals_across_rois.append(pvals_across_thr)
        null_model_fit_across_rois.append(null_model_fit)
        sim_across_rois.append(sim_across_thr)
        
    pearson_pvals=np.array(pvals_across_rois)
    fit=np.array(null_model_fit_across_rois)
    corr_to_retain=np.where(fit>=0.75)[0]
    np.save('pvals_pearson_whole_brain.npy',pearson_pvals[corr_to_retain,:,0].flatten())
    np.save('fit_null_empirical_whole_brain.npy',fit)
    rej,corr_pvals=fdrcorrection(pearson_pvals[corr_to_retain,:,0].flatten(),alpha=0.05,method='p')
    
if __name__ == "__main__":
    main()     
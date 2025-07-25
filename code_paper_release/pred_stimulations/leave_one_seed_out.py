#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 11:22:30 2022

@author: ludovicocoletta
"""

import glob
import numpy as np
import os
from utils.utils import flag_duplicated_seeds_based_file_name
from utils.utils import remove_outlier_from_seeds_maps
from utils.utils import leave_one_seed_out
import re
import nibabel as nib
import pandas as pd

def main():
    
    folders_of_int=sorted(glob.glob('01_SBA_SPEECH_ARREST*'))
    #folders_of_int=['01_SBA_AMODAL_ANOMIA','01_SBA_SPEECH_ARREST','01_SBA_SENSORIAL']
    
    #thr_dice=0.25
    
    filename_mni='/home/ludovicocoletta/Documents/IntraOpMap_2022/07_rsfmri_APSS/IntraOpMap_RestingState/templates/GM_cort_subcort_mask_2mm.nii.gz'
    img = nib.load(filename_mni) 
    aff = img.affine
    head = img.header
    
    os.makedirs('maps_2mm_cleaned',exist_ok=True)

    for ii in folders_of_int:
        
        os.chdir(ii)
        
        seed_maps_positive=sorted(glob.glob('group_level_randomise/'+'*bin.nii.gz'))
        
        unique_indices=flag_duplicated_seeds_based_file_name(seed_maps_positive)
        
        seed_maps_no_duplicates=[seed_maps_positive[ii] for ii in unique_indices]
                          
        seed_maps_no_duplicates_no_out=seed_maps_no_duplicates.copy()       
        filenames=[iii.split('/')[-1] for iii in seed_maps_no_duplicates_no_out]
        all_coords=[re.findall("-?[0-9]+_-?[0-9]+_-?[0-9]+", iii)[0] for iii in filenames]
        
        path_to_stim_points=[glob.glob('seeds_2mm/'+'*'+iii+'*2mm*WM*gz')[0] for iii in all_coords]
        seed_maps_negative=[glob.glob('group_level_randomise_inverted/'+'*'+iii+'*bin*gz')[0] for iii in all_coords]
        
        accuracy=leave_one_seed_out(seed_maps_no_duplicates_no_out,seed_maps_negative,path_to_stim_points)
        
        os.chdir('..')
        
        pd.DataFrame(accuracy).to_csv('maps_2mm_cleaned/'+'_'.join(ii.split('_')[2:])+'.csv')
        
if __name__ == "__main__":
    main()  

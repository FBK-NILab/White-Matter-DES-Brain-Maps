#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 10:41:11 2023

@author: ludovicocoletta
"""

import nibabel as nib 
import numpy as np
import pandas as pd
import glob


def main():

    # SUBJ OF INTEREST
    df=pd.read_csv("participants.tsv", sep="\t")
    df_no_nan=df.dropna()
    df_final=df_no_nan[df_no_nan['aphasia_type_from_wab'] != 8] # no aphasia diagnosis
    
    unique_subjects=df_final['subject'].unique().tolist()
    
    # Relationship between DES frequency value and symptoms' severity during chronic phase (at least 6 months after stroke)
    # We will probably need to control for age at scan as well, and stroke's onset time
    # We take the first valid time point
    list_analysis_1=np.zeros((len(unique_subjects),19),dtype='object')
      
    # we fill list 1
    
    for sub_index, sub in enumerate(unique_subjects):
        
        if sub=="'M2036'":
            continue
        
        df_sel=df_final[df_final['subject']==sub]
        df_sel=df_sel.sort_values(by=['Days_POS_MRI'])
        time_between_scan_beh=df_sel['MRI_daysPOS_minus_WAB_daysPOS'].to_numpy().astype(float)
        days_post_mri=df_sel['Days_POS_MRI'].to_numpy().astype(float)
        ind_of_int=np.where((time_between_scan_beh<=10)&(days_post_mri>=180))[0]
        
        if ind_of_int.shape[0]>=1:
            
            list_analysis_1[sub_index,0]=df_sel['participant_id'].iloc[ind_of_int[0]]
            list_analysis_1[sub_index,2]=float(df_sel['Days_POS_MRI'].iloc[ind_of_int[0]])
            list_analysis_1[sub_index,3]=float(df_sel['MRI_daysPOS_minus_WAB_daysPOS'].iloc[ind_of_int[0]])
            list_analysis_1[sub_index,4:19]=(df_sel[['information_content','fluency_rating','comprehension_yes_no',
                                                     'comprehension_auditory_wor','comprehension_sequential_c','object_naming',
                                                     'word_fluency','sentence_completion','responsive_speech',
                                                     'spontaneous_speech_rating','comprehension_subscore','repetition_subscore',
                                                     'naming_subscore','wab_r_aq','aphasia_type_from_wab']].iloc[ind_of_int[0]]).tolist()
            
            if df_sel.iloc[0,2] != '#NUM!':
                list_analysis_1[sub_index,1]=float(df_sel['age_at_stroke'].iloc[ind_of_int[0]])
     
    indices_to_remove=np.unique(np.where(list_analysis_1[:,1]==0)[0])   
    list_1_cleaned=np.delete(list_analysis_1,indices_to_remove,0)
    
    df_for_stat=df[df['participant_id'].isin(list(list_1_cleaned[:,0]))]
    age_paper_stat=df_for_stat['age_at_stroke'].to_numpy().astype(float)
    age_paper_stat.mean()
    age_paper_stat.std(ddof=1)
           
    # NETWORKS

    func_of_int=['SEMANTIC', 'PHONOLOGICAL', 'MOVEMENT_ARREST', 'SPEECH_ARREST', 'TALOZZI_BOSTON']
        
    rel_vol_intersection_hubs=np.zeros((len(list_1_cleaned),len(func_of_int)+1))
    
    for sub_index, sub in enumerate(list_1_cleaned):
        
        print(sub_index)
        sub_id=sub[0].split('_')[0]+'_ses-'+sub[0].split('_')[1]
        lesion_mask_data=nib.load(glob.glob('lesion_MNI_1mm/'+sub_id+'_lesion_mask_'+'*flirt*.nii.gz')[0]).get_fdata()
        rel_vol_intersection_hubs[sub_index,-1]=np.count_nonzero(lesion_mask_data)/1000.

        for func_index,func in enumerate(func_of_int):
            
            if func=='TALOZZI_BOSTON':
                talozzi_normative_map=nib.load('comparison_talozzi/boston_raw_bin_up.nii.gz').get_fdata()
                inters=talozzi_normative_map*lesion_mask_data
                rel_vol_intersection_hubs[sub_index,func_index]=(np.count_nonzero(inters)/np.count_nonzero(talozzi_normative_map))*100
            else:
                func_clusters=nib.load(glob.glob('lesion_MNI_1mm/'+sub_id+'_lesion_mask_'+func+'_cluster_size.nii.gz')[0]).get_fdata()
                func_clusters[func_clusters!=0]=1
                
                func_hubs=nib.load('wm_maps/'+func+'_union_norm_thr_bin.nii.gz').get_fdata()
    
                
                dummy_mul=func_hubs*func_clusters
                
                if np.unique(dummy_mul).shape[0]==1:
                    continue
                else:
                    rel_vol_intersection_hubs[sub_index,func_index]=(np.count_nonzero(func_clusters)/np.count_nonzero(func_hubs))*100
    
    df_to_export=pd.DataFrame(np.hstack(
        (rel_vol_intersection_hubs,
         list_1_cleaned[:,-2].reshape((-1,1)),
         list_1_cleaned[:,-3].reshape((-1,1)),
         list_1_cleaned[:,-4].reshape((-1,1))
        )),columns=['SEM_hubs_coverage',
                    'PHONO_hubs_coverage',
                    'MOV_ARREST_hubs_coverage',
                    'SPEECH_ARREST_hubs_coverage',
                    'TALOZZI_hubs_coverage',
                    'tot_lesion_volume',
                    'wab_r_aq',
                    'naming_subscore',
                    'repetition_subscore']
                    )
    df_to_export.to_csv('intersection_with_hubs.csv',index=False)
       
if __name__ == "__main__":
    main() 
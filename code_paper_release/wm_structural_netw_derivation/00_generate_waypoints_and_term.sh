#!/bin/bash

path_to_wm=/home/ludovicocoletta/Documents/REMAP_subcortical/WM_FC/04_LNM_waypoints/seg/WM_re.nii.gz
path_to_cortical_thr=/home/ludovicocoletta/Documents/IntraOpMap_2022/04_preprocessed_rsfmri/00_draft/maps_2mm_cleaned_and_thr/leave_one_seed_out_results/report
path_to_subcortical_thr=/home/ludovicocoletta/Documents/REMAP_subcortical/WM_FC/02_ALL_LNM_MAPS/maps_2mm_cleaned/report
path_to_cortical_nets=/home/ludovicocoletta/Documents/IntraOpMap_2022/04_preprocessed_rsfmri/00_draft/maps_2mm_cleaned_and_thr/netw_no_thr
path_to_subcortical_nets=/home/ludovicocoletta/Documents/REMAP_subcortical/WM_FC/02_ALL_LNM_MAPS/with_repeated_elements
path_5htt=/home/ludovicocoletta/Documents/REMAP_subcortical/WM_FC/04_LNM_waypoints/seg/5htt_MNI.nii.gz

mkdir -p way_points_terminations

for net in SEMANTIC PHONOLOGICAL SPEECH_ARREST

do

  path_to_cortical_map_of_int=$(echo $path_to_cortical_nets/${net}_POSITIVE*gz)
  cortical_thr_of_int=$(cat $path_to_cortical_thr/${net}_hub_thr.txt)
  
  path_to_subcortical_map_of_int=$(echo $path_to_subcortical_nets/${net}_union_randomise.nii.gz)
  subcortical_thr_of_int=$(cat $path_to_subcortical_thr/${net}_hub_thr.txt)
  
  #thr cortical maps, upsample to 1mm, and create gm-wm interface mask
  
  fslmaths $path_to_cortical_map_of_int -thr $cortical_thr_of_int -bin ${net}_cortical_mask.nii.gz
  
  flirt -in ${net}_cortical_mask.nii.gz -ref ${net}_cortical_mask.nii.gz -applyisoxfm 1.0 -interp nearestneighbour -nosearch -out ${net}_cortical_mask.nii.gz
  
  5tt2gmwmi \
    -mask_in ${net}_cortical_mask.nii.gz \
    $path_5htt \
    $PWD/way_points_terminations/${net}_cortical_mask_gw_wm.nii.gz
    
  fslmaths $PWD/way_points_terminations/${net}_cortical_mask_gw_wm.nii.gz -bin $PWD/way_points_terminations/${net}_cortical_mask_gw_wm.nii.gz
    
  #thr subcortical maps, upsample to 1mm, and and mask
  
  fslmaths $path_to_subcortical_map_of_int -thr $subcortical_thr_of_int -bin ${net}_subcortical_mask.nii.gz
  flirt -in ${net}_subcortical_mask.nii.gz -ref ${net}_subcortical_mask.nii.gz -applyisoxfm 1.0 -interp nearestneighbour -nosearch -out ${net}_subcortical_mask.nii.gz
  fslmaths ${net}_subcortical_mask.nii.gz -mul $path_to_wm $PWD/way_points_terminations/${net}_subcortical_mask.nii.gz
  
  rm $PWD/*mask.nii.gz
  
done

  
    

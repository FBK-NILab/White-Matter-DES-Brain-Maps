#!/bin/bash

path_tck=$PWD/HCP_tck #edit this
echo $path_tck/*/*/*/track.tck | tr " " "\n" > tcklist.txt
path_to_mask=$FSLDIR/data/standard/MNI152_T1_1mm_brain_mask.nii.gz

for path_to_tck_file in $(cat tcklist.txt)
do  
    dir_name=$(dirname $path_to_tck_file)
    
    for func in ANOMIA SEMANTIC PHONOLOGICAL SPEECH_ARREST
    do
    

        mkdir -p ${dir_name}/wm_hubs_${func}
             
        time tckedit \
            -include $PWD/way_points_terminations/${func}_cortical_mask_gw_wm.nii.gz \
            -include $PWD/way_points_terminations/${func}_subcortical_mask.nii.gz \
            $path_to_tck_file \
            ${dir_name}/tck_filt.tck
                   
        time tckmap \
            -force \
            -precise \
            -template $FSLDIR/data/standard/MNI152_T1_1mm_brain.nii.gz \
            ${dir_name}/tck_filt.tck \
            ${dir_name}/wm_hubs_${func}/${func}_wh_hubs.nii.gz
            
        min_max_val_single_image=( $(fslstats ${dir_name}/wm_hubs_${func}/${func}_wh_hubs.nii.gz -k $path_to_mask -R) )
        val_range_single_image=$(echo "scale=32;${min_max_val_single_image[1]}-${min_max_val_single_image[0]}" | bc)
        fslmaths -dt double ${dir_name}/wm_hubs_${func}/${func}_wh_hubs.nii.gz -sub ${min_max_val_single_image[0]} -div $val_range_single_image ${dir_name}/wm_hubs_${func}/${func}_wh_hubs_norm.nii.gz
                        
        rm ${dir_name}/tck_filt.tck
        rm ${dir_name}/wm_hubs_${func}/${func}_wh_hubs.nii.gz


    done

    
done


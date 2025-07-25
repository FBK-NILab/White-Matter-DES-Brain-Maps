#!/bin/bash

path_tracto=$PWD/HCP_tck #edit this
path_to_mask=$FSLDIR/data/standard/MNI152_T1_1mm_brain_mask.nii.gz

for func in SEMANTIC PHONOLOGICAL ANOMIA AMODAL_ANOMIA SPEECH_ARREST VERBAL_APRAXIA MOVEMENT_ARREST

do
    echo $path_tracto/*/*/*/wm_hubs_${func}/*gz | tr " " "\n" > ${func}_hubs.txt
    
    fslmerge -t ${func}_as_ts.nii.gz $(cat ${func}_hubs.txt)
    fslmaths ${func}_as_ts.nii.gz -Tmean ${func}_mean.nii.gz
    
    rm ${func}_as_ts.nii.gz
    
    struct=(`cat "${func}_hubs.txt"`)
    for (( i = 0 ; i < ${#struct[@]} ; i++))
    do
        printf "%s " "${struct[$i]} -add"          
    done > ${func}_file_list_all_images.txt

    fslmaths -dt double $(cat ${func}_file_list_all_images.txt) 0 ${func}_union.nii.gz -odt double

    min_max_val_single_image=( $(fslstats ${func}_union.nii.gz -k $path_to_mask -R) )
    val_range_single_image=$(echo "scale=32;${min_max_val_single_image[1]}-${min_max_val_single_image[0]}" | bc)
    fslmaths -dt double ${func}_union.nii.gz -sub ${min_max_val_single_image[0]} -div $val_range_single_image ${func}_union_norm.nii.gz
     
done

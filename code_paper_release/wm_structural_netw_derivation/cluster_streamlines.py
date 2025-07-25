#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 12:20:13 2023

@author: ludovicocoletta
"""

from dipy.io.streamline import load_tck, save_tck
from dipy.io.stateful_tractogram import Space, StatefulTractogram
from dipy.segment.clustering import QuickBundles
import time
import glob
from multiprocessing import Pool

def cluster_tck(in_tck):
    start_time=time.time()
    data=load_tck(in_tck, reference, bbox_valid_check=True)
    streamlines = data.streamlines
    qb = QuickBundles(threshold=20.)
    clusters = qb.cluster(streamlines)
    out_name=in_tck.replace('tck_filt','tck_clust')
    trac=StatefulTractogram(clusters.centroids,reference=reference,space=Space.RASMM)
    save_tck(trac,out_name)
    return time.time() -start_time
    

def log_result(result):
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    times.append(result)
    
def main():
    
    global reference    
    #reference='/home/ludovico/Tools/fsl/data/standard/MNI152_T1_1mm.nii.gz'
    reference='/home/ludovicocoletta/Tools/fsl/data/standard/MNI152_T1_1mm.nii.gz'
    global times
    tcks=sorted(glob.glob('HCP_tck/*/*/*/*/*_tck_filt.tck'))
    pool = Pool(32)
    times=pool.map(cluster_tck, tcks)
    pool.close()
    

if __name__ == "__main__":
        main()  

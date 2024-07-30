# Author: Antti Kiviaho
# Date: 14.2.2024
# 
# Updated 1.3.2024
# Concatenating Visium sections into a pseudobulk for later testing

import os 
os.chdir('/lustre/scratch/kiviaho/prostate_spatial/')

import scanpy as sc
import numpy as np
import pandas as pd
from tqdm import tqdm

from utils import load_from_pickle,get_sample_ids_reorder


samples = get_sample_ids_reorder()
# Create an object to save to
pseudobulk_df = pd.DataFrame()

for s in tqdm(samples,desc='Processing samples'):

    # Copy the slide
    slide = sc.read_h5ad('./data/visium_with_regions/'+s+'_with_regions.h5ad')

    # Calculate a sum over the spots to get to a pseudobulk
    s_pseudobulk_df = pd.DataFrame(slide.layers['counts'].sum(axis=0).T,index=slide.var.index,columns=[s])
    
    # Merge to create a dataframe with samples as columns, genes as index
    pseudobulk_df = pd.merge(pseudobulk_df,s_pseudobulk_df,left_index=True,right_index=True,how='outer')

# Fill empty
pseudobulk_df = pseudobulk_df.fillna(0)

pseudobulk_df.to_csv('data/spatial_pseudobulk_unnormalized.csv')

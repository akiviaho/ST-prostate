import os 
os.chdir('/lustre/scratch/kiviaho/prostate_spatial/')

import scanpy as sc
import numpy as np
import pandas as pd

from tqdm import tqdm
from utils import get_sample_ids_reorder, save_to_pickle
samples = get_sample_ids_reorder()

import warnings
warnings.filterwarnings("ignore")

if __name__ == '__main__':

######### Scanpy scoring method ###########

    gene_set_df = pd.read_excel('gene_sets_of_interest.xlsx',header=None).drop(columns=0).set_index(1).T
    custom_gene_set_scanpy_scores = {}

    # Calculate the scanpy scores for custom gene sets
    for sample in tqdm(samples, desc="Processing sample", unit="sample"):
        slide = sc.read_h5ad('./data/visium_with_regions/'+sample+'_with_regions.h5ad')

        for col in gene_set_df.columns:
            sc.tl.score_genes(slide,gene_set_df[col].dropna(),score_name=col,random_state=2531035)

        custom_gene_set_scanpy_scores[sample] = slide.obs[gene_set_df.columns].copy()

        # Save the dict object
        save_to_pickle(custom_gene_set_scanpy_scores,'./data/spatial_scanpy_score_results.pkl')

    print('scoring done!')
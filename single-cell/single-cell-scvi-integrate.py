from utils import load_from_pickle, save_to_pickle
import scanpy as sc
import anndata as ad
import scib
import scvi
import seaborn as sns
import torch
import os
from matplotlib import pyplot as plt
from datetime import datetime

# Author: Antti Kiviaho
# Date: 27.2.2023
#
# A script for running scvi integration on single cell datasets
# Dong 2020, Chen 2021, Cheng 2022, Chen 2022, Song 2022, Wong 2022, Hirz 2023

if __name__ == "__main__":

    current_date = datetime.today().strftime('%Y%m%d')

    ###############

    # Load the data
    adata = load_from_pickle('./sc-reference/normalized_sc_7_datasets.pickle')
    adata = ad.concat(adata)
    adata.obs_names_make_unique() # Some duplicate indices persists

    # Preprocess and scale
    adata.obs.dataset = adata.obs.dataset.astype('category')
    scib.preprocessing.scale_batch(adata,batch='dataset')
    print('Scaling done...')

    adata.raw = adata

    adata = scib.preprocessing.hvg_batch(adata,batch_key='dataset',target_genes=2000,flavor='seurat',adataOut=True)
    print('HVGs calculated...')

    print('CUDA is available: ' + str(torch.cuda.is_available()))
    print('GPUs available: ' + str(torch.cuda.device_count()))

    print('Initiating training on GPU ...')
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="dataset")
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")

    vae.train(use_gpu=True)

    vae.save('scvi_model_'+current_date)
    adata.obsm["X_scVI"] = vae.get_latent_representation()

    sc.pp.neighbors(adata, use_rep="X_scVI",random_state=745634)
    sc.tl.umap(adata,random_state=745634)
    sc.tl.leiden(adata, key_added="VI_clusters")
    print('NN graph, UMAP & Leiden ready...')

    save_to_pickle(adata,'scvi_integrated_7_sc_datasets_'+current_date+'.pkl')
    print('SCVI integration saved...')

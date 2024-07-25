import os
os.chdir('/lustre/scratch/kiviaho/prostate_spatial')

import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=float32,force_device=True'
import cell2location

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

from utils import load_from_pickle

results_folder = './c2l-results/'
date = '20240125'
# create paths and names to results folders for reference regression and cell2location models
run_name = results_folder + 'cell2location_map_'+ date + '/'



if __name__ == '__main__':

    # Load the single-cell cell type reference: export estimated expression in 'cell type'
    # cell2location_map_20230908/reference_signatures is the correct one
    adata_ref = sc.read_h5ad('c2l-results/cell2location_map_20230908/reference_signatures/sc_reference_signatures.h5ad')

    if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
        inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_ref.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_ref.uns['mod']['factor_names']]].copy()
    inf_aver.columns = adata_ref.uns['mod']['factor_names']
    
    del adata_ref


    # Load visium data and set it up properly (raw, unnormalized counts)
    adata_vis_individually = load_from_pickle('./data/normalized_no_pathology_filter_visium_data.pkl')

    # ADDED 25.1.2024 – concatenate also the available ARNEO samples
    adata_vis_arneo = load_from_pickle('./arneo/data/normalized_arneo_visium_data.pkl')
    adata_vis_individually.update(adata_vis_arneo)

    adata_vis = sc.concat(adata_vis_individually)
    del adata_vis_individually
    
    # find shared genes and subset both anndata and reference signatures
    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    # prepare anndata for cell2location model
    cell2location.models.Cell2location.setup_anndata(adata=adata_vis,
                                                      batch_key='sample_id',
                                                      layer='counts')

    print('adata_vis prior to training:')
    adata_vis

    # create and train the model
    mod = cell2location.models.Cell2location(
        adata_vis, cell_state_df=inf_aver,
        # the expected average cell abundance: tissue-dependent
        # hyper-prior which can be estimated from paired histology:
        N_cells_per_location=21,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection:
        detection_alpha=20
    )
    mod.view_anndata_setup()

    # Train the model
    mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=30851, # this number is for Tampere + ARNEO ADT/APA 4 batches, three batches Tampere only 34000
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True,
         )

    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_vis = mod.export_posterior(
        adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
    )

    # Save anndata object with results
    adata_file = 'visium_adata_with_c2l_mapping_'+date+'.h5ad'
    adata_vis.write(adata_file)
    adata_file

    # Save model
    mod.save(run_name, overwrite=True)

    ## PLOTTING ##
    # plot ELBO loss history during training, removing first 100 epochs from the plot
    mod.plot_history(1000)
    plt.legend(labels=['full data training'])
    plt.savefig('c2l_mapping_model_training_ELBO_'+date+'.png')
    plt.clf()


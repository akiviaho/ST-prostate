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


if __name__ == '__main__':
    
    #########

    # Change model setup in terms of layer & labels !

    # Change these paths when re-running
    results_folder = './c2l-results/cell2location_map_20230908'
    adata_ref_path = './single_cell_reference_with_nmf_derived_annotations_20230908.h5ad'
    ########

    # create paths and names to results folders for reference regression and cell2location models
    ref_run_name = results_folder + 'reference_signatures/'


    # Check if paths exist, and create them if not
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
    if not os.path.exists(ref_run_name):
        os.makedirs(ref_run_name)

    adata_ref = sc.read_h5ad(adata_ref_path)  

    # prepare anndata for the regression model
    cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        #    layer='counts',
                            # 10X reaction / sample / batch
                            batch_key='dataset',
                            # cell type, covariate used for constructing signatures
                            labels_key='final_annotation',
                            categorical_covariate_keys=['sample']
                        )


    # create the regression model
    from cell2location.models import RegressionModel
    mod = RegressionModel(adata_ref)

    # view anndata_setup as a sanity check
    mod.view_anndata_setup()
    mod.train(max_epochs=250, use_gpu=True)

    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_ref = mod.export_posterior(
        adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
    )

        # Save model
    mod.save(ref_run_name, overwrite=True)

    # Save anndata object with results
    adata_file = ref_run_name + 'sc_reference_signatures.h5ad'
    adata_ref.write(adata_file)
    adata_file

    mod.plot_history(20)
    plt.savefig(results_folder+'c2l_single_cell_reference_training_history.png')
    plt.clf()

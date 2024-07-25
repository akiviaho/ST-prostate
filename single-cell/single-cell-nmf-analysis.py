import os 
os.chdir('/lustre/scratch/kiviaho/prostate_spatial')
import numpy as np
import pandas as pd
import scanpy as sc
from itertools import product
import nimfa
import argparse


def save_to_pickle(obj,filename):
    import pickle
    with open(filename, 'wb') as handle:
        pickle.dump(obj, handle, protocol=pickle.HIGHEST_PROTOCOL)

def nmf_preprocessing(dat,n_var_genes=2000):

    dat.X = dat.layers['counts'].copy()
    sc.pp.filter_genes(dat,min_counts=10)
    sc.pp.normalize_total(dat)
    sc.pp.scale(dat)
    sc.pp.highly_variable_genes(dat,n_top_genes=n_var_genes,flavor='seurat_v3',
                                subset=True, layer='counts')

    # Replace negative entries with zero
    dat.X[dat.X < 0] = 0

    return(dat)

def calculate_matrix_orders(arr):
    # get the indices that would sort each row of the array
    sort_indices = np.argsort(-arr, axis=1)

    # create an array to mark the sorted order
    sorted_order = np.empty_like(sort_indices)
    rows, cols = np.indices(arr.shape)
    sorted_order[rows, sort_indices] = cols

    # replace each entry in the original array with its index in the sorted order
    result = sorted_order.astype(int)

    return result


# Helper function to calculate the Jaccard index of two lists
def jaccard(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    intersection = set1.intersection(set2)
    union = set1.union(set2)
    return len(intersection)/len(union)


################# Params ########################

parser = argparse.ArgumentParser()
parser.add_argument('--filename', type=str, help='Name of the file to be used')
parser.add_argument('--n_var_genes', type=int,default=2000, help='The target number of HVGs')
parser.add_argument('--min_genes_in_initial_module', type=int,default=5, help='The number of genes required in the initial NMF module')
parser.add_argument('--required_jaccard_index', type=float,default=0.05, help='The required Jaccard overlap for two gene modules')
args = parser.parse_args()

filename = args.filename
n_hvgs = args.n_var_genes
min_genes_in_initial_module = args.min_genes_in_initial_module
required_jaccard_index = args.required_jaccard_index
min_cells = 100
nmf_comp_range = np.arange(5,11)[::-1] # Changed in the latest run

if __name__ == '__main__':

    # Download data
    ctype_dat = sc.read_h5ad(filename)

    # Initialize dict objects
    nmf_sample_dict = {}
    dq_samples = {}

    ################# First part: calculating NMFs ########################

    for s in np.unique(ctype_dat.obs['sample']):

        if len(ctype_dat.obs[ctype_dat.obs['sample']==s]) >= min_cells:

            # Subset the data to only include a single sample
            dat = ctype_dat[ctype_dat.obs['sample'] == s]
            
            # Do a sample-specific preprocessing
            dat = nmf_preprocessing(dat,n_hvgs)

            for n_comps in nmf_comp_range:

                # Perform non-negative matrix factorization using nsNMF
                f = nimfa.Nsnmf(dat.X, rank=n_comps)
                f_fit = f()
                    
                # Extract the resulting matrices into variables W and H
                W = np.array(f_fit.basis())
                H = np.array(f_fit.coef())

                res = {'samples':dat.obs,'genes':dat.var,'sample_weights':W,'gene_weights':H}

                gene_w = res['gene_weights'].T
                genes = list(res['genes'].index)

                mat_1 = calculate_matrix_orders(gene_w.T).T
                mat_2 = calculate_matrix_orders(gene_w)

                # Iterate through the factors (columns)
                genes_by_factors = {}
                for factor in range(gene_w.shape[1]):
                    valid_genes = list()
                    # Iterate through the genes, starting from the highest weighted gene of this factor
                    for i in np.argsort(mat_1[:,factor]): 
                        if mat_2[i,factor] ==0: # Is this the factor the gene effects the most?
                            
                            # If yes, add it to the list of genes
                            valid_genes.append(genes[i])
                        else:
                            # if not, stop adding genes into the list, move on to the next factor
                            break 
                    # Append the list of valid genes into a dictionary under the appropriate key
                    genes_by_factors['factor'+str(factor)] = valid_genes

                
                all_lists_have_at_least = True

                for lst in genes_by_factors.values():
                    if len(lst) < min_genes_in_initial_module:
                        all_lists_have_at_least = False
                        break

                if all_lists_have_at_least:
                    print(s+": valid factors (min. "+str(min_genes_in_initial_module)+" genes) when using n_comps="+str(n_comps))
                    nmf_sample_dict[s] = genes_by_factors
                    break
        else:
            dq_samples[s] = len(ctype_dat.obs[ctype_dat.obs['sample']==s])
        
    for k in list(dq_samples.keys()):
        print(k + ' not processed for too few cells (' + str(dq_samples[k])+')')

    ################# Second part: Observing module overlaps ########################

    # A dictionary to keep track of overlapping factors
    overlapping_factors = {}

    # Go through each sample factor combination
    for (sample1, factor1), (sample2, factor2) in product([(sample, factor) for sample in nmf_sample_dict for factor in nmf_sample_dict[sample]], repeat=2):
        # Ignore if it's the same sample and factor combination
        if sample1 == sample2 and factor1 == factor2:
            continue
        # Calculate the Jaccard index of the two lists
        jaccard_index = jaccard(nmf_sample_dict[sample1][factor1], nmf_sample_dict[sample2][factor2])
        # If the overlap is more than 5%, add it to the overlapping_factors dictionary
        if jaccard_index > required_jaccard_index:
            overlapping_factors.setdefault((sample1, factor1), set()).add((sample2, factor2))
            
    # A set to keep track of factors with more than two overlaps
    valid_factors = set()

    # Go through each factor
    for key, value in overlapping_factors.items():
        # If the factor overlaps with at least one other factors, add it to the valid_factors set
        if len(value) >= 2:
            valid_factors.add(key)

    # A dictionary to keep track of how many overlapping factors each gene has
    gene_overlaps = {}

    # Go through each valid factor
    for factor in valid_factors:
        # Get the sample and factor
        sample, factor_name = factor
        # Go through each gene in the factor
        for gene in nmf_sample_dict[sample][factor_name]:
            # Add one to the gene's overlaps count
            gene_overlaps[gene] = gene_overlaps.get(gene, 0) + 1
            
    # Create an empty 2D array to serve as the adjacency matrix
    adj_matrix = [[0 for gene2 in gene_overlaps.keys()] for gene1 in gene_overlaps.keys()]

    # Loop through each valid factor
    for factor in valid_factors:
        # Get the sample and factor
        sample, factor_name = factor
        
        # Loop through each gene in the factor
        for gene1 in nmf_sample_dict[sample][factor_name]:
            # Loop through each other gene in the same factor
            for gene2 in nmf_sample_dict[sample][factor_name]:
                # Ignore if it's the same gene
                if gene1 == gene2:
                    continue
                # Increment the value for this pair of genes in the adjacency matrix
                adj_matrix[list(gene_overlaps.keys()).index(gene1)][list(gene_overlaps.keys()).index(gene2)] += 1

    adj_df = pd.DataFrame(adj_matrix,columns=list(gene_overlaps.keys()),index=list(gene_overlaps.keys()))

    # Save the adjacency matrix to a csv
    save_file = filename.replace('.h5ad','_nmf_derived_gene_adjacencies_mod_'+
                                 str(min_genes_in_initial_module)+'_jac'
                                 +str(required_jaccard_index)+'.csv')
    adj_df.to_csv(save_file)



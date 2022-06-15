#Run this script in an interactive session on the HPC. Use the conda env: /gpfs/data/fisherlab/conda_envs/scVelo
#Requires PHATE and PAGA initialized umap coordinates in 'umap' file
#'umap' denotes seurat computed umap. 'paga_umap' denotes paga intitialized umap coordianates (computed via scanpy_functions.py). 'phate' denotes PHATE 
#coordiantes (computed via scanpy_functions.py)

#Much of this script still needs to be parameterized. e.g. cellranger directories

import anndata
import pandas as pd
import numpy as np
import matplotlib as plt
import argparse
import math

def main(args):

	integrate_looms(args.sample_obs_file, args.emb_file, args.clus_map_file, args.out_dir, args.sample_ids_file)

def integrate_looms(sample_obs_file, emb_file, clus_map_file, out_dir, sample_ids_file):
	
	sample_obs = pd.read_csv(sample_obs_file)
	embeddings = pd.read_csv(emb_file)
	seurat_clusters = pd.read_csv(clus_map_file)
	samples = pd.read_csv(sample_ids_file)

	sample_list = samples.iloc[:,0].tolist()

	looms = []
	for i in sample_list:
		looms.append(load_looms(i, sample_obs))
		integrated_loom = looms[0].concatenate(looms[1:], index_unique = None)

	#Add dimmensionality reduction coordinates
	loom_index = pd.DataFrame(integrated_loom.obs.index)
	loom_index = loom_index.rename(columns = {0:'CellID'})
	embeddings = embeddings.rename(columns = {'Unnamed: 0':'CellID'})
	embeddings_ordered = loom_index.merge(embeddings, on = "CellID")
	
	n_emb = len(embeddings_ordered.columns)/2
	n_emb = math.floor(n_emb)
	ind = 0
	for i in range(1,n_emb+1):
		integrated_loom.obsm['_'.join(embeddings_ordered.columns[i+ind].split("_")[:-1])] = embeddings_ordered.iloc[:,i+ind:i+2+ind].values
		ind += 1

	#Add clusters
	seurat_clusters = seurat_clusters.rename(columns = {'Unnamed: 0':'CellID'})
	seurat_clusters_ordered = loom_index.merge(seurat_clusters, on = "CellID")
	seurat_clusters_ordered = seurat_clusters_ordered.iloc[:,1:]
	integrated_loom.uns['seurat_clusters'] = seurat_clusters_ordered.astype('category').values
	integrated_loom.obs['seurat_clusters'] = seurat_clusters_ordered.astype('category').values

	integrated_loom.write_loom(filename='%sintegrated_loom' % (out_dir), write_obsm_varm=True)

def load_looms(sample, sample_obs):
	cellID = sample_obs[sample_obs['x'].str.contains("%s:" % (sample))]
	loom = anndata.read_loom("/gpfs/data/giannarellilab/Mike_G/scRNAseq_52021/cellranger_rerun/%s/%s.loom" % (sample, sample))
	loom = loom[np.isin(loom.obs.index, cellID)]
	loom.var_names_make_unique()
	return(loom)

def parseArguments():
	
	parser = argparse.ArgumentParser(prog="Make_int_loom", description='')
	required = parser.add_argument_group('required arguments')

	required.add_argument('-cell_obs', nargs='?', required=True, help='Cell barcode to sample mappings' , dest='sample_obs_file')
	required.add_argument('-emb', nargs='?', required=True, help='Cell barcode to embedding coordinate mappings. Ive been including 4 different embeddings. This script will add any number of embeddings to the object that are supplied in the embedding file', dest='emb_file')
	required.add_argument('-clusters', nargs='?', required=True, help='Cell barcode to cluster mapping', dest='clus_map_file')
	required.add_argument('-o', '--output_dir', nargs='?', required=True, help='Output directory for integrated loom', dest='out_dir')
	required.add_argument('-sample_ids', nargs='?', required=True, help='ids of samples to integrate. Must mach Cellranger directory names' , dest='sample_ids_file')

	return parser.parse_args()

args = parseArguments()
main(args)
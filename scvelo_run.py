#Run this script in an interactive session on the HPC. Use the conda env: /gpfs/data/fisherlab/conda_envs/scVelo
import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt

sample_obs = pd.read_csv("/gpfs/data/giannarellilab/Mike_G/scRNAseq_82021/NK_cells/Cellcycle_regressed/cellID_obs.csv")
umap = pd.read_csv("/gpfs/data/giannarellilab/Mike_G/scRNAseq_82021/NK_cells/Cellcycle_regressed/cell_embeddings.csv")
clusters = pd.read_csv("/gpfs/data/giannarellilab/Mike_G/scRNAseq_82021/NK_cells/Cellcycle_regressed/clusters.csv")

sample_list = ['CV7209', 'CV7221D', 'CV7221HC', 'CV7223', 'CV7231A', 'CV7286', 'CV7291', 'CV7293', 'CV7318', 'CV7320', 'CV7322', 'CV7323POS', 'CV7327', 'CV7330', 'CV7339', 'CV7341', 'CV7342', 'CV7343', 'CV7344', 'CV7346', 'CV7347', 'CV7292']

def load_looms(sample, sample_obs):
	cellID = sample_obs[sample_obs['x'].str.contains("%s:" % (sample))]
	loom = anndata.read_loom("/gpfs/data/giannarellilab/Mike_G/scRNAseq_52021/cellranger_rerun/%s/%s.loom" % (sample, sample))
	loom = loom[np.isin(loom.obs.index, cellID)]
	loom.var_names_make_unique()
	return(loom)

looms = []
for i in sample_list:
	looms.append(load_looms(i, sample_obs))
integrated_loom = looms[0].concatenate(looms[1:], index_unique = None)

#Add umap coordinates
loom_index = pd.DataFrame(integrated_loom.obs.index)
loom_index = loom_index.rename(columns = {0:'CellID'})
umap = umap.rename(columns = {'Unnamed: 0':'CellID'})
umap_ordered = loom_index.merge(umap, on = "CellID")
umap_ordered = umap_ordered.iloc[:,1:]
integrated_loom.obsm['X_umap'] = umap_ordered.values

#Add clusters
clusters = clusters.rename(columns = {'Unnamed: 0':'CellID'})
clusters_ordered = loom_index.merge(clusters, on = "CellID")
clusters_ordered = clusters_ordered.iloc[:,1:]
integrated_loom.uns['clusters'] = clusters_ordered.astype('category').values
integrated_loom.obs['clusters'] = clusters_ordered.astype('category').values

p = ["#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC"]


#scvelo
scv.set_figure_params()
scv.pp.filter_and_normalize(integrated_loom, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(integrated_loom, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(integrated_loom, n_jobs = 40)
scv.tl.velocity(integrated_loom, mode='dynamical')
scv.tl.velocity_graph(integrated_loom)
scv.pl.velocity_embedding(integrated_loom,arrow_length=3, arrow_size=2, dpi=1000, basis = 'umap', save = 'embedding.png', color = 'clusters', palette = p, legend_loc='right')
scv.pl.velocity_embedding_stream(integrated_loom, basis = 'umap', save = 'stream.png', color = 'clusters', palette = p, min_mass = 2, legend_loc='right')
scv.tl.velocity_pseudotime(integrated_loom)
scv.pl.scatter(integrated_loom, color='velocity_pseudotime', cmap='gnuplot', save = 'pseudotime.png', legend_loc='right')

scv.tl.latent_time(integrated_loom)
scv.pl.scatter(integrated_loom, color='latent_time', color_map='gnuplot', size=80, save = 'latenttime.png', legend_loc='right')

top_genes = integrated_loom.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(integrated_loom, var_names=top_genes, sortby='latent_time', col_color='clusters', n_convolve=100, save = 'latenttime_heatmap.png')

top_genes = integrated_loom.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(integrated_loom, basis=top_genes[:20], ncols=5, frameon=True, save = 'top_genes.png')

scv.pl.velocity(integrated_loom, ['NCAM1', 'FCGR3A', 'SELL', 'GNLY', 'ILR7A', 'GZMB'], ncols=2, save = 'interesting_genes.png')

scv.tl.velocity_confidence(integrated_loom)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(integrated_loom, c=keys, cmap='coolwarm', perc=[5, 95], save = 'length_confidence.png')

scv.tl.velocity_clusters(integrated_loom)
scv.pl.scatter(integrated_loom, color='velocity_clusters', save = 'velocity_clusters.png')

scv.tl.score_genes_cell_cycle(integrated_loom)
scv.pl.scatter(integrated_loom, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95], save = 'cell_cycle_scores.png')

scv.tl.terminal_states(integrated_loom)
scv.pl.scatter(integrated_loom, color=['root_cells', 'end_points'], save = 'infered_terminal_states.png')

scv.tl.rank_dynamical_genes(integrated_loom, groupby='clusters')
df = scv.get_df(integrated_loom, 'rank_dynamical_genes/names')
for cluster in ['Transitional 1', 'CD56 bright', 'transitional 2', 'Activated 1']:
	scv.pl.scatter(integrated_loom, df[cluster][:5], ylabel=cluster, frameon=True, save = 'top5_likelihood_%s.png' % (cluster))
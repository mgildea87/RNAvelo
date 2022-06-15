#Run this script in an interactive session on the HPC. Use the conda env: /gpfs/data/fisherlab/conda_envs/scVelo

import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
import argparse
import cellrank as cr
import os
from scipy import sparse


def main(args):
	path = '%sRNA_velocity/' % (args.out)
	isExist = os.path.exists(path)
	if not isExist:
		os.makedirs(path)

	seurat_cluster_cols = pd.read_csv(args.cluster_colors)
	scc = seurat_cluster_cols['x'].tolist()

	if not args.skip_velocity:
		integrated_loom = anndata.read_loom(args.loom)
		integrated_loom = velocity(integrated_loom, args.seurat)
		integrated_loom.write_h5ad(filename='%sintegrated_loom_velo.h5ad' % (path))
	if args.skip_velocity:
		integrated_loom = anndata.read_h5ad(filename='%sintegrated_loom_velo.h5ad' % (path))
	
	if not args.skip_scVelo_plots:
		scVelo_plots(integrated_loom, path, scc, args.embedding)

	if not args.skip_CellRank:
		integrated_loom = CellRank(integrated_loom)
	
	if not args.skip_CellRank_plots:	
		CellRankPlots(integrated_loom, path, scc, args.embedding)

def velocity(integrated_loom, seurat):
	
	scv.pp.filter_and_normalize(integrated_loom, min_shared_counts=100, n_top_genes=2000)

	seurat = anndata.read_h5ad(seurat)
	seurat = seurat[seurat.obs_names.isin(integrated_loom.obs_names)]
	seurat = seurat[seurat.obs_names.sort_values()]
	integrated_loom = integrated_loom[integrated_loom.obs_names.sort_values()]
	seurat.var_names = seurat.var['SCT_features']
	seurat = seurat[:, seurat.var_names.isin(integrated_loom.var_names)]
	seurat = seurat[:, seurat.var_names.sort_values()]
	integrated_loom = integrated_loom[:,integrated_loom.var_names.sort_values()]
	integrated_loom = integrated_loom[:, integrated_loom.var_names.isin(seurat.var_names)]
	integrated_loom.X = seurat.X

	scv.pp.filter_and_normalize(integrated_loom, min_shared_counts=0, n_top_genes=2000)
	scv.pp.moments(integrated_loom, n_pcs=30, n_neighbors=30)
	scv.tl.recover_dynamics(integrated_loom, n_jobs = 40)
	scv.tl.velocity(integrated_loom, mode='dynamical')
	scv.tl.velocity_graph(integrated_loom)
	scv.tl.latent_time(integrated_loom)
	scv.tl.velocity_pseudotime(integrated_loom)
	scv.tl.velocity_confidence(integrated_loom)
	scv.tl.terminal_states(integrated_loom)
	scv.tl.rank_dynamical_genes(integrated_loom, groupby='seurat_clusters')
	scv.tl.paga(integrated_loom, groups='seurat_clusters')
	return(integrated_loom)

#scvelo
def scVelo_plots(integrated_loom, path, scc, embedding):
	path = '%sscVelo_outs/' % (path)
	isExist = os.path.exists(path)

	if not isExist:
		os.makedirs(path)

	scv.pl.velocity_embedding(integrated_loom,arrow_length=3, arrow_size=2, dpi=1000, basis = embedding, save = '%s%s_embedding.png' % (path, embedding), color = 'seurat_clusters', palette = scc, legend_loc='right')
	scv.pl.velocity_embedding(integrated_loom,arrow_length=3, arrow_size=2, dpi=1000, basis = embedding, save = '%s%s_embedding.png' % (path, embedding), color = 'seurat_clusters', palette = scc, legend_loc='right')

	scv.pl.velocity_embedding_stream(integrated_loom, basis = embedding, save = '%sstream_%s_embedding.png' % (path, embedding), color = 'seurat_clusters', palette = scc, min_mass = 1, legend_loc='right', density = 2, linewidth = 0.5)

	scv.pl.velocity_embedding_grid(integrated_loom, basis = embedding, save = '%s%s_grid.png' % (path, embedding), color = 'seurat_clusters', palette = scc, min_mass = 1, legend_loc='right')

	scv.pl.scatter(integrated_loom, color='velocity_pseudotime', cmap='gnuplot', save = '%spseudotime_%s.png' % (path, embedding), legend_loc='right', basis = embedding)

	scv.pl.scatter(integrated_loom, color='latent_time', color_map='gnuplot', size=80, save = '%slatenttime_%s.png' % (path, embedding), legend_loc='right', basis = embedding)

	top_genes = integrated_loom.var['fit_likelihood'].sort_values(ascending=False).index[:300]
	scv.pl.heatmap(integrated_loom, var_names=top_genes, sortby='latent_time', col_color='seurat_clusters', palette = scc, n_convolve=100, save = '%slatenttime_heatmap.png' % (path))

	top_genes = integrated_loom.var['fit_likelihood'].sort_values(ascending=False).index
	scv.pl.scatter(integrated_loom, basis=top_genes[:20], ncols=5, frameon=True, save = '%stop_genes.png' % (path), palette = scc, color = 'seurat_clusters')

	keys = 'velocity_length', 'velocity_confidence'
	scv.pl.scatter(integrated_loom, c=keys, cmap='coolwarm', perc=[5, 95], save = '%slength_confidence_%s.png' % (path, embedding), basis = embedding)

	df = integrated_loom.obs.groupby('seurat_clusters')[keys].mean().T
	np.savetxt(fname = '%sMean_velocity_confidence_per_cluster.csv' % (path), X = df, delimiter = ',')

	scv.pl.scatter(integrated_loom, color=['root_cells', 'end_points'], save = '%sinfered_terminal_states_%s.png' % (path, embedding), basis = embedding)

	df = scv.get_df(integrated_loom, 'rank_dynamical_genes/names')
	np.savetxt(fname = '%sranked_dynamical_genes_per_cluster.csv' % (path), X = df, delimiter = ',', fmt = '%s')
	
	for cluster in df.columns:
		scv.pl.scatter(integrated_loom, df[cluster][:20], ncols = 5, ylabel=cluster, frameon=True, palette = scc, color = 'seurat_clusters', save = '%stop20_likelihood_%s.png' % (path, cluster))

	scv.pl.velocity_graph(integrated_loom, threshold=.1, color='seurat_clusters', palette = scc, save = '%sConnectivity_graph_%s.png' % (path, embedding), basis = embedding)

	scv.pl.paga(integrated_loom, basis=embedding, size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, save = '%sdirected_PAGA_%s.png' % (path, embedding))
	pd.DataFrame.sparse.from_spmatrix(integrated_loom.uns['paga']['transitions_confidence']).to_csv('%sPAGA_transitions_confidence.csv' % (path))

def CellRank(integrated_loom):
	
	cr.tl.terminal_states(integrated_loom, cluster_key="seurat_clusters")
	cr.tl.initial_states(integrated_loom, cluster_key="seurat_clusters")
	cr.tl.lineages(integrated_loom)
	cr.tl.lineage_drivers(integrated_loom)
	return(integrated_loom)

def CellRankPlots(integrated_loom, path, scc, embedding):
	path = '%sCellRank_outs/' % (path)
	isExist = os.path.exists(path)
	if not isExist:
		os.makedirs(path)
	cr.pl.terminal_states(integrated_loom, palette = scc, discrete=True,save = '%sterminal_states_%s.png' % (path, embedding), basis = embedding)
	cr.pl.initial_states(integrated_loom, palette = scc,discrete=True, save = '%sinitial_states_%s.png' % (path, embedding), basis = embedding)
	cr.pl.terminal_states(integrated_loom, palette = scc, discrete=False,save = '%sterminal_states_bycell_%s.png' % (path, embedding), basis = embedding)
	cr.pl.initial_states(integrated_loom, palette = scc,discrete=False, save = '%sinitial_states_bycell_%s.png' % (path, embedding), basis = embedding)
	cr.pl.lineages(integrated_loom, palette = scc,same_plot=False, save = '%slineages_%s.png' % (path, embedding), basis = embedding)
	cr.pl.lineages(integrated_loom,palette = scc, same_plot=False, save = '%slineages_all_%s.png' % (path, embedding), basis = embedding)
#	cr.pl.cluster_fates(integrated_loom,mode="paga_pie",cluster_key="seurat_clusters",basis="UMAP", palette = scc, title="directed PAGA", save = '%scellfates_agg_PAGA_%s.png' % (path, embedding))
#	cr.pl.cluster_fates(integrated_loom,mode="heatmap",cluster_key="seurat_clusters",basis="UMAP", palette = scc, title="directed PAGA", save = '%scellfates_agg_hm_%s.png' % (path, embedding))
#	cr.pl.lineage_drivers(integrated_loom, lineage="9_Treg Activated", save = '%stop20%s.png' % (path, embedding), n_genes = 20)

def parseArguments():
	
	parser = argparse.ArgumentParser(prog="scvelo_run", description='')
	required = parser.add_argument_group('required arguments')
	velocity = parser.add_argument_group('velocity options')
	CellRank = parser.add_argument_group('CellRank options')
	scVelo_plots = parser.add_argument_group('scVelo plots options')
	CellRankPlots = parser.add_argument_group('CellRank plots options')

	required.add_argument('-l', '--loom', nargs='?', required=False, help='Absolute path to integrated loom file. Not required if using --skip_velocity' , dest='loom')
	required.add_argument('-c', '--cluster_colors', nargs='?', required=False, help='Absolute path to file with cluster color palette. Not required if using bot --skip_scVelo_plots and --skip_CellRank_plots', dest='cluster_colors')
	required.add_argument('-o', '--output_dir', nargs='?', required=True, help='Absolute path to directory where outputs should be written' , dest='out')
	required.add_argument('-e', '--embedding', default='UMAP', required=False, help='Embedding to use for plots. Generally "UMAP". Default is "UMAP"', dest='embedding')
	
	velocity.add_argument('--skip_velocity', action='store_true', required=False, help='If RNA velocity has already been computed for the --loom object via scVelo. Disregard to run scVelo RNA velocity analysis. When velocity is enabled the resulting loom file will be written as integrated_loom_velo. Use this as input along with the skip_velocity argument.', dest='skip_velocity')
	velocity.add_argument('-s', '--seurat', nargs='?', required=False, help='Seurat h5ad file. For adding integrated data to .X in integrated loom. Not required if using --skip_velocity', dest='seurat')

	scVelo_plots.add_argument('--skip_scVelo_plots', action='store_true', required=False, help='Skip scVelo plots', dest='skip_scVelo_plots')

	CellRank.add_argument('--skip_CellRank', action='store_true', required=False, help='Skip CellRank analysis', dest='skip_CellRank')

	CellRankPlots.add_argument('--skip_CellRank_plots', action='store_true', required=False, help='Skip CellRank analysis', dest='skip_CellRank_plots')

	return parser.parse_args()

args = parseArguments()
main(args)
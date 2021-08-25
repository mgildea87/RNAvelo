# RNA velocity analysis workflow
Note: this requires a lot more work. This is just the bare bones I've used on a specific project. The goal is to enable someone else in the group to run RNA velocity analysis with a list of .bam files and a Seurat object. Will get to this when I can. At the moment there's a lot of harcoding in here for a specific project. Additionally, this version is built for integrated datasets e.g. 2 or more scRNA-seq samples that have been integrated. Should enable single samples to be run as well.

The purpose of this work flow is to run an RNA velocity analysis via scVelo (https://scvelo.readthedocs.io/, https://github.com/theislab/scvelo) on scRNA-seq data. 

## Required input files
	
	1. .bam files for each scRNA-seq sample.
		Cellranger automatically outputs sorted and indexed .bam files
	2. Cell barcodes
		Cell ranger automatically outputs this. It is located here: (filtered_feature_bc_matrix/barcodes.tsv.gz) in each Cellranger folder
	3. A processed Seurat object (filtering -> dimmensionality reduction)
		From this object the following will be generated:
			a. A list of cell ID barcodes (this can be any subset of cells)
			b. A table of cell ID barcodes vs embedding coordinates (generally from a umap)
			c. A table of Cell ID barcodes vs cluster or annotation ID's
			d. A table of Cell ID barcodes vs any other meta data that may be useful in the RNA velocity visualizations
			e. A list of color hex codes for coloring cells in the scvelo plots


## Workflow steps:

	1. Generate loom files for each samples from .bam files
		This is accomplished via Velocyto (http://velocyto.org/). This object will contain count matrices for spliced and unspliced reads per gene per cell.
	2. Integrate the loom files and add meta data (umap coordinates, cluster IDs)
	3. Run scVelo dynamical modeling 
		This can of course be changed. There are a lot of parameters in each of the functions
	4. Output a bunch of plots and data?



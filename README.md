# RNA velocity analysis workflow
Note: this requires a lot more work. This is just the bare bones I've used on a specific project. The goal is to enable someone else in the group to run RNA velocity analysis with a list of .bam files and a Seurat object. Will get to this when I can. At the moment there's a lot of harcoding in here for a specific project.

The purpose of this work flow is to run an RNA velocity analysis via scVelo (https://scvelo.readthedocs.io/, https://github.com/theislab/scvelo) on scRNA-seq data. 

## Required input files
	
	1. .bam files for each scRNA-seq sample.
		Cellranger automatically outputs sorted and indexed .bam files
	2. A processed Seurat object
		From this object the following will be generated:
			a. A list of cell ID barcodes
			b. A table of cell ID barcodes vs embedding coordinates (generall from a umap)
			c. A table of Cell ID barcodes vs cluster or annotation ID's
			d. A table of Cell ID barcodes vs any other meta data that may be useful in the RNA velocity visualizations
			e. A list of color hex codes for coloring cells in the scvelo plots


This work flow encompasses the following steps:

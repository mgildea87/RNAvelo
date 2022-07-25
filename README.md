# RNA velocity analysis workflow
Note: this requires some more work. The goal is to enable someone else in the group to run RNA velocity analysis with a list of .bam files and a Seurat object. Will get to this when I can. At the moment there's some harcoding in here for a specific project. Additionally, this version is built for integrated datasets e.g. 2 or more scRNA-seq samples that have been integrated. I should enable single samples to be run as well.

The purpose of this work flow is to run an RNA velocity analysis via scVelo and CellRank (DOI:10.1038/s41587-020-0591-3,https://scvelo.readthedocs.io/, https://github.com/theislab/scvelo, DOI:10.1038/s41592-021-01346-6, https://github.com/theislab/cellrank) on scRNA-seq data. 

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
		Run Velocyto for each sample. This requires .bam files. Velocyto counts reads that overlap exon-exon and intron-exon junctions and outputs a count matrix for spliced and unspliced counts in a format called .loom. 
			a.	Velocyto.sh submits the velocyto jobs contained in velocyto_commands.txt. The parameters are pretty much default. I haven’t played around with it much. Remember where you save each sample’s .loom file. I generally save them in the same folder as that sample’s Cellranger data.
	2.	Next you need to output several things from your R analysis.
		a. R items requried (if using R for analysis)
			i.	Seurat integrated object as .h5ad
        	ii.	Sample:cellbarcodes reformatted to match those in the .loom files
        	iii.	Sample:cellbarcode to cluster mappings. E.g. this cell belongs to cluster X. You can use any Seurat metadata column
        	iv.	Sample:cellbarcode embeddings. This is a table containing each sample:cellbarcode and a 2 columns (x, y) for each embedding. I generally have 4 (PCA, UMAP, PHATE, PAGA initialized UMAP). The script will handle any number of them automatically.
        	v.	Cluster colors. This is a list of color codes from Seurat so the scVelo and CellRank plots use the same palette
		b.	I have included 2 functions in CVRCFunc (the package I’ve been putting useful functions into) that will automatically output all these files given a seurat object as input (and some other parameters). https://github.com/mgildea87/CVRCFunc. You can install the package or copy paste the functions intro your R.
        	i.	ExportSeuratMeta() this outputs ii-v above
        	ii.	ExportSeurath5ad() this outputs the Seurat integrated assay as .h5ad. This function also converts the cellbarcodes to match the .loom barcodes. Eventually we will replace the raw count matrix that is in the .loom file by default with the integrated count matrix.
	3.	Run Make_int_loom.py this script integrates all of the .loom files into a single .loom object. It also adds all of the metadata (ii-iv) from step 2.
		a.	There are several required arguments which are files output by step 2.b.i. If you run python ‘Make_int_loom.py –help’ The required parameters and their flags will be listed. 
		b.	Within the script, the location of the individual sample .loom files is hardcoded. So you will need to modify.
	4.	scvelo_run.py. This script runs scvelo and CellRank and outputs a bunch of plots.
		a.	Importantly, this script requires the Seurat .h5ad integrated assay file output in step 2.b.ii. It will reformat the integrated count matrix to match the integrated_loom and replace the raw count matrix that is in the integrated_loom by default.
		b.	If you run python ‘scvvelo_run.py –help’ The required/optional parameters and their flags will be listed. There are a bunch. I tried to modularize the functions so you don’t have to repeat the computationally heavy steps to make new plots. E.g. see --skip_velocity
		c.	This script is where step 2.v is used to add the same color palette to the plots as in Seurat.

## Notes.
	1.	I made this mistake, but when you do the integration in Seurat you should make it keep all genes in the integrated assay. RNA velocity will pick genes that may not be in your integrated count matrix (top 3000 or whatever). The script will handle that no problem but it will reduce the number of genes being used. RNA velocity picks the top 2000 variable genes (by spliced and unspliced counts I believe). I only included the top 3000 when I did the Seurat integration and when I combined this with the RNA velocity genes, only ~1K overlapped and were used in the CellRank model.
	2.	All of these scripts (except the R functions) have to be run in a conda environment because most of these packages are python packages. You could install all the packages yourself in your python package directory but I would caution against this (as the authors for most of these packages do!). It’s possible there will be conflicts etc. I know everything works in the conda environment. Anyway, to load it: ‘module load miniconda3/cpu/4.9.2’ then ‘conda activate /gpfs/data/fisherlab/conda_envs/scVelo/’. Just don’t install/remove/update anything!
	3.	I run Make_int_loom.py and scvelo_run.py in an interactive session. Or you could write a .sh pipeline for these scripts and submit as a job. Up to you!


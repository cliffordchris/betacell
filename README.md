# betacell README


## Overview: 

The code found in this repo is used to perform QC on raw single cell data obtained from Balboa et al. and to produce a png image.
The image produced describes the 3 variables: (# of RNA features, RNA counts, and mitochondrial gene expression) across the different transcriptomic datasets analyzed.

## Analysis:

The analysis performed is detailed and outlined at: <https://satijalab.org/seurat/archive/v3.0/pbmc3k_tutorial.html>

## Data:

This project repo is used to study publically available transcriptomic data available from:
<https://singlecell.broadinstitute.org/single_cell/study/SCP1526/functional-metabolic-and-transcriptional-maturation-of-human-pancreatic-islets-derived-from-stem-cells?#study-summary>  
This data contains scRNA data from in vitro and in vivo studies of stem cell derived islets and mature human islets.

From the data authors: *"Single-cell transcriptomics of SC-islets in vitro and throughout 6 months of engraftment in mice revealed a continuous maturation trajectory culminating in a transcriptional landscape closely resembling that of primary islets. Our thorough evaluation of SC-islet maturation highlights their advanced degree of functionality and supports their use in further efforts to understand and combat diabetes."*

*The code used to generate our first image is also available in this repo in the "Raw data" folder.*

## Installation:

1) To generate the 'QC plot.png' download the Raw Data in the "Raw data" file, and place it in a RStudio project directory with "main.R".
2) Install the following dependencies:

			R version 4.2.2
		Package:			Version:
		1 devtools			2.4.5
		2 dplyr			1.1.0
		3 remotes			2.4.2
		4 Seurat			4.3.0
		5 SeuratData			0.2.2	   (code to download this library at the top of "main.R", NOT in CRAN repo)
		6 SeuratDisk			0.0.0.9020 (code to download this library at the top of "main.R", NOT in CRAN repo)

3) Run uncommented code, by default the output image will be named "QC plot.png"

## Folder Structure

		<Main>
			.gitattributes				stores lfs file directories (large files)
			.gitignore					stores ignore file directories
			QC plot.png					Output image of QC data for our sample scRNA sequencing data, generated by "main.R"
			main.R 						Main code to process scRNA data, and output "QC plot.png" image.
			<Raw data>
				barcodes.tsv.gz				cell ids
				endocrine.counts.mtx.gz		Endocrine UMI counts
				features.tsv.gz				gene ids

## Citation:

*Functional, metabolic and transcriptional maturation of stem cell derived beta cells.*  
Diego Balboa, Tom Barsby, Väinö Lithovius, Jonna Saarimäki-Vire, Muhmmad Omar-Hmeadi,
Oleg Dyachok, Hossam Montaser, Per-Eric Lund, Mingyu Yang, Hazem Ibrahim, Anna Näätänen,
Vikash Chandra, Helena Vihinen, Eija Jokitalo, Jouni Kvist, Jarkko Ustinov, Anni I. Nieminen,
Emilia Kuuluvainen, Ville Hietakangas, Pekka Katajisto, Joey Lau, Per-Ola Carlsson, Sebastian Barg,
Anders Tengholm, Timo Otonkoski  
bioRxiv 2021.03.31.437748; doi: https://doi.org/10.1101/2021.03.31.437748

*Seurat*. Satijalab https://satijalab.org/seurat/authors.html#citation
# GMMchi_scs_pipeline

GMMchi is the python package for postprocessing single cell RNA-seq data specifically on mixtures of similar cell type. GMMchi_scs_pipeline leverages GMMchi to remove cluster and visulization bias due to highly varying factors such as library size, library complexity as well as highly expressed genes including housekeeping genes, mitochondrial and ribosomal gene. 

Mixtures of similar cell types can be extremely challenging in seperating during clustering and visualization due to the high degree of homology between the expression profile across each cell. In addition, there failure in removing confounding factors may lead to serious misinterpretation in the data. GMMchi_scs_pipeline presents the first step in attempting to achieve a better seperation in mixtures of similar cell types. The package contains the necessary tools for postprocessing, dimensionality-reduction as well as visualization. The output returned by the postprocessing function of this package also provides a simple format if the user wishes to use other downstream visualization or clustering techinques. 

Part of this package is the leverage of GMMchi, another pipeline we built that can be downloaded and found [here](https://github.com/jeffliu6068/GMMchi). 

### Download Package

Download the GMMchi_scs_pipeline package by:
```
pip install git+https://github.com/jeffliu6068/GMMchi_scs_pipeline.git
```
or 
```
pip install GMMchi_scs_pipeline
```

### Import

Once installed, import the package by: 

```
import GMMchi_scs_pipeline
```
## Intuition: How GMMchi_scs_pipeline Works in Postprocessing scRNA-seq Data

Here, we outline each of the step included in the postprocessing pipeline:

1) Remove empty barcodes
2) Remove doublets 
3) Remove non-expressing genes
4) Remove barcodes (cells) with low library complexity
5) Remove poor quality cells via the GMMchi-based filter using housekeeping genes
6) Remove barcodes (cells) with high level of mitochondrial rna expression
7) Normalize data

# Available Tools in the GMMchi_scs_pipeline Package

## Postprocess Input scRNA-seq Data

GMMchi_scs.GMMchi_scs_pipeline is the one-step postprocessing pipeline that takes in a dataframe with barcodes (row) x genes (columns) and returns a postprocessed dataframe with the same format

```
#run through your raw scRNA-seq data through the GMMchi single cell postprocessing pipeline
import GMMchi_scs_pipeline as GMMchi_scs #load the library

postprocessed_df = GMMchi_scs.GMMchi_scs_pipeline(cell_lines_scRNA)
```
## Map/Visualize Postprocessed Data Using UMAP

GMMchi_scs.UMAP_graph takes in a dataframe with barcdoes (row) x genes (columns) and returns a dataframe with barcodes (row) x UMAP features (columns) ready for visualization. For more information on UMAP please see [here](https://umap-learn.readthedocs.io/en/latest/basic_usage.html)

```
#map the postprocessed data with UMAP (dimensionality reduction technique)
UMAP_df = GMMchi_scs.UMAP_graph(postprocessed_df)
```
## Label UMAP for Downstream Visualization

GMMchi_scs.Label_graph is a built-in function that takes in a dataframe with barcodes (row) x features (columns). We've simplified the method so that users can quickly visualize their genes of interest easily. 

```
#use this if you just want the cells to be colored if the cell is expressing a gene above threshold 
Label_graph(postprocessed_df, UMAP_df, label_list=['ALPI']) 

# use this if you want the cells to be colored in according to the level of expression of the gene
Label_graph(postprocessed_df, UMAP_df, label_list=['ALPI'], boolean_visualization=False) 
```
# Working Example

Please find a working example in the example folder

## Authors

* **Ta-Chun (Jeff) Liu** - [jeffliu6068](https://github.com/jeffliu6068)
* **Sir Walter Fred Bodmer FRS FRSE** - *Supervision*

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration: Thank you for all that has contributed ideas and expertise to make this possible. Let's advance science together. 

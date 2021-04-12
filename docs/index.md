# VeloViz

<!-- badges: start -->
[![R build status](https://github.com/JEFworks/veloviz/workflows/R-CMD-check/badge.svg)](https://github.com/JEFworks/veloviz/actions)
<!-- badges: end -->

`VeloViz` creates RNA-velocity-informed 2D embeddings for single cell transcriptomics data.

<img src="https://github.com/JEFworks-Lab/veloviz/blob/package/docs/img/schematic_for_website.png?raw=true"/>


The overall approach is detailed in the [preprint](https://www.biorxiv.org/content/10.1101/2021.01.28.425293v1).

## Installation

To install `VeloViz`, we recommend using `devtools`:

``` r
require(devtools)
devtools::install_github('JEFworks-Lab/veloviz')
```

## Example

Below is a short example showing how to create a VeloViz embedding on sc-RNAseq data.   

``` r
library(veloviz)
# load built in scRNA-seq data
clusters = pancreas$clusters # cell type annotations
pcs = pancreas$pcs # PCs used to make other embeddings (UMAP,tSNE..)
vel = pancreas$vel # RNA velocity

curr = vel$current # current transcriptional state
proj = vel$projected # projected transcriptional state

# build VeloViz graph
veloviz = buildVeloviz(
curr = curr, proj = proj,
normalize.depth = TRUE,
use.ods.genes = TRUE,
alpha = 0.05,
pca = TRUE,
nPCs = 20,
center = TRUE,
scale = TRUE,
k = 5,
similarity.threshold = 0.25,
distance.weight = 1,
distance.threshold = 0.5,
weighted = TRUE,
seed = 0,
verbose = FALSE
)

# extract VeloViz embedding
emb.veloviz = veloviz$fdg_coords

# compare to other embeddings

par(mfrow = c(2,2))
#PCA
emb.pca = pcs[,1:2]
plotEmbedding(emb.pca, groups=pancreas$clusters, main='PCA')

#tSNE
set.seed(0)
emb.tsne = Rtsne::Rtsne(pcs, perplexity=30)$Y 
rownames(emb.tsne) = rownames(pcs)
plotEmbedding(emb.tsne, groups=pancreas$clusters, main='tSNE',
xlab = "t-SNE X", ylab = "t-SNE Y")

#UMAP
set.seed(0)
emb.umap = uwot::umap(pcs, min_dist = 0.5)
rownames(emb.umap) <- rownames(pcs)
plotEmbedding(emb.umap, groups=pancreas$clusters, main='UMAP',
xlab = "UMAP X", ylab = "UMAP Y")

#VeloViz
plotEmbedding(emb.veloviz, groups=clusters[rownames(emb.veloviz)], main='veloviz')

```
<img src="https://github.com/JEFworks-Lab/veloviz/blob/package/docs/img/readme_example.png?raw=true"/>

## Tutorials
[scRNA-seq data preprocessing and visualization using VeloViz](pancreas)  
[MERFISH cell cycle visualization using VeloViz](merfish)  
[Understanding VeloViz parameters](simulation) \
[Visualizing the VeloViz graph using UMAP](umap) \
[VeloViz with dynamic velocity estimates from scVelo](scVeloVignette)

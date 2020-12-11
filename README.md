
# VeloViz

<!-- badges: start -->
<!-- badges: end -->

`VeloViz` creates an RNA-velocity-informed 2D embedding for single cell transcriptomics data.

![](https://github.com/JEFworks/veloviz/blob/package/fig1a.png)

The overall approach is detailed in the following publication: [add link]

## Installation

To install `VeloViz`, we recommend using `devtools`:

``` r
require(devtools)
devtools::install_github('JEFworks/veloviz')
```

## Example

Below is a short example showing how to create a VeloViz embedding on sc-RNAseq data. In depth tutorials and examples using other data sets are available in our vignettes.  

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

# visualize VeloViz embedding
emb.veloviz = veloviz$fdg_coords
plotEmbedding(emb.veloviz, groups=clusters[rownames(emb.veloviz)], main='veloviz')

```

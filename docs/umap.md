Visualizing the VeloViz graph using UMAP
====================================

In this example, we will use Veloviz to produce a velocity-informed 2D
embedding and will then pass the computed nearest neighbor data into
UMAP. This is useful to users who want to use common algorithms such as
t-SNE and UMAP to layout their embedding. We will use the pancreas
endocrinogenesis dataset in this example.

First, load libraries:

    library(veloviz)
    library(velocyto.R)
    library(uwot)

Next, get pancreas data from VeloViz package:

    data(pancreas)

    spliced = pancreas$spliced
    unspliced = pancreas$unspliced
    clusters = pancreas$clusters
    pcs = pancreas$pcs

    #choose colors based on clusters for plotting later
    cell.cols = rainbow(8)[as.numeric(clusters)]
    names(cell.cols) = names(clusters)

Compute velocity:

    #cell distance in PC space
    cell.dist = as.dist(1-cor(t(pcs))) # cell distance in PC space

    vel = gene.relative.velocity.estimates(spliced,
                                           unspliced,
                                           kCells = 30,
                                           cell.dist = cell.dist,
                                           fit.quantile = 0.1)

Embed using UMAP with PCs as inputs:

    set.seed(0)
    emb.umap <- uwot::umap(pcs, min_dist = 0.5)
    rownames(emb.umap) <- rownames(pcs)

    plotEmbedding(emb.umap, colors = cell.cols, 
                  main = 'UMAP', xlab = "X", ylab = "Y")

![](umap_files/figure-markdown_strict/unnamed-chunk-5-1.png)

Now, build VeloViz graph:

    curr <- vel$current 
    proj <- vel$projected

    veloviz.graph <- buildVeloviz(
      curr = curr, 
      proj = proj,
      normalize.depth = TRUE,
      use.ods.genes = TRUE,
      alpha = 0.05,
      pca = TRUE,
      nPCs = 20,
      center = TRUE,
      scale = TRUE,
      k = 20,
      similarity.threshold = 0,
      distance.weight = 1,
      distance.threshold = 0,
      weighted = TRUE,
      seed = 0,
      verbose = FALSE
    )

    emb.veloviz <- veloviz.graph$fdg_coords
    plotEmbedding(emb.veloviz, colors = cell.cols, 
                  main = 'VeloViz with F-R', xlab = "X", ylab = "Y")

![](umap_files/figure-markdown_strict/unnamed-chunk-6-1.png)

Now, use UMAP to embed the velocity informed graph constructed using
VeloViz:

    veloviz.nnGraph <- asNNGraph(veloviz.graph) #converts veloviz igraph object to a format that UMAP understands 

    set.seed(0)
    emb.umapVelo <- uwot::umap(X = NULL, nn_method = veloviz.nnGraph, min_dist = 1)
    rownames(emb.umapVelo) <- rownames(emb.veloviz)
    plotEmbedding(emb.umapVelo, colors = cell.cols, 
                  main = 'VeloViz with UMAP', xlab = "X", ylab = "Y")

![](umap_files/figure-markdown_strict/unnamed-chunk-7-1.png)

    par(mfrow = c(1,3))

    plotEmbedding(emb.umap, colors = cell.cols, 
                  main = 'UMAP', xlab = "X", ylab = "Y")
    plotEmbedding(emb.veloviz, colors = cell.cols, 
                  main = 'VeloViz with F-R', xlab = "X", ylab = "Y")
    plotEmbedding(emb.umapVelo, colors = cell.cols, 
                  main = 'VeloViz with UMAP', xlab = "X", ylab = "Y")

![](umap_files/figure-markdown_strict/unnamed-chunk-8-1.png)

Let’s try it when there is a gap in the data. First, load the data from
the VeloViz package:

Next, get pancreas data from VeloViz package:

    data("pancreasWithGap")

    spliced = pancreasWithGap$spliced
    unspliced = pancreasWithGap$unspliced
    clusters = pancreasWithGap$clusters
    pcs = pancreasWithGap$pcs

    #choose colors based on clusters for plotting later
    cell.cols = rainbow(8)[as.numeric(clusters)]
    names(cell.cols) = names(clusters)

Compute velocity:

    #cell distance in PC space
    cell.dist = as.dist(1-cor(t(pcs))) # cell distance in PC space

    vel = gene.relative.velocity.estimates(spliced,
                                           unspliced,
                                           kCells = 30,
                                           cell.dist = cell.dist,
                                           fit.quantile = 0.1)

Embed using UMAP with PCs as inputs:

    set.seed(0)
    emb.umap <- uwot::umap(pcs, min_dist = 0.5)
    rownames(emb.umap) <- rownames(pcs)

    plotEmbedding(emb.umap, colors = cell.cols, 
                  main = 'UMAP', xlab = "X", ylab = "Y")

![](umap_files/figure-markdown_strict/unnamed-chunk-12-1.png)

Now, build VeloViz graph:

    curr <- vel$current 
    proj <- vel$projected

    veloviz.graph <- buildVeloviz(
      curr = curr, 
      proj = proj,
      normalize.depth = TRUE,
      use.ods.genes = TRUE,
      alpha = 0.05,
      pca = TRUE,
      nPCs = 20,
      center = TRUE,
      scale = TRUE,
      k = 20,
      similarity.threshold = 0,
      distance.weight = 1,
      distance.threshold = 0,
      weighted = TRUE,
      seed = 0,
      verbose = FALSE
    )

    emb.veloviz <- veloviz.graph$fdg_coords
    plotEmbedding(emb.veloviz, colors = cell.cols, 
                  main = 'VeloViz with F-R', xlab = "X", ylab = "Y")

![](umap_files/figure-markdown_strict/unnamed-chunk-13-1.png)

Now, use UMAP to embed the velocity informed graph constructed using
VeloViz:

    veloviz.nnGraph <- asNNGraph(veloviz.graph) #converts veloviz igraph object to a format that UMAP understands 

    set.seed(0)
    emb.umapVelo <- uwot::umap(X = NULL, nn_method = veloviz.nnGraph, min_dist = 1)
    rownames(emb.umapVelo) <- rownames(emb.veloviz)
    plotEmbedding(emb.umapVelo, colors = cell.cols, 
                  main = 'VeloViz with UMAP', xlab = "X", ylab = "Y")

![](umap_files/figure-markdown_strict/unnamed-chunk-14-1.png)

    par(mfrow = c(1,3))

    plotEmbedding(emb.umap, colors = cell.cols, 
                  main = 'UMAP', xlab = "X", ylab = "Y")
    plotEmbedding(emb.veloviz, colors = cell.cols, 
                  main = 'VeloViz with F-R', xlab = "X", ylab = "Y")
    plotEmbedding(emb.umapVelo, colors = cell.cols, 
                  main = 'VeloViz with UMAP', xlab = "X", ylab = "Y")

![](umap_files/figure-markdown_strict/unnamed-chunk-15-1.png)

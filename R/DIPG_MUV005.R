########
## scvelo analysis of DIPG patient MUV005
## @author: Jean Fan
## @date: August 17, 2020
########

######## Get data
dir <- '~/Desktop/DIPG/data/Pontine/'
emat <- read.csv(paste0(dir, 'spliced_reads.txt.gz'), sep="\t")
nmat <- read.csv(paste0(dir, 'unspliced_reads.txt.gz'), sep="\t")
#smat <- read.csv(paste0(dir, 'spanning_reads.txt.gz'), sep="\t")

library(Matrix)
emat <- Matrix(as.matrix(emat), sparse=TRUE)
nmat <- Matrix(as.matrix(nmat), sparse=TRUE)
#smat <- Matrix(as.matrix(smat), sparse=TRUE)

## parse out sample name
samples <- sapply(colnames(emat), function(x) strsplit(x, '[.]')[[1]][1])
table(samples)

## use one sample
vi <- samples == "MUV5"
emat.sub <- emat[,vi]
nmat.sub <- nmat[,vi]
#smat.sub <- smat[,vi]
samples.sub <- factor(samples[vi])
table(samples.sub)

######## Cluster analysis
library(MERINGUE)
cd <- emat.sub + nmat.sub
dim(cd)
mat <- normalizeCounts(cd)

## fast PCA with RSpectra
m <- t(mat)
m <- m - rowMeans(m)
set.seed(0)
pca <- RSpectra::svds(
  A    = m,
  k    = 20,
  opts = list(
    center = FALSE, scale = FALSE, maxitr = 2000, tol = 1e-10
  )
)
names(pca)

## look at elbow plot to check what is reasonable number of pcs
val <- pca$d
plot(val, type="l")
points(val)
N <- 10
abline(v=N, col='red')
pcs <- pca$u[, 1:N]
rownames(pcs) <- colnames(mat)
colnames(pcs) <- paste0('PC', 1:N)
head(pcs)

loading <- pca$v
rownames(loading) <- colnames(m)
head(loading)

pcs <- m %*% loading[,1:N]
pcs <- as.matrix(pcs)

## 2D visualization
library(uwot)
emb <- umap(pcs, min_dist = 0.5)
rownames(emb) <- rownames(pcs)

## visualize genes
par(mfrow=c(3,3))
genes <- c('MKI67', 'PCNA', 'CD44', 'APOE', 'ALDOC', 'PDGFRA', 'OLIG1', 'MBP', 'CLDN11')
genes %in% rownames(mat)
genes <- intersect(genes, rownames(mat))
sapply(genes, function(g) {
  gexp <- mat[g,]
  plotEmbedding(emb, colors=gexp, main=g)
})

## diffusion map
set.seed(0)
library(destiny)
foo <- as.matrix(pcs)
dm <- DiffusionMap(foo, k=30)
emb.dm <- eigenvectors(dm)[,c(1,2)]
par(mfrow=c(3,3))
sapply(genes, function(g) {
  gexp <- mat[g,]
  plotEmbedding(emb.dm, colors=gexp, main=g)
})

####### velocity
## restrict to protein coding genes
library(biomaRt)
ensembl=useMart("ensembl")
ensemblHuman = useDataset("hsapiens_gene_ensembl",mart=ensembl)
humanProteinCodingGenes = getBM(attributes=c("ensembl_gene_id","external_gene_name","description"),
                                filters='biotype', values=c('protein_coding'), mart=ensemblHuman)
head(humanProteinCodingGenes)
genes.have <- intersect(humanProteinCodingGenes$external_gene_name, rownames(mat))
length(genes.have)

library(velocyto.R)
cell.dist <- as.dist(1-cor(t(pcs))) ## cell distance in PC space
fit.quantile <- 0.05
rvel.cd <- gene.relative.velocity.estimates(emat.sub[genes.have,], nmat.sub[genes.have,],
                                            deltaT=1, kCells=10,
                                            cell.dist=cell.dist,
                                            fit.quantile=fit.quantile)
names(rvel.cd)
head(sort(rvel.cd$gamma, decreasing=TRUE), n=10)

## Plot a few genes
gene.relative.velocity.estimates(emat.sub, nmat.sub,
                                 kCells = 10,
                                 fit.quantile = fit.quantile,
                                 old.fit=rvel.cd,
                                 show.gene='YY1', ## oligo maturation TF
                                 cell.emb=emb.dm)
## Plot on embedding
show.velocity.on.embedding.cor(scale(emb), rvel.cd, n=100,
                               scale='sqrt',
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=2)
show.velocity.on.embedding.cor(scale(emb.dm), rvel.cd, n=100,
                               scale='sqrt',
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=2)

####### scvelo
library(reticulate)
conda_list()
use_condaenv("r-velocity", required = TRUE)
scv <- import("scvelo")
scv$logging$print_version()

ad <- import("anndata", convert = FALSE)
dat <- ad$AnnData(
  X=t(emat+nmat),
  obs=colnames(emat),
  var=rownames(emat),
  layers=list('spliced'=t(emat), 'unspliced'=t(nmat))
  )
dat$var
dat$obs
dat$layers

scv$pp$filter_genes(dat)
scv$pp$moments(dat)
scv$tl$recover_dynamics(dat)
## save since takes awhile
dat$write_h5ad('data/MUV005.h5ad')

## still giving me issues
#scv$tl$velocity(dat, mode='dynamical')
#scv$tl$velocity_graph(dat)
#scv$pl$velocity_embedding_stream(dat, basis='pca')

## look at top dynamic genes
topgenes <- dat$var["fit_likelihood"]
topgenes_vals <- py_to_r(topgenes$values)
names(topgenes_vals) <- rownames(cd)[as.numeric(py_to_r(topgenes$index$tolist()))]
head(sort(topgenes_vals, decreasing=TRUE), n=30)

g <- 'MAP2'
gene.relative.velocity.estimates(emat.sub, nmat.sub,
                                 kCells = 10,
                                 fit.quantile = fit.quantile,
                                 old.fit=rvel.cd,
                                 show.gene=g,
                                 cell.emb=emb.dm)

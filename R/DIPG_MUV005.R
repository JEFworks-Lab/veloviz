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
library(MUDAN) ## helper https://github.com/JEFworks/MUDAN
cd <- emat.sub + nmat.sub
dim(cd)

## restrict to protein coding genes (for now...need to follow up on how data was processed)
## unclear why some genes have no variance
## also spanning reads should be subset of unspliced but doesn't seem to be
library(biomaRt)
ensembl=useMart("ensembl")
ensemblHuman = useDataset("hsapiens_gene_ensembl",mart=ensembl)
humanProteinCodingGenes = getBM(attributes=c("ensembl_gene_id","external_gene_name","description"),
                                filters='biotype', values=c('protein_coding'), mart=ensemblHuman)
head(humanProteinCodingGenes)
genes.have <- intersect(humanProteinCodingGenes$external_gene_name, rownames(cd))
length(genes.have)
cd <- cd[genes.have,]

cd <- cleanCounts(cd, max.lib.size = Inf, min.reads = 1)
mat <- normalizeCounts(cd)
matnorm <- normalizeVariance(cd, plot=TRUE, details=TRUE)
odgenes <- rownames(mat)[matnorm$ods] ## overdispersed genes
length(odgenes)

## fast PCA with RSpectra
m <- t(mat[odgenes,])
m <- t(t(m) - colMeans(m))
colMeans(m) ## double check mean 0
v <- apply(m, 2, var)
m2 <- t(t(m)/sqrt(v))
apply(m2, 2, var) ## double check var 1
set.seed(0)
pca <- RSpectra::svds(
  A    = m2,
  k    = 100,
  opts = list(
    center = FALSE, scale = FALSE, maxitr = 2000, tol = 1e-10
  )
)
names(pca)

## look at elbow plot to check what is reasonable number of pcs
val <- pca$d
plot(val, type="l")
points(val)
N <- 30
abline(v=N, col='red')
pcs <- pca$u[, 1:N]
rownames(pcs) <- colnames(mat)
colnames(pcs) <- paste0('PC', 1:N)
head(pcs)

## can save loadings for projection later
loading <- pca$v
rownames(loading) <- colnames(m)
pcs <- m2 %*% loading[,1:N]
pcs <- as.matrix(pcs)
dim(pcs)

## 2D visualization
library(uwot)
emb <- umap(pcs, min_dist = 0.5)
rownames(emb) <- rownames(pcs)

## visualize genes
par(mfrow=c(4,3))
genes <- c('MKI67', 'PCNA', 'TOP2A', ## proliferation markers
           'CD44', 'APOE', 'ALDOC', ## astrocyte markers
           'PDGFRA', 'CSPG4', 'OLIG1', ## oligodendrocyte progenitor markers
           'MBP', 'CLDN11', 'PLP1' ## mature oligodendrocyte markers
           )
genes %in% rownames(mat)
genes <- intersect(genes, rownames(mat))
sapply(genes, function(g) {
  gexp <- mat[g,]
  plotEmbedding(emb, colors=gexp, main=g)
})

## graph based clustering
com <- getComMembership(pcs, k=30, method=igraph::cluster_louvain)
table(com)
plotEmbedding(emb, groups=com, mark.clusters=TRUE)

## manual annotation
annot <- as.character(com)
annot[com == 1] <- 'SC' ## stem cell
annot[com  %in% c(2)] <- 'AS' ## astrocytic
annot[com %in% c(4,5,6)] <- 'OC' ## oligodendrocytic
annot[com == 3] <- 'mOC' ## myelinating/mature oligodendritic?
names(annot) <- names(com)
plotEmbedding(emb, groups=annot, mark.clusters=TRUE)

## diffusion map
set.seed(0)
library(destiny)
foo <- as.matrix(pcs)
dm <- DiffusionMap(foo, k=30)
emb.dm <- eigenvectors(dm)[,c(1,2)]
par(mfrow=c(4,3))
sapply(genes, function(g) {
  gexp <- mat[g,]
  plotEmbedding(emb.dm, colors=gexp, main=g)
})
plotEmbedding(emb.dm, groups=annot, mark.clusters=TRUE)
cell.cols <- MUDAN:::fac2col(annot)

####### velocity
library(velocyto.R)
cell.dist <- as.dist(1-cor(t(pcs))) ## cell distance in PC space
fit.quantile <- 0.05
## restrict to protein coding genes
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
                                 cell.emb=emb.dm,
                                 cell.colors = cell.cols)
gene.relative.velocity.estimates(emat.sub, nmat.sub,
                                 kCells = 10,
                                 fit.quantile = fit.quantile,
                                 old.fit=rvel.cd,
                                 show.gene=genes[8],
                                 cell.emb=emb.dm,
                                 cell.colors = cell.cols)

## Plot on embedding
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

## create AnnData object
ad <- import("anndata", convert = FALSE)
dfobs <- data.frame(clusters=annot)
dfvar <- data.frame(genes.have)
rownames(dfvar) <- genes.have
dat <- ad$AnnData(
  X=t(emat.sub[genes.have,]+nmat.sub[genes.have,]),
  obs=dfobs,
  var=dfvar,
  layers=list('spliced'=t(emat.sub[genes.have,]), 'unspliced'=t(nmat.sub[genes.have,])),
  obsm=list('X_diffusion'=emb.dm, 'X_umap'=emb) ## can add more embeddings
  )
## double check
dat$var
dat$obs
dat$layers
dat$obsm

## run scvelo dynamic model
scv$pp$filter_genes(dat)
scv$pp$moments(dat)
scv$tl$recover_dynamics(dat)
## save since takes awhile
dat$write_h5ad('data/MUV005.h5ad')
#dat = scv$read('data/MUV005.h5ad')

## plot (creates pop up window)
scv$tl$velocity(dat, mode='dynamical')
scv$tl$velocity_graph(dat)
scv$pl$velocity_embedding_stream(dat, basis='diffusion')

## latent time
scv$tl$latent_time(dat)
lt <- py_to_r(dat$obs$latent_time)
names(lt) <- py_to_r(dat$obs$names)[names(lt)]
head(lt)
hist(lt)
foo <- scale(exp(lt))[,1] ## transform?
range(foo)
foo[foo < -1] <- -1
foo[foo > 1] <- 1
plotEmbedding(emb.dm, col=foo)

## top dynamic genes
topgenes <- dat$var["fit_likelihood"]
topgenes_vals <- py_to_r(topgenes$values)
names(topgenes_vals) <- py_to_r(dat$var_names$values)
head(topgenes_vals)
head(sort(topgenes_vals, decreasing=TRUE), n=30)

g <- 'GPRC5C'
gene.relative.velocity.estimates(emat.sub, nmat.sub,
                                 kCells = 10,
                                 fit.quantile = fit.quantile,
                                 old.fit=rvel.cd,
                                 show.gene=g,
                                 cell.emb=emb.dm,
                                 cell.colors = cell.cols)

## group by clusters
scv$tl$rank_dynamical_genes(dat, groupby='clusters')
df = scv$get_df(dat, 'rank_dynamical_genes/names')
head(df)

g <- 'CELA3A'
gene.relative.velocity.estimates(emat.sub, nmat.sub,
                                 kCells = 10,
                                 fit.quantile = fit.quantile,
                                 old.fit=rvel.cd,
                                 show.gene=g,
                                 cell.emb=emb.dm,
                                 cell.colors = cell.cols)
g <- 'MMP16'
gene.relative.velocity.estimates(emat.sub, nmat.sub,
                                 kCells = 10,
                                 fit.quantile = fit.quantile,
                                 old.fit=rvel.cd,
                                 show.gene=g,
                                 cell.emb=emb.dm,
                                 cell.colors = cell.cols)

## differential kinetics
genes
scv$tl$differential_kinetic_test(dat, var_names=genes, groupby='clusters')
## insufficient cells since only one sample

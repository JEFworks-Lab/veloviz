#### SETUP ####
#use velocyto2 environment 
setwd("/Users/lylaatta/Documents/GitHub/veloviz/tutorial_lyla")
library(MUDAN)
library(sva)
library(RSpectra)
library(velocyto.R)
library(Rtsne)
library(destiny)
library(umap)

#### LOAD DATA ####
cellExp = read.csv("../data/S12_cell_gexp.csv", header = TRUE)
rownames(cellExp) = cellExp$X
cellExp = cellExp[,2:ncol(cellExp)]

nucExp = read.csv("../data/S14_nuc_gexp.csv", header = TRUE)
rownames(nucExp) = nucExp$X
nucExp = nucExp[,2:ncol(nucExp)]

cytoExp = cellExp - nucExp

geneInfo = read.csv("../data/S1_codebook.csv", header = TRUE, stringsAsFactors = FALSE)
longGenes = geneInfo[1:9050,'name'] #non-overlapping probes
shortGenes = geneInfo[9279:10278,]$name #overlapping probes
blankGenes = geneInfo[grepl('Blank',geneInfo$name),]$name #blank controls 

#only use long genes for analysis 
testGenes = as.character(longGenes)
testCell = cellExp[testGenes,]
testNuc = nucExp[testGenes,]
testCyto = cytoExp[testGenes,]

#### CLEANUP ####
cd = testCell

#batch annotation 
batch = sapply(colnames(cd), function(x) strsplit(x, '_')[[1]][1])
batch = factor(batch)

#filter genes with low expression 
genesToKeep = rowMeans(cd) > 1
cd = cd[genesToKeep,]

#batch correct
mat.bc = ComBat(as.matrix(cd),batch[colnames(cd)])
mat.bc[mat.bc<0] = 0 #set negative expression to 0

#CPM normalize 
mat.cpm = normalizeCounts(mat.bc)

#Normalize variance 
mat.norm = normalizeVariance(mat.cpm, details = TRUE, plot = FALSE)

#Over dispersed genes, log normalize
mat = log10(as.matrix(t(mat.norm$mat[mat.norm$ods,]))+1)
#saveRDS(mat,file = "mat_u2os.rds")

#### DIM REDUCTION ####
#PCA 
pca = svds(A = mat, k = 50, opts = list(center = TRUE, scale = FALSE, maxitr = 2000, tol = 1e-10))
#saveRDS(pca,file = "pca.rds")
pcaInfo = list()
pcaInfo[['pca']] = pca
pcaInfo[['cellNames']] = colnames(mat.bc)
pcaInfo[['geneNames']] = colnames(mat)
#saveRDS(pcaInfo,file = "pcaInfo.rds")

val = pca$d
plot(val)
points(val)

N = 30 
abline(v=13, col = 'red')

pcs = pca$u[,1:N]
rownames(pcs) = colnames(mat.bc)
colnames(pcs) = paste0('PC',1:N)

#low dimension embedding 
emb.test = pcs[,1:2]


#community detection 
com = getComMembership(pcs, k=300, method = igraph::cluster_louvain)
cluster.label = factor(com)
cell.color = MUDAN:::fac2col(cluster.label)

par(mfrow = c(1,1),mar=rep(5,4))
plotEmbedding(emb.test, groups = com, mark.clusters = TRUE, show.legend = TRUE, xlab = 'PC1', ylab = 'PC2', 
              main = 'clustering', verbose = FALSE)
plot(scale(emb.test), pch=16, xlab = 'PC1', ylab = 'PC2', col = cell.color)
plotEmbedding(emb.test, groups=batch, show.legend = TRUE,
              xlab = 'PC1', ylab = 'PC2', main = 'batch', verbose = FALSE)

## Plot cell-cycle and house keeping genes noted in the manuscript
par(mfrow=c(3,3), mar=rep(1,4))
gs <- c('MCM5','SNN','UNG','MCM6','DSCC1','BCL2L1','CCNF','KIF2C','PPIE') 
invisible(lapply(gs, function(g) {
  gexp = scale(log10(mat.cpm[g,]+1))[,1]
  gexp[gexp > 2] <- 2
  gexp[gexp < -2] <- -2
  plotEmbedding(emb.test, main=g, col=gexp, verbose=FALSE)
}))

cluster.label = factor(com)
cell.color = MUDAN:::fac2col(cluster.label)

#### VELOCITY ####

subcells = names(batch)
emat = as.matrix(testCyto[,subcells])
nmat = as.matrix(testNuc[,subcells])

cell.dist = as.dist(1-cor(t(pcs[subcells,])))
fit.quantile = 0.05
# rvel.cd = gene.relative.velocity.estimates(emat, nmat, deltaT = 1, kCells = 30, cell.dist = cell.dist, fit.quantile = fit.quantile)
# saveRDS(rvel.cd,file = "rvelcd.rds")
rvel.cd = readRDS("rvelcd.rds")

gene.relative.velocity.estimates(emat, nmat, kCells = 30,
                                 fit.quantile = fit.quantile,
                                 old.fit=rvel.cd,
                                 show.gene='KIF2C',
                                 cell.emb=emb.test,
                                 cell.colors=cell.color)

#plot velocity projections on PCs
show.velocity.on.embedding.cor(scale(emb.test), rvel.cd, n=100, scale='sqrt', cell.colors=cell.color,
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1)


#plot velocity projections on tsne embedding 
set.seed(1)
emb.tsne = Rtsne(pcs, is_distance = FALSE, perplexity = 10, num_threads = 1, verbose = FALSE)$Y
rownames(emb.tsne) = rownames(pcs)

show.velocity.on.embedding.cor(scale(emb.tsne), rvel.cd, n=100, scale='sqrt', cell.colors=cell.color,
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=2)


#plot velocity projections on diffusion mapping 
emb.destiny = DiffusionMap(pcs)
emb.destiny = eigenvectors(emb.destiny)[,1:2]

show.velocity.on.embedding.cor(scale(emb.destiny), rvel.cd, n=100, scale='sqrt', cell.colors=cell.color,
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=2)

#plot velocity projection on umap embedding 

emb.umap = umap(pcs)
plot(scale(emb.umap$layout), pch=16, xlab = 'UMAP X', ylab = 'UMAP Y', col = cell.color)

show.velocity.on.embedding.cor(scale(emb.umap$layout), rvel.cd, n=100, scale='sqrt', cell.colors=cell.color,
                               cex=1, arrow.scale=1, show.grid.flow=TRUE,
                               min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=2)






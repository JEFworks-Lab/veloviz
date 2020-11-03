########### quantifying trajectory contiguity across simulated gaps in the visualized embedding
library(velocyto.R)
library(reticulate)
library(Matrix)
library(MUDAN)
conda_list()
use_condaenv("r-velocity", required = TRUE)
scv = import("scvelo")

adata = scv$datasets$pancreas()

spliced = as.matrix(t(adata$layers['spliced']))
unspliced = as.matrix(t(adata$layers['unspliced']))
cells = adata$obs_names$values
genes = adata$var_names$values
colnames(spliced) = colnames(unspliced) = cells
rownames(spliced) = rownames(unspliced) = genes

clusters = adata$obs$clusters #extract clusters
names(clusters) = adata$obs_names$values
col = rainbow(length(levels(clusters)),s = 0.8, v = 0.8)
cell.cols = col[clusters] #color according to cluster
names(cell.cols) = names(clusters)

## old embedding
emb_umap = adata$obsm['X_umap'] #extract umap embedding
rownames(emb_umap) <- names(clusters)
par(mfrow = c(1,1))
plotEmbedding(emb_umap, groups=clusters, xlab = "UMAP X", ylab = "UMAP Y", mark.clusters=TRUE)
#plot(emb_umap)

## downsample for faster runs
set.seed(0)
good.cells <- sample(rownames(emb_umap), nrow(emb_umap)/5)
length(good.cells )

plotEmbedding(emb_umap[good.cells ,], groups=clusters, xlab = "UMAP X", ylab = "UMAP Y", mark.clusters=TRUE)
emb_umap = emb_umap[good.cells ,]
clusters = clusters[good.cells]
cell.cols = cell.cols[good.cells]
spliced = spliced[,good.cells]
unspliced = unspliced[,good.cells]

# remove cells to create trajectory gap
x = emb_umap[,1]
vi <- x > -5 & x < 0
table(vi)
good.cells <- rownames(emb_umap)[!vi]
plotEmbedding(emb_umap[good.cells,], groups=clusters, xlab = "UMAP X", ylab = "UMAP Y", mark.clusters=TRUE)
emb_umap = emb_umap[good.cells ,]
clusters = clusters[good.cells]
cell.cols = cell.cols[good.cells]
spliced = spliced[,good.cells]
unspliced = unspliced[,good.cells]

plot(emb_umap)
## cells before gap
x = emb_umap[,1]
cells.before.gap <- rownames(emb_umap)[x < -5 & x > -6]
points(emb_umap[cells.before.gap,], col='red')
## cells after gap
cells.after.gap <- rownames(emb_umap)[x > 0 & x < 1]
points(emb_umap[cells.after.gap,], col='blue')

############# Velocity analysis
## filter
gexpS = log10(rowSums(spliced)+1)
gexpU = log10(rowSums(unspliced)+1)
par(mfrow = c(2,2))
hist(gexpS, breaks = 500, main = "log spliced counts")
hist(gexpU, breaks = 500, main = "log unspliced counts")

## less stringent?
goodGenes = genes[gexpS > 1 & gexpU > 1]
spliced = spliced[goodGenes,]
unspliced = unspliced[goodGenes,]
dim(spliced)
dim(unspliced)

### analyze
allExp = spliced + unspliced #use combined spliced and unspliced counts
all.cpm = MUDAN::normalizeCounts(allExp) #cpm normalize
all.norm = MUDAN::normalizeVariance(all.cpm,details = TRUE, plot = FALSE) #variance stabilize
odsGenes = goodGenes[all.norm$ods] #overdispersed genes
all.logODS = log10(as.matrix(t(all.norm$mat[all.norm$ods,]))+1) #keep overdispersed genes and log normalize

## PCA
pca = RSpectra::svds(A = all.logODS, k=50, opts = list(center = TRUE, scale = FALSE, maxitr = 2000, tol = 1e-10))
cell.dist = as.dist(1-cor(t(pca$u))) #cell distance in PC space

library(velocyto.R)
vel = gene.relative.velocity.estimates(spliced,unspliced,deltaT=1,kCells=30,cell.dist=cell.dist,fit.quantile=0.1,mult=100)
names(vel)

############# Veloviz
#getwd()
#source('../R/main.R')
mat = vel$conv.emat.norm + vel$conv.nmat.norm
panc_k30 <- veloviz(mat,
                    vel$current, vel$projected,
                    cell.cols[colnames(mat)],
                    k = 15, weighted=TRUE,
                    similarity.threshold = -0.5, seed = 0)

panc_k30$fdg_coords = scale(panc_k30$fdg_coords, center=FALSE, scale=TRUE)

plotEmbedding(panc_k30$fdg_coords, col=cell.cols)
points(panc_k30$fdg_coords[cells.before.gap,], col='red')
points(panc_k30$fdg_coords[cells.after.gap,], col='blue')

x1 = colMeans(panc_k30$fdg_coords[cells.before.gap,])
x2 = colMeans(panc_k30$fdg_coords[cells.after.gap,])
points(rbind(x1, x2), col='black', pch=16, cex=2)
dist(rbind(x1, x2))

library(uwot)
set.seed(0)
emb <- umap(pca$u[,1:10], min_dist = 0.5)
rownames(emb) <- colnames(allExp)
head(emb)
emb = scale(emb, center=FALSE, scale=TRUE)

plotEmbedding(emb, groups=clusters, xlab = "UMAP X", ylab = "UMAP Y", mark.clusters=FALSE)
points(emb[cells.before.gap,], col='red')
points(emb[cells.after.gap,], col='blue')

y1 = colMeans(emb[cells.before.gap,])
y2 = colMeans(emb[cells.after.gap,])
points(rbind(y1, y2), col='black', pch=16, cex=2)
dist(rbind(y1, y2))


############ Repeat
## not the most efficient due to redoing normalizations, but try it out for now
results <- do.call(rbind, lapply(1:20, function(i) {
  print(i)
  panc_k30 <- veloviz(mat,
                      vel$current, vel$projected,
                      cell.cols[colnames(mat)],
                      k = 15, weighted=TRUE,
                      similarity.threshold = -0.5,
                      seed = i,
                      plot = FALSE)

  panc_k30$fdg_coords = scale(panc_k30$fdg_coords, center=FALSE, scale=TRUE)

  plotEmbedding(panc_k30$fdg_coords, col=cell.cols)
  points(panc_k30$fdg_coords[cells.before.gap,], col='red')
  points(panc_k30$fdg_coords[cells.after.gap,], col='blue')

  x1 = colMeans(panc_k30$fdg_coords[cells.before.gap,])
  x2 = colMeans(panc_k30$fdg_coords[cells.after.gap,])
  points(rbind(x1, x2), col='black', pch=16, cex=2)
  r1 = dist(rbind(x1, x2))
  print(r1)

  set.seed(i)
  emb <- umap(pca$u[,1:10], min_dist = 0.5)
  rownames(emb) <- colnames(allExp)
  head(emb)
  emb = scale(emb, center=FALSE, scale=TRUE)

  plotEmbedding(emb, groups=clusters, xlab = "UMAP X", ylab = "UMAP Y", mark.clusters=FALSE)
  points(emb[cells.before.gap,], col='red')
  points(emb[cells.after.gap,], col='blue')

  y1 = colMeans(emb[cells.before.gap,])
  y2 = colMeans(emb[cells.after.gap,])
  points(rbind(y1, y2), col='black', pch=16, cex=2)
  r2 = dist(rbind(y1, y2))
  print(r2)

  return(c(r1, r2))
}))
colnames(results) <- c('veloviz', 'umap')
head(results)

par(mfrow=c(2,1))
hist(results[,1],main='veloviz distance', xlim=c(0, max(results)+1))
hist(results[,2],main='umap distance', xlim=c(0, max(results)+1), breaks=20)


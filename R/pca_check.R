########### PCA check

########### Get MERFISH data
getwd()
dir <- 'data/'
list.files(dir)
cell_gexp <- as.matrix(read.csv(paste0(dir, 'S12_cell_gexp.csv.gz'), row.names=1))
print(cell_gexp[1:5,1:5])
nuc_gexp <- as.matrix(read.csv(paste0(dir, 'S14_nuc_gexp.csv.gz'), row.names=1))
print(nuc_gexp[1:5,1:5])
cyto_gexp <- cell_gexp - nuc_gexp
print(cyto_gexp[1:5,1:5])

gene_info <- read.csv(paste0(dir, 'S1_codebook.csv.gz'), header=FALSE, stringsAsFactors = FALSE)
long.genes <- gene_info[2:9051,1]
short.genes <- gene_info[9280:10279,1]
bad.genes <- gene_info[,1][grepl('Blank', gene_info[,1])]

length(long.genes)
length(short.genes)
length(bad.genes)

test.genes <- long.genes ## use subset of genes as in original paper
cell_gexp <- cell_gexp[test.genes,]
nuc_gexp <- nuc_gexp[test.genes,]
cyto_gexp <- cyto_gexp[test.genes,]

spliced = cyto_gexp
unspliced = nuc_gexp

## normalize
allExp = spliced + unspliced #use combined spliced and unspliced counts
all.cpm = MUDAN::normalizeCounts(allExp) #cpm normalize
all.norm = MUDAN::normalizeVariance(all.cpm, details = TRUE, plot = FALSE) #variance stabilize
odsGenes = rownames(allExp)[all.norm$ods] #overdispersed genes
length(odsGenes)

################### double check that variance has already been normalized
## to take into consideration the dependency between expression magnitude and variance
par(mfrow=c(1,2))
## before
plot(
  log10(apply(all.cpm, 1, mean)),
  log10(apply(all.cpm, 1, var)),
  pch="."
)
points(
  log10(apply(all.cpm, 1, mean))[odsGenes],
  log10(apply(all.cpm, 1, var))[odsGenes],
  col='red',
  pch="."
)
## after
varNormMat <- all.norm$mat
plot(
  log10(apply(varNormMat, 1, mean)),
  log10(apply(varNormMat, 1, var)),
  pch="."
)
points(
  log10(apply(varNormMat, 1, mean))[odsGenes],
  log10(apply(varNormMat, 1, var))[odsGenes],
  col='red',
  pch="."
)
## so future PCA will not be dominated by highly expressed genes

####################### compare PCAs
all.logODS <- log10(varNormMat[odsGenes,]+1) ##overdispsersed genes only

pca0 = prcomp(t(all.logODS), center = TRUE, scale = FALSE) ##slow
names(pca0)
pca0$x[1:5,1:5]
pca0$rotation[1:5,1:5]

pca1 = svds(A = t(all.logODS), k=50, opts = list(center = TRUE, scale = FALSE, maxitr = 2000, tol = 1e-10))
names(pca1)
m = all.logODS - apply(all.logODS, 1, mean) ## center
apply(m, 1, mean) ## check
pca1.x =  t(m) %*% pca1$v
pca1.x[1:5,1:5]
pca1$v[1:5,1:5]

plot(pca0$x[,1], pca1.x[,1]) ## should be correlated, signs may flip
plot(pca0$x[,1], pca1$u[,1]) ## should be correlated, scale will differ
plot(pca1$u[,1], pca1.x[,1]) ## should be correlated, scale will differ

## what happens if we fail to center
pca1.x.wrong =  t(all.logODS) %*% pca1$v
pca1.x.wrong[1:5,1:5]

plot(pca0$x[,1], pca1.x.wrong[,1])
plot(pca1$u[,1], pca1.x.wrong[,1])
## still correlated, but values just no longer centered
## so PCs are not zero centered
## not as bad as I expected for some reason

## center myself
m = all.logODS - apply(all.logODS, 1, mean)
pca2 = svds(A = t(m), k=50, opts = list(center = FALSE, scale = FALSE, maxitr = 2000, tol = 1e-10))
pca2.x =  t(m) %*% pca2$v
pca2.x[1:5,1:5]
pca2$v[1:5,1:5]

plot(pca1.x[,1], pca2.x[,1]) ## should be nearly exactly the same (some stochasticity due to approximations)
table(pca1.x[,1] == pca2.x[,1])

####################### Rest of analysis
### velocity
vel <- gene.relative.velocity.estimates(spliced, unspliced,
                                        deltaT=1, kCells=30,
                                        cell.dist=cell.dist)
vel$kCells
curr = vel$current
proj = vel$projected

### TODO: given velocity object, clean up into functions to
### 1. normalize both curr and projected variance
### 2. project both curr and projected into same PC space
### 3. test for optimal parameters?
### 4. visualize

## normalize variance on current
## use same normalizing factor for projected
foo <- normalizeVariance(curr, details=TRUE, plot=TRUE)
names(foo)
head(foo$df)
## double check
ods.genes <- rownames(foo$mat)[foo$ods]

par(mfrow=c(2,2))
plot(foo$df$m, foo$df$v)
points(foo$df[ods.genes,]$m, foo$df[ods.genes,]$v, col='red')
plot(foo$df$m, foo$df$res)
points(foo$df[ods.genes,]$m, foo$df[ods.genes,]$res, col='red')

## variance normalize both
## the current and projected transcriptional states
## by same scaling factor
curr.sub <- curr[ods.genes,] * foo$df[ods.genes,]$gsf
table(curr.sub[,1] == foo$mat[ods.genes,1]) ## double check
proj.sub <- proj[ods.genes,] * foo$df[ods.genes,]$gsf

plot(apply(curr.sub, 1, mean), apply(curr.sub, 1, var))
plot(apply(proj.sub, 1, mean), apply(proj.sub, 1, var))

curr.sub <- log10(curr.sub+1)
proj.sub <- log10(proj.sub+1)

library(RSpectra)
pca = svds(A = t(curr.sub), k=50, opts = list(center = TRUE, scale = FALSE, maxitr = 2000, tol = 1e-10))
range(pca$u)
var = pca$d
plot(var)

nPCs <- 20

m <- t(curr.sub)
m <- m - apply(m, 1, mean)
range(apply(m, 1, mean))
pca.curr <- m %*% pca$v[,1:nPCs]
rownames(pca.curr) <- colnames(curr.sub)
head(pca.curr)

m <- t(proj.sub)
m <- m - apply(m, 1, mean)
range(apply(m, 1, mean))
pca.proj <- m %*% pca$v[,1:nPCs]
rownames(pca.proj) <- colnames(proj.sub)
head(pca.proj)

colnames(pca.curr) <- colnames(pca.proj) <- paste0('PC', 1:ncol(pca.curr))

## double check zero centered
par(mfrow=c(2,2))
hist(pca.curr)
hist(pca.proj)

################ our visualization
library(igraph)
library(matie)
library(RANN)
library(Rcpp)
getwd()
setwd('graphViz')
#source('projectedNeighbors.R')
source("projectedNeighbors_weightedCD.R")

gsim.gg = graphViz(t(pca.curr), t(pca.proj), 30,
                   cell.colors=cell.colors,
                   similarity_threshold=0,
                   distance_weight = 1,
                   weighted=FALSE,
                   plot = FALSE,
                   return_graph = TRUE)
par(mfrow=c(1,1))
plot(gsim.gg$fdg_coords, main = "FDG: vertex coordinates", col=cell.colors, pch=16)

par(mfrow=c(2,2), mar=rep(1,4))
m <- log10(all.cpm+1)
gs = c('MCM5','MCM6','CCNF','KIF2C')
invisible(lapply(gs, function(g){
  gexp = scale(m[g,])[,1]
  gexp[gexp > 2] <- 2
  gexp[gexp < -2] <- -2
  plotEmbedding(gsim.gg$fdg_coords,
                       main=g,
                       col=gexp,
                       verbose=FALSE)
}))

par(mfrow=c(1,1))
g <- gsim.gg$graph
V(g)$label = NA
V(g)$size = 2
E(g)$arrow.size = 1
E(g)$color = rgb(0,0,0,0.01)
plot(g)

deltaExp = vel$deltaE
scores = consistency(gsim.gg$fdg_coords, deltaExp, 20, FALSE)
mean(scores)

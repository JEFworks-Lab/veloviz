#' Get projected neighbors
#'
#' @export
#'
projectedNeighbors = function(observed,projected,k,distance_metric="L2",similarity_metric="cosine",distance_weight = 1, distance_threshold = 1, similarity_threshold = -1){
  #observed: genes/PCs (rows) x cells (columns) matrix of observed cells 
  #projected: genes/PCs (rows) x cells (columns) matrix of projected states of cells in observed (same order)
  #k: number of nearest neighbors 
  #distance_metric: "L1" or "L2" to calculate cell-cell distance, d
  #similarity_metric: "cosine" or "pearson" correlation to calculate difference x velocity similarity, sim. NOTE: pearson similarity behaves weird with two dimensions 
  #distance_weight: weight w in inverse distance component of composite distance: CD = [1/(1+w*d)]*[- sim]
  #distance threshold: [0,1], proportion of largest inverse distances to include when determining which edges to include. e.g. 0.1 --> keeps edges with top 10% of inverse distances 
  #similarity_threshold: minimum similarity between velocity vector and cell_i-->nn_i vector for edge to be included. default: all edges included. 
  
  observed = t(observed)
  projected = t(projected)
  n = nrow(observed)
  
  nn_idx = matrix(NA,nrow=n, ncol = 1) #projected nearest neighbor 
  knn_idx = matrix(NA,nrow=n, ncol = k) #projected k nearest neighbors 
  all_dists = matrix(NA,nrow = n, ncol = n)
  all_invDist = matrix(NA,nrow = n, ncol = n)
  all_negVectSim = matrix(NA,nrow = n, ncol = n)
  
  for (i in seq(1,nrow(observed))){
    cell_i = observed[i,] #current cell whose projected neighbor we want to find
    proj_i = projected[i,] #projected state of cell i 
    obs_exc_i = observed[-i,] #cell_i can't be its own projected neighbor - remove from consideration
    
    nn_dists_i = matrix(NA, nrow=nrow(obs_exc_i), ncol = 1) #distances of all other observed cells from cell_i
    nn_invDist_i = matrix(NA, nrow=nrow(obs_exc_i), ncol = 1) #inverse distance component
    nn_negVectSim_i = matrix(NA, nrow=nrow(obs_exc_i), ncol = 1) #negative similarity component
    
    #calculate distances between cell_i and all other cells j 
    cell_i_dists = pwiseDists(cell_i,proj_i,as.matrix(obs_exc_i),distance_metric,similarity_metric,distance_weight)
    nn_dists_i = cell_i_dists[,"CompositeDistance"]
    nn_invDist_i = cell_i_dists[,"InverseDistance"]
    nn_negVectSim_i = cell_i_dists[,"NegativeSimilarity"]
    
    #add cell_i's dists to all_dists 
    all_dists[i,-i] = nn_dists_i
    all_invDist[i,-i] = nn_invDist_i
    all_negVectSim[i,-i] = nn_negVectSim_i
    
    #find index of minimum distance between cell_i and other cells - those will become k nearest neighbots knn_i 
    min_dist_idx = which(nn_dists_i==min(nn_dists_i))
    #find indices of k nearest neighbors 
    k_min_dists_idx = order(nn_dists_i)[1:k]
    #exclude neighbor if similarity below threshold 
    sim_nn_k = nn_negVectSim_i[k_min_dists_idx] #similarities of minimum distance neighbors 
    new_k_idx = k_min_dists_idx[sim_nn_k <= (-1*similarity_threshold)] #indices of minimum distance neighbors excluding those whose similarities are above threshold 
    
    #correct min_idx - will be off by 1 if neighbor idx is greater than cell idx 
    min_dist_idx[which(min_dist_idx>=i)] = min_dist_idx[which(min_dist_idx>=i)] + 1
    #correct min k idx 
    k_min_dists_idx[which(k_min_dists_idx>=i)] = k_min_dists_idx[which(k_min_dists_idx>=i)] + 1
    if (length(new_k_idx)>0){  #### CHANGED HERE
      new_k_idx[which(new_k_idx>=i)] = new_k_idx[which(new_k_idx>=i)] + 1
    }
    
    #add current cell's neighbors to all cell neighbors matrix
    nn_idx[i] = min_dist_idx
    #knn_idx[i,] = k_min_dists_idx
    if (length(new_k_idx)>0){
      knn_idx[i,c(1:length(new_k_idx))] = new_k_idx
    }
  }
  
  #threshold distances 
  n_keep = ceiling(distance_threshold*length(all_invDist)) # number of edges to keep based on prop in distance_threshold
  dist_cutoff = sort(all_invDist,decreasing = TRUE)[n_keep] # distance cutoff for edges to keep 
  inv_dist_knn = t(sapply(c(1:n), function(x) all_invDist[x,knn_idx[x,]]))
  
  dist_comp = list()
  dist_comp[["invDist"]] = all_invDist
  dist_comp[["negVectSim"]] = all_negVectSim
  ## adding edge weights 
  edge_weights = t(sapply(c(1:n), function(x) all_dists[x,knn_idx[x,]]))
  edge_weights_dt = edge_weights
  edge_weights_dt[inv_dist_knn<dist_cutoff] = NA
  knn_idx_dt = knn_idx
  knn_idx_dt[inv_dist_knn<dist_cutoff] = NA
  
  out = list()
  out[['NNs']] = nn_idx
  out[['kNNs']] = knn_idx_dt
  out[['edge_weights']] = edge_weights_dt
  out[['all_dists']] = all_dists
  out[["dist_comp"]] = dist_comp
  out[["dist_thresh"]] = dist_cutoff
  return(out)
}

#' Visualize as velocity informed force directed graph
#'
#' @export
#'
graphViz = function(observed, projected, k, distance_metric = "L2", similarity_metric = "cosine", distance_weight = 1, distance_threshold = 1, similarity_threshold = -1, weighted = TRUE, remove_unconnected = TRUE, cell.colors = NA, title = NA, plot = TRUE, return_graph = FALSE){
  #observed, projected, k, distance_metric, similarity_metric, similarity_threshold: same arguments needed for projected neighbors
  #cell.colors: list of length nCells with colors corresponding to cluster IDs
  #return_graph: logical indicating whether to return graph object g and fdg coordinates fdg
  #
  
  ncells = ncol(observed)
  if (is.na(title)){
    title = ""
  }
  
  #find projected neighbors 
  nns = projectedNeighbors(observed, projected, k, distance_metric, similarity_metric, distance_weight, distance_threshold, similarity_threshold)
  print("Done finding neighbors")
  
  #make edge list 
  edgeList = matrix(nrow = 0, ncol = 2)
  edgeWeights = c()
  for (n in seq(1:k)){
    edgeList = rbind(edgeList, cbind(seq(1,ncells),nns$kNNs[,n]))
    edgeWeights = c(edgeWeights, nns$edge_weights[,n])
  }
  edgeList = na.omit(edgeList)
  edgeWeights = na.omit(edgeWeights)
  
  #make graph 
  #initialize empty graph with all cells 
  g = make_empty_graph(n=ncells,directed = TRUE)
  #add edges defined in edgeList
  edgeList = as.vector(t(edgeList)) #changing to required format for add_edges
  g = add_edges(g,edges = edgeList)
  #g = graph_from_edgelist(edgeList,directed = TRUE)    #old edgeList format
  #add edge weights if specified 
  if (weighted){
    print("calculating weights")
    E(g)$weight = max(edgeWeights) - edgeWeights
    #E(g)$weight = abs(edgeWeights)
    #print(abs(edgeWeights)[1:10])
  }
  #add vertex colors corresponding to cluster
  V(g)$color = cell.colors
  V(g)$size = 2
  E(g)$arrow.size = 0.5
  vertex.names = colnames(observed)
  #V(g)$name = vertex.names
  
  if (gsize(g)==0){
    print("WARNING: graph has no edges. Try lowering the similarity threshold.")
  }
  
  if (remove_unconnected){
    unconnected.vertices = which(degree(g)==0)
    g = delete.vertices(g,unconnected.vertices)
    if (length(unconnected.vertices)>0){
      vertex.names = vertex.names[-unconnected.vertices]
    }
    
  }
  
  #make force directed graph 
  fdg = layout_with_fr(g,dim=2)
  colnames(fdg) = c("C1","C2")
  rownames(fdg) = vertex.names
  print("Done making graph")
  
  if (plot){
    #plot both graphs 
    #par(mfrow = c(1,2))
    #plot(g)
    #plot.igraph(g,layout = fdg,vertex.label = NA, vertex.size = 5, vertex.color = adjustcolor(col = V(g)$color, alpha.f = 0.1), edge.color = "black") #####
    #plot(scale(fdg), col = cell.colors, pch = 16, main = paste("FDG cell coordinates: \n", title))
    plot(scale(fdg), col = V(g)$color, pch = 16, main = paste(title), cex = 0.5)
    
    #plot velocity on FDG embedding 
    # show.velocity.on.embedding.cor(scale(fdg),vel, n=100, scale='sqrt', cell.colors=cell.colors,cex=1, arrow.scale=1,
    #                                show.grid.flow=TRUE, min.grid.cell.mass=0.5, grid.n=30, arrow.lwd=1, main = paste("FDG embedding: ",title))
    #text(scale(fdg)+0.1,labels = seq(1,dim(fdg)[1]))
  }
  
  if (return_graph){
    out = list()
    out[['graph']] = g
    out[['fdg_coords']] = fdg
    out[['projected_neighbors']] = nns
    return(out)
  }
  
}

#' Consistency score calculation
#'
#' @export
#'
consistency = function(fdg.coords,delta.exp,nNeighbors,plot.hist = FALSE){
  #fdg.coords: coordinates of FDG in which consistency score is to be evaluated
  #delta.exp: change in gene expression based on velocity (part of velocyto output)
  #nNeighbors: number of nearest neighbors in FDG to find and calculate correlation with


  ncells = nrow(fdg.coords)
  cell.names = row.names(fdg.coords)
  neighbors = nn2(fdg.coords,k=nNeighbors+1)$nn.idx[,2:(nNeighbors+1)] # nNeighbors nearest neighbors to each cell: nCells x nNeighbors
  #neighbor.cors = t(sapply(seq(1:ncells), function(i) sapply(seq(1:nNeighbors), function(n) cor(delta.exp[,i],delta.exp[,neighbors[i,n]]))))  #### convert to c++
  neighbor.cors = pwiseCors(as.matrix(delta.exp),neighbors,nNeighbors)
  rownames(neighbor.cors) = cell.names
  cell.consistency = rowMeans(neighbor.cors)

  if (plot.hist){
    hist(cell.consistency,breaks = 100)
  }
  return(cell.consistency)
}
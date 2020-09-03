myDist = function(cell_i,proj_i,nn_i,distance_metric="L2",similarity_metric="cosine"){
  #cell_i: current cell
  #proj_i: projected state of cell i based on velocity
  #nn_i: putative nearest neighbor to cell_i 
  #distance_metric: "L1" or "L2" 
  #similarity_metric" "cosine" or "pearson" 
  
  #distance to minimize is between proj_i and nn_i: d = nn_i - proj_i 
  #angle to minimize is between velocity (v = proj_i - cell_i) and cell_i --> nn_i (n = nn_i - cell_i)
  
  d = nn_i - proj_i
  v = proj_i - cell_i 
  n = nn_i - cell_i 
  
  #distance  
  if (distance_metric=="L2"){
    invDist = 1/stats::dist(rbind(nn_i,proj_i),method = "euclidean")
  } else if (distance_metric=="L1"){
    invDist = 1/stats::dist(rbind(nn_i,proj_i),method = "manhattan")
  }
  
  #vector similarity 
  if (similarity_metric=="cosine"){
    negVectSim = -1*(v %*% n)/(sqrt(sum(v^2))*sqrt(sum(n^2))) 
  } else if (similarity_metric=="pearson"){
    negVectSim = 1 - cor(v,n, method = "pearson")
  }
  
  dists = list()
  dists[["invDist"]] = invDist
  dists[["negVectSim"]] = negVectSim
  dists[["myDist"]] = invDist*negVectSim
  return(dists)
}

projectedNeighbors = function(observed,projected,k,distance_metric="L2",similarity_metric="cosine",similarity_threshold = -1){
  #observed: genes x cells matrix of observed cells 
  #projected: genes x cells matrix of projected states of cells in observed (same order)
  #k: number of nearest neighbors 
  #distance_metric: "L1" or "L2" 
  #similarity_metric: "cosine" or "pearson". NOTE" pearson similarity behaves weird with two dimensions 
  #similarity_threshold: minimum similarity between velocity vector and cell_i-->nn_i vector for edge to be included. default all included. 
  
  observed = t(observed)
  projected = t(projected)
  n = nrow(observed)
  
  nn_idx = matrix(NA,nrow=nrow(observed), ncol = 1) #projected nearest neighbor 
  knn_idx = matrix(NA,nrow=nrow(observed), ncol = k) #projected k nearest neighbors 
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
    
    #calculate distances between cell_i and all other cells j [USE SAPPLY HERE???]
    for (j in seq(1,nrow(obs_exc_i))){
      nn_j = obs_exc_i[j,] 
      dist_j = myDist(cell_i,proj_i,nn_j,distance_metric,similarity_metric) #distance between cell_i and current cell_j
      nn_dists_i[j] = dist_j$myDist
      nn_invDist_i[j] = dist_j$invDist
      nn_negVectSim_i[j] = dist_j$negVectSim
    }
    
    
    #add cell_i's dists to all_dists 
    all_dists[i,-i] = nn_dists_i
    all_invDist[i,-i] = nn_invDist_i
    all_negVectSim[i,-i] = nn_negVectSim_i
    
    #find index of minimum distance between cell_i and other cells --> nn_i 
    min_dist_idx = which(nn_dists_i==min(nn_dists_i))
    #find indices of k nearest neighbors 
    k_min_dists_idx = order(nn_dists_i)[1:k]
    #exclude neighbor if similarity below threshold 
    sim_nn_k = nn_negVectSim_i[k_min_dists_idx] #similarities of minimum distance neighbors 
    new_k_idx = k_min_dists_idx[sim_nn_k <= (-1*similarity_threshold)] #indices of minimum distane neighbors excluding those whose similarities are above threshold 
    
    #correct min_idx
    min_dist_idx[which(min_dist_idx>=i)] = min_dist_idx[which(min_dist_idx>=i)] + 1
    #correct min k idx 
    k_min_dists_idx[which(k_min_dists_idx>=i)] = k_min_dists_idx[which(k_min_dists_idx>=i)] + 1
    new_k_idx[which(new_k_idx>=i)] = new_k_idx[which(new_k_idx>=i)] + 1
    
    #add current cell's neighbors to all cell neighbors matrix
    nn_idx[i] = min_dist_idx
    #knn_idx[i,] = k_min_dists_idx
    if (length(new_k_idx)>0){
      knn_idx[i,c(1:length(new_k_idx))] = new_k_idx
    }
    
  }
  dist_comp = list()
  dist_comp[["invDist"]] = all_invDist
  dist_comp[["negVectSim"]] = all_negVectSim
  out = list()
  out[['NNs']] = nn_idx
  out[['kNNs']] = knn_idx
  out[['all_dists']] = all_dists
  out[["dist_comp"]] = dist_comp
  return(out)
}

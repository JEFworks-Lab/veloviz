#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector myDist(NumericVector cell_i, NumericVector proj_i, NumericVector nn_i, String distance_metric, String similarity_metric, float distance_weight){
  //cell_i: current cell
  //proj_i: projected state of cell i based on velocity
  //nn_i: putative nearest neighbor to cell_i 
  //distance_metric: "L1" or "L2" 
  //similarity_metric" "cosine" or "pearson" 
  
  //distance to minimize is between proj_i and nn_i: d = nn_i - proj_i 
  //angle to minimize is between velocity (v = proj_i - cell_i) and cell_i --> nn_i (n = nn_i - cell_i)
  
  int dim = cell_i.size();
  NumericVector d = nn_i - proj_i;
  NumericVector v = proj_i - cell_i;
  NumericVector n = nn_i - cell_i;
  
  //inverse distance between nn_i and proj
  float dist = 0;
  if (distance_metric=="L2"){
    for(int i = 0; i<dim; ++i){
      dist = dist + pow(d[i],2);
    }
    dist = pow((1+(distance_weight*sqrt(dist))),-1);
  } else if (distance_metric=="L1"){
    for(int i = 0; i<dim; ++i){
      dist = dist + abs(d[i]);
    }
    dist = pow((1+(distance_weight*dist)),-1);
  }
  
  //negative similarity between velocity and (difference between cell_i and nn_i)
  float sim = 0;
  if (similarity_metric=="cosine"){
    float dot = 0;
    float sum_v_sqr = 0;
    float sum_n_sqr = 0;
    for(int i = 0; i<dim; ++i){
      dot = dot + v[i]*n[i];
      sum_v_sqr = sum_v_sqr + pow(v[i],2);
      sum_n_sqr = sum_n_sqr + pow(n[i],2);
    }
    sim = (-1 * dot)/(sqrt(sum_v_sqr)*sqrt(sum_n_sqr));
  } else if (similarity_metric=="pearson"){
    float sum_v = 0;
    float sum_v_sqr = 0;
    float sum_n = 0;
    float sum_n_sqr = 0;
    float sum_vn = 0;
    for(int i = 0; i<dim; ++i){
      sum_v = sum_v + v[i];
      sum_n = sum_n + n[i];
      
      sum_v_sqr = sum_v_sqr + pow(v[i],2);
      sum_n_sqr = sum_n_sqr + pow(n[i],2);
      
      sum_vn = sum_vn + (v[i]*n[i]);
    }
    float corr = (((dim*sum_vn) - (sum_v*sum_n))/(sqrt((dim*sum_v_sqr)-pow(sum_v,2))*(sqrt((dim*sum_n_sqr)-pow(sum_n,2)))));
    sim = 1-corr;
  }
  
  
  float compDist = dist*sim;
  NumericVector myDists = NumericVector::create(compDist, dist, sim);
  
  
  
  return myDists;
  
}


// [[Rcpp::export]]
NumericMatrix pwiseDists(NumericVector cell_i, NumericVector proj_i, NumericMatrix obs_exc_i, String distance_metric, String similarity_metric, float distance_weight){
  //calculates myDists (compDist, dist, sim) between cell_i and all other cells j (in obs_exc_i)
  
  int n_neighbors = obs_exc_i.nrow();
  NumericMatrix allDists(n_neighbors,3);
  
  for(int j = 0; j<n_neighbors; ++j){
    NumericVector curr_cell_j = obs_exc_i(j,_);
    NumericVector currDist = myDist(cell_i,proj_i,curr_cell_j,distance_metric,similarity_metric,distance_weight);
    allDists(j,_) = currDist;
    
  }
  
  colnames(allDists) = CharacterVector::create("CompositeDistance", "InverseDistance", "NegativeSimilarity");
  
  return allDists;
}


// [[Rcpp::export]]
NumericMatrix pwiseCors(NumericMatrix deltaExp, NumericMatrix neighborIdx, int nNeighbors){
  // calculates correlation between velocity vectors (from deltaExp) of each cell and it nNeighbors nearest neighbors on the FDG graph whose indices are given in neighborIdx
  // deltaExp: nGenes x nCells matrix of projected change in gene expression (from velocyto)
  // neighborIdx: nCells x nNeighbors matrix of indices of each cell's nearest neighbors on the FDG embedding
  // nNeighbors: number of nearest neighbors to calculate each cells correlation
  
  
  int nCells = deltaExp.ncol();
  int nGenes = deltaExp.nrow();
  NumericMatrix neighborCors(nCells, nNeighbors);
  
  for(int i=0; i<nCells; ++i){
    NumericVector velocity_i = deltaExp(_,i); //current cell's velocity vector
    
    for(int n=0; n<nNeighbors; ++n){
      int currentNeighborIdx = neighborIdx(i,n) - 1;
      NumericVector velocity_n = deltaExp(_,currentNeighborIdx); //current neighbor's velocity vector
      
      float sum_i = 0;
      float sum_i_sqr = 0;
      float sum_n = 0;
      float sum_n_sqr = 0;
      float sum_in = 0;
      
      for(int g = 0; g<nGenes; ++g){
        sum_i = sum_i + velocity_i[g];
        sum_n = sum_n + velocity_n[g];
        
        sum_i_sqr = sum_i_sqr + pow(velocity_i[g],2);
        sum_n_sqr = sum_n_sqr + pow(velocity_n[g],2);
        
        sum_in = sum_in + (velocity_i[g]*velocity_n[g]);
      }
      float corr = (((nGenes*sum_in) - (sum_i*sum_n))/(sqrt((nGenes*sum_i_sqr)-pow(sum_i,2))*(sqrt((nGenes*sum_n_sqr)-pow(sum_n,2)))));
      neighborCors(i,n) = corr;
    }
    
  }
  
  return neighborCors;
}



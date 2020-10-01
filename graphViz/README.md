# graphViz latest scripts and outputs
\
**projectedNeighbors.R**  \
contains functions for creating and visualizing velocity based embedding  \
`myDist`: composite distance function based on Euclidean distance and cosine similarity. Can also compute using L1 distance and correlation.  \
`projectedNeighbors`: finds k nearest neighbors using `myDist` given observed and projected states, similarity threshold. Outputs array where columns are identified nearest neighbors for each cell (row).  \
`graphViz`: creates graph based on projected neighbors identified by `projectedNeighbors` and finds fdg layout. Projects velocities onto fdg embedding.   \
`consistency`: calculates cell consistency score given an embedding and velocity vectors.  \
**graphVizC.cpp**  
`myDist`: c++ implementation of `myDist` in `projectedNeighbors.R` described above.  
`pwiseDists`: calculates distance between a cell and a matrix of other cells using `myDist`. Used by `projectedNeighbors` above.  
\
**graph_pancViz.Rmd**  \
Visualization of graph based velocity embedding using pancreas data from scVelo. Looks at effects of changing parameters: k, simThresh, L1 vs L2 distance, cosine similarity vs correlation.\
*Output in ./outputs/graph_pancViz.html*\
\
**graphall_u2osViz.Rmd**  \
Visualization of graph based velocity embedding using U2OS data. Looks at effects of changing parameters: k, simThresh, L1 vs L2 distance, cosine similarity vs correlation.\
*Output in ./outputs/graphall_u2osViz.html*\
\
**graph_neuroViz.Rmd**  
Visualization of graph based velocity embedding using dentate gyrus data from scVelo. Looks at effects of changing parameters: k, simThresh, L1 vs L2 distance, cosine similarity vs correlation.\
*Output in ./outputs/graph_neuroViz.html*  

**projectedNeighbor_vs_velocityGraph.Rmd**  
Comparison of embedding made using our method vs using velocity graph from scVelo, which is based on similarity matrix.  
*Output in ./outputs/projectedNeighbor_vs_velocityGraph.html*   

**addingEdgeWeights.Rmd**  
Adding edge weights to graph used to make force directed embedding and comparing consistency scores between embeddings made with and without edge weights.  
*Output in ./outputs/addingEdgeWeights.html*  

**disconnectedTrajectories_panc.Rmd**  
Looking at the effect on graph based visualization of removing intermediate cell types in the endocrine development trajectory.  
*Output in ./outputs/disconnectedTrajectories.html*  

**2020_09_RotationPresSims.Rmd**\
Toy data simulations showing rationale behind composite distance.

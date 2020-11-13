# graphViz latest scripts and outputs
\
**projectedNeighbors.R**  \
contains functions for creating and visualizing velocity based embedding  \
`myDist`: composite distance function based on Euclidean distance and cosine similarity. Can also compute using L1 distance and pearson correlation.  \
`projectedNeighbors`: finds k nearest neighbors using `myDist` given observed and projected states, similarity threshold. Outputs array where columns are identified nearest neighbors for each cell (row).  \
`graphViz`: creates graph based on projected neighbors identified by `projectedNeighbors` and finds fdg layout. Projects velocities onto fdg embedding.   \
`consistency`: calculates cell consistency score given an embedding and velocity vectors.  
**projectedNeighbors_weightedCD.R**  
Same as above but with weighted distance added.  
**projectedNeighborsR.R**  
Same as above but no C++ implementation.  

\
**graphVizC.cpp**  
`myDist`: c++ implementation of `myDist` in `projectedNeighbors.R` described above.  
`pwiseDists`: calculates distance between a cell and a matrix of other cells using `myDist`. Used by `projectedNeighbors` above.  
`pwiseCor`: calculates correlation between velocity vectors of cell and its nearest neighbors in the FDG embedding. Used by `consistency` above.  
**graphVizCD_weighted**  
same as above but with weighted distance added.  
\
**graph_pancViz.Rmd**  \
Visualization of graph based velocity embedding using pancreas data from scVelo. Looks at effects of changing parameters: k, simThresh, L1 vs L2 distance, cosine similarity vs correlation.\
*Outputs*  
./outputs/graph_pancViz_subset_unweighted.html: original visualization on subset of cells using unweighted graph.  
./outputs/graph_pancViz_allCells_unweighted.html: visualization on all cells using unweighted graph.  
./outputs/graph_pancViz_subset_weighted.html: visualization on all cells using weighted graph.  

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
*Outputs*
./outputs/disconnectedTrajectories_panc_subset_weighted.html  
./outputs/disconnectedTrajectories_panc_allCells_unweighted.html  

**simCycle**  
Looking at unwrapping effect of cycle simulation.  
*Outputs*  
./outputs/simCycle.html  

**simTraj.Rmd**  
Simulating branching trajectories and comparing velocity based embedding to other embeddings when plotting all cells and when removing intermediates.  
*Outputs*  
./outputs/simTraj.html  
*Figures*  
Two branch points:  
./figures/two_branch_all.svg  
Two branch points, missing intermediates (two different layouts):  
./figures/two_branch_noInt.svg  
./figures/two_branch_noInt2.svg  
Three branch points:  
./figures/three_branch_all.svg  
Three branch points, missing intermediates (two different layouts):  
./figures/three_branch_noInt.svg  
./figures/three_branch_noInt2.svg  

**projectedNeighbors_vs_cellrankGraph.Rmd**  
Comparing to graph made from CellRank for all data vs data with missing intermediates.  
*Outputs*
./outputs/projectedNeighbors_vs_cellrankGraph_missing_fev.html  
./outputs/projectedNeighbors_vs_cellrankGraph_missing_ngnhigh.html  
./outputs/projectedNeighbors_vs_cellrankGraph_missing_somebeta.html

**testingWeightedCD.Rmd**  
Trying different distance weights on cycle simulation.  
*Outputs*  
./outputs/testingWeightedCD.html  

**testingThresholdDist.Rmd**  
Adding percentile distance threshold to projectedNeighbors and testing on cycle simulation  
*Outputs*  
./outputs/testingThresholdDist_cyclesim.html  

**2020_09_RotationPresSims.Rmd**\
Toy data simulations showing rationale behind composite distance.

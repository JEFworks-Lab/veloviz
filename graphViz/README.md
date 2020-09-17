# graphViz latest scripts and outputs \
\
**projectedNeighbors.R**  \
contains functions for creating and visualizing velocity based embedding  \
*myDist:* composite distance function based on Euclidean distance and cosine similarity. Can also compute using L1 distance and correlation.  \
*projectedNeighbors:* finds k nearest neighbors using composite distance function given observed and projected states, similarity threshold. Outputs array where columns are identified nearest neighbors for each cell (row).  \
*graphViz:* creates graph based on projected neighbors identified by `projectedNeighbors` and finds fdg layout. Projects velocities onto fdg embedding.   \
*consistency:* calculates cell consistency score given an embedding and velocity vectors.  \
\
**graph_pancViz.Rmd**  \
Visualization of graph based velocity embedding using pancreas data from scVelo. Looks at effects of changing parameters: k, simThresh, L1 vs L2 distance, cosine similarity vs correlation.\
*Output in ./outputs/graph_pancViz.html*\
\
**graphall_u2osViz.Rmd**  \
Visualization of graph based velocity embedding using U2OS data. Looks at effects of changing parameters: k, simThresh, L1 vs L2 distance, cosine similarity vs correlation.\
*Output in ./outputs/graphall_u2osViz.html*\
\
**2020_09_RotationPresSims.Rmd**\
Toy data simulations showing rationale behind composite distance.

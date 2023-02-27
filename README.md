# SparseGTCCANet

The Sparse Generalized Tensor Canonical Correlation Network Inference (SGTCCA-Net) is a novel multi-omics network analysis pipeline that is designed for extract higher-order/pair-wise correlation exists in the multi-omics data. Below is the workflow of this novel network analysis pipeline:


![Figures](https://user-images.githubusercontent.com/36642788/221464373-dcc8d03d-b410-484c-9d0f-1d3c89028d20.png)

In general, this network analysis pipeline can be partitioned into multiple steps:

- Molecular feature weights extraction with Sparse Generalized Tensor Canonical Correlation Analysis.
- Global network construction based on molecular feature weights.
- Network prunning with PageRank algorithm and NetSHy summarization score.
- Network edge-cut and visualization.

Step 1 and 2 is similar to SmCCNet workflow, but step 3 above is novel for SGTCCA-Net pipeline. It is a stepwise approach with a combination of different evaluation metrics. Assume in global network there are $p$ molecular features, then this approach is given by:

- 1. Calculate PageRank score for all $p$ molecular features, and rank them accordingly.
- 2. Starting from a minimally possible network size $m_1$ (defined by users), include the top $m_1$ molecular features, calculate NetSHy summarization scores.
- 3. Iterate the following steps until reaching the maxmimally possible network size $m_2$ (defined by users):
  - Include one more molecular feature into the network based on PageRank score, then calculate NetSHy score for this updated network.
  - Calculate the correlation between this network summarization score and phenotype, call it $\rho$. 
- 4. Find minimal network size $k$ with $\rho$ within the 80% range (this value can be defined by user) of the maximum correlation w.r.t. phenotype.
- 5. For each network size that is greater than $k$ and within the 80% rankge of the maximum correlation w.r.t. phenotype, calculate the correlation between summarization score at each network size and $k$.
- 6. Based on step 5, setting up and correlation threshold (defined by user), such that only candidate network size with correlation (calculated in step 5) that pass the threshold will be kept.
- 7. Based on the final network size candidates, select the maxmimum/minimum of them depending on user's preference.

After network prunning, the network is densely connected, which is not desirable since some molecular features may not connect with each other. Therefore, in the last  step (before network visualization), we use the actual correlation between molecular features to filter out edges with two nodes that are actually weakly correlated, this cut-off can be set up by user (usually I would recommend setting up this threshold within the range of 0.1-0.2).

For network visualization, the recommended way to do it is through Cytoscape (software can be downloaded here: https://cytoscape.org/). To allow the communication between Cytocape and R, RCy3 package can be used (https://www.bioconductor.org/packages/release/bioc/html/RCy3.html).

All the source code for SGTCCA-Net has been included under the code folder along with an example script. Most of the parameters are used in simulation studies and don't need to be changed. The algorithm is totally parallizable and this feature has been included in the code. 
  

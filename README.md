# SparseGTCCANet

The Sparse Generalized Tensor Canonical Correlation Network Inference (SGTCCA-Net) is a novel multi-omics network analysis pipeline that is designed for extract higher-order/pair-wise correlation exists in the multi-omics data. Below is the workflow of this novel network analysis pipeline:


![Figures](https://user-images.githubusercontent.com/36642788/221464373-dcc8d03d-b410-484c-9d0f-1d3c89028d20.png)

In general, this network analysis pipeline can be partitioned into multiple steps:

- Molecular feature weights extraction with Sparse Generalized Tensor Canonical Correlation Analysis.
- Global network construction based on molecular feature weights.
- Network prunning with PageRank algorithm and NetSHy summarization score.
- Network edge-cut and visualization.

Step 1 and 2 is similar to SmCCNet workflow, but step 3 is novel for SGTCCA-Net pipeline. It is a stepwise approach with a combination of different evaluation metrics. Assume in global network there are $p$ molecular features, then this approach is given by:

- 1. Calculate PageRank score for all $p$ molecular features, and rank them accordingly.
- 2. Starting from a minimally possible network size $m_1$ (defined by users), include the top $m_1$ molecular features, calculate NetSHy summarization scores.
- 3. Iterate the following steps until reaching the maxmimally possible network size $m_2$ (defined by users):
  - Include one more molecular feature into the network based on PageRank score, then calculate NetSHy score for this updated network.
  - Calculate the correlation between this network summarization score and (1) phenotype and (2) summarization score at minimal network size.
- 4. 
  

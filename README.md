# SparseGTCCANet

The Sparse Generalized Tensor Canonical Correlation Network Inference (SGTCCA-Net) is a novel multi-omics network analysis pipeline that is designed for extract higher-order/pair-wise correlation exists in the multi-omics data. Below is the workflow of this novel network analysis pipeline:


![Figures](https://user-images.githubusercontent.com/36642788/221464373-dcc8d03d-b410-484c-9d0f-1d3c89028d20.png)

In general, this network analysis pipeline can be partitioned into multiple steps:

- Molecular feature weights extraction with Sparse Generalized Tensor Canonical Correlation Analysis.
- Global network construction based on molecular feature weights.
- Network prunning with PageRank algorithm and NetSHy summarization score.

Step 1 and 2 is similar to SmCCNet workflow, but step 3 is novel for SGTCCA-Net pipeline. It is a stepwise approach with a combination of different evaluation metrics. Assume in global network there are $p$ molecular features, then this approach is given by:

- Calculate PageRank score for all $p$ molecular features, and rank them accordingly.
- Starting from a minimally possible network size $m_1$ (defined by users), 

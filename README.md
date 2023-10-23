# Kernel Eigenfunction Embeddings for High-Dimensional Data


We propose a kernel-spectral embedding algorithm for learning low-dimensional nonlinear structures from high-dimensional and noisy observations, where the datasets are assumed to be sampled from an intrinsically low-dimensional manifold and corrupted by high-dimensional noise. The algorithm employs an adaptive bandwidth selection procedure which does not rely on prior knowledge of the underlying manifold. The obtained low-dimensional embeddings can be further utilized for downstream purposes such as data visualization, clustering and prediction. Our method is theoretically justified and practically interpretable. We establish the convergence of the final embeddings to their noiseless counterparts when the dimension and size of the samples are comparably large, and characterize the effect of the signal-to-noise ratio on the rate of convergence and phase transition. We also prove convergence of the embeddings to the eigenfunctions of an integral operator defined by the kernel map of some reproducing kernel Hilbert space capturing the underlying nonlinear structures. Numerical simulations and analysis of three real datasets show the superior empirical performance of the proposed method, compared to many existing methods, on learning various manifolds in diverse applications.

The method is based on the paper:

Ding, X., and Ma, R. (2023) Learning Low-Dimensional Nonlinear Structures from High-Dimensional Noisy Data: An Integral Operator Approach. **Annals of Statistics** 51(4), 1744-1769. https://arxiv.org/pdf/2203.00126.


# Content

The folder /data contains the example datasets analyzed in the paper.

The folder /code contains R scripts for analyzing the example datasets.

# System Requirements

The meta-visualization package requires only a standard computer with enough RAM to support the operations defined by a user. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB
CPU: 4+ cores, 3.3+ GHz/core

The R implementation of the method is tested under R version 4.1.1, and requires the R packages: `rARPACK`,`Rfast`,`rARPACK`,`dimRed`,`cluster`,`ggplot2`,`lle`,`Seurat`,`TSCAN`,`SCORPIUS`.

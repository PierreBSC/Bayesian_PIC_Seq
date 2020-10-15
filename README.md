# Bayesian PIC-Seq

This repository contains a set of Stan-based scripts allowing to analyze PIC-seq (Physically Interacting Cell - Sequencing) data. Those scripts have been tested on the data from the Nature Biotech paper "Dissecting cellular crosstalk by sequencing physically interacting cells".

Briefly this PIC-seq processing script take as an input the gene expression profile of PIC/doublets together with the mean expression profiled of each single-cell cluster. The fitting relies on [Stan](https://mc-stan.org/) to create an extremely efficient Hamiltonian Monte Carlo (HMC) sampler able to estimate the contribution of each single-cell cluster to the gene expression profile of each PIC/doublet.

Those scripts are relying on several R packages : **rstan**, **Ternary** and **pheatmap**.

Analysis of the conventionnal single-cell data can be performed by any package of choice such as [Seurat](https://satijalab.org/seurat/), [Pagoda2](https://github.com/kharchenkolab/pagoda2) or [MetaCell](https://tanaylab.github.io/metacell/).

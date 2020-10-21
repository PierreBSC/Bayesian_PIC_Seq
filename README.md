 :# Bayesian PIC-Seq

This repository contains a set of Stan-based scripts allowing to analyze PIC-seq (Physically Interacting Cell - Sequencing) data. Those scripts have been tested on the data from the Nature Biotech paper "Dissecting cellular crosstalk by sequencing physically interacting cells".

Briefly this PIC-seq processing script take as an input the gene expression profile of PIC/doublets together with the mean expression profiled of each single-cell cluster. The fitting relies on [Stan](https://mc-stan.org/) to create an extremely efficient Hamiltonian Monte Carlo (HMC) sampler able to estimate the contribution of each single-cell cluster to the gene expression profile of each PIC/doublet.

Those scripts are relying on several R packages : **rstan**, **Matrix**, **Ternary**, **foreach** and **doParallel**.

Analysis of the conventionnal single-cell data can be performed by any package of choice such as [Seurat](https://satijalab.org/seurat/), [Pagoda2](https://github.com/kharchenkolab/pagoda2) or [MetaCell](https://tanaylab.github.io/metacell/).

All the core functions are stored in the **Bayesian_PIC_seq.R** script.

# Example of use :

Here we will use the data provided by the authors of the paper "Dissecting cellular crosstalk by sequencing physically interacting cells" : they are fully available on the GEO website (GEO accession number : GSE135382). We focused on the in-vitro data where mouse monocytes derived Dendritic Cells (moDCs) where put in contact with T-cells. PIC and single-cells were sequenced 3h, 20h and 48h post co-culture. 

The analysis of the single-cells is done using the Pagoda2 [Pagoda2](https://github.com/kharchenkolab/pagoda2) pipeline. Data are first filtered to remove low quality cells/PIC and non expressed genes 

```r
##data_raw object : Sparse Matrix describing the gene expression of PIC and single cells 

par(las=1,bty="l")
lib_size = colSums(data_raw)

hist(log10(lib_size),n=50,col="grey")
abline(v=log10(350),lty=2,col="red")
boxplot(log10(lib_size)~Gating,outline=F,col=c("red3","green3","orange3"),xlab="Gating",ylab="Total UMIs (Log10)")

gene_size = rowSums(data_raw)
hist(log10(gene_size),n=30,col="grey",xlim=c(0,5))
abline(v=log10(50),lty=2,col="red")

# data_count : filtered expression data, only cells with more than 500 UMIs and genes with more than 50 UMIs are kepts
data_count = data_raw[gene_size>50,lib_size>500]

#Filtering the Time and Experimental condition vectors
Gating_count = Gating[lib_size>500]
Time_point_count = Time_point[lib_size>500]

```

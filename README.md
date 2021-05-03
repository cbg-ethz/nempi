# Introduction

If many genes are perturbed in a population of cells, this can lead to
diseases like cancer. The perturbations can happen in different ways,
e.g. via mutations, copy number abberations or methylation. However,
not all perturbations are observed in all samples.

Nested Effects Model-based perturbation inference (NEM$\pi$) uses
observed perturbation profiles and gene expression data to infer
unobserved perturbations and augment observed ones. The causal 
network of the perturbed genes
(P-genes) is modelled as an adjacency matrix $\phi$ and the genes with
observed gene expression (E-genes) are modelled with the attachment
$\theta$ with $\theta_{ij}=1$, if E-gene $j$ is attached to 
S-gene $i$. If E-gene $j$ is attached to P-gene $i$, $j$ shows an effect
for a perturbation of P-gene $i$. Hence, $\phi\theta$ predicts gene 
expression profiles, which can be compared to the real 
data. NEM$\pi$ iteratively infers a network $\phi$ based on 
gene expression profiles and a perturbation profile, and the 
perturbation profile based on a network $\phi$.

Install:
--------

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("nempi")
```

Most recent (devel) version:

```r
install.packages("devtools")

library(devtools)

install_github("cbg-ethz/nempi")

library(nempi)
```
For the reproduction of the publication see the script in the other directory.

## Reference

Pirkl M, Beerenwinkel N (2021). "Inferring perturbation profiles of cancer samples." Bioinformatics. https://doi.org/10.1093/bioinformatics/btab113.
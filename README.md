# CANMOD v2
CanMod is a computational pipeline to identify gene regulatory modules containing 6 main steps:

  - Step 1: Cluster the differentially expressed (DE) genes based on GO-term similarity and generate Gene Clusters (GCs)
  - Step 2: Identify regulators (i.e., miRNAs or transcription factors) for each DE gene
  - Step 3: Cluster the regulators based on shared target similarity and generate Regulator Clusters (RCs)
  - Step 4: Generate modules using the regulators and the target GCs
  
  Iteratively refining modules running:
  - Step 5: Refine each module using biclustering
  - Step 6: Merge modules that share high similarity

Run CANMOD with the command `Rscript CANMOD_v2.R` for sample run, or `Rscript CANMOD_v2.R *data_name*`  where *data_name.rda* is the file name containing required data under `data` folder. *data_name*.rda will have the following variables:

  - `mRNA`: A data frame containing gene expression data (Rows are all DE genes and columns are samples).
  - `miRNA`: A data frame containing MicroRNA expression data (Rows are all DE microRNAs and columns are samples).
  - `methyl`: A data frame containing DNA methylation data (Rows are all DE genes and columns are samples).
  - `cnv`: A data frame containing Copy number aberration data (Rows are all DE genes and columns are samples).
! Samples in these four variables will be in the same order.

---

### User Options
- Hyperparameters (lines 34-36):
  - `drop.thr`: Stopping criteria for iterating Step 5 and 6. (default is 0.05)
  - `GC_thr`: Gene cluster (GC) threshold (default is 0.45)
  - `RC_thr`: Regulator cluster (RC) threshold (default is 0.10)
- To input user-defined GCs instead of CANMOD-generated GCs, provide the following files:
  - `user.defined.gene.cluster.list.rda`: File containing the list of GCs
  - `user.defined.similarities.rda`: File containing the similarity matrix between genes in GCs

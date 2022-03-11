# CANMOD v2
CanMod is a computational pipeline to identify gene regulatory modules containing 6 main steps:

  - Step 1: Cluster the differentially expressed (DE) genes genes based on their GO-term similarity
  - Step 2: Identify regulators (i.e., miRNAs or transcription factors (TFs)) for each DE gene
  - Step 3: Cluster the regulators based on shared target similarity
  - Step 4: Generate modules using the regulators and the target gene clusters
  Iteratively refining modules running:
  - Step 5: Refine each modules so that the expression of regulators and target genes in each module is correlated
  - Step 6: Merge modules which share a large amount of similar targets

Source code to run CANMOD could be found in CANMOD_v2.R.

Run CANMOD with the command `Rscript CANMOD_v2.R *data_name*`  where *data_name* is one of [BRCA, KIRC, LUSC, PRAD, UCEC] for the corresponding preprocessed TCGA datasets of CANMOD.

---

### User Options
- Hyperparameters (lines 34-36):
  - `drop.thr`: Stopping criteria for iterating Step 5 and 6. (default is 0.05)
  - `GC_thr`: Gene cluster threshold (default is 0.45)
  - `RC_thr`: Regulator cluster threshold (default is 0.10)
- User-defined Gene Clusters for user defined data:
  - `user.defined.gene.cluster.list.rda`: File containing the list of Gene Clusters
  - `user.defined.similarities.rda`: File containing the similarity matrix between genes in Gene Clusters

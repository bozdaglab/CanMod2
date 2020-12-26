require(clusterProfiler)
require(biomaRt)

# Enrichment of inferred modules
# inputs should be entrez gene id
GO_enrichment_ofModule = function(module, genes){
  result <- enrichGO(module, universe = as.character(genes), OrgDb = 'org.Hs.eg.db', ont="BP", qvalueCutoff = 0.01)
  return(result$Description)
}

# inputs should be entrez gene id
KEGG_enrichment_ofModule = function(module, genes){
  result <- enrichKEGG(module, universe = as.character(genes), qvalueCutoff = 0.01, organism = "hsa", pAdjustMethod = "BH")
  return(result$Description)
}

# inputs should be entrez gene id
CH_enrichment_ofModule = function(module, genes){
  t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, entrez_gene)
  em <- enricher(gene = module, TERM2GENE=t2g, universe = genes, qvalueCutoff = 0.01)
  return(em$Description)
}
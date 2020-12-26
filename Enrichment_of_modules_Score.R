install.packages("msigdbr")
library("biomaRt")
library("clusterProfiler")
library("org.Hs.eg.db")
library("msigdbr")
# How to use the function

  # Generate universe as entrez ONCE using following three lines
  universe.gene = rownames(variable.selection.data$mRNA)
  ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  universe = unlist(unname(getBM(attributes = c('entrezgene_id'), filters = 'hgnc_symbol', values=universe.gene, mart=ensembl)))
  
  All_GO_enrichment_Module = function(module, genes){
    resulth <- enrichGO(module, universe = as.character(genes), OrgDb = 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.99, qvalueCutoff = 0.99)
    return (resulth@result$ID)
  }
  
  # For each module, run following two lines
  go.enrich.module.all<-lapply(1:10,function(i){
    module<-dtc.final.target.list[[i]]
    module = unlist(unname(getBM(attributes = c('entrezgene_id'), filters = 'hgnc_symbol', values=module, mart=ensembl)))
    return (All_GO_enrichment_Module(module, universe))
    
  })
  all.go.enrich.terms<-unique(unlist(go.enrich.module.all))
  
  Selected_GO_enrichment_Module = function(module, genes){
    resultq <- enrichGO(module, universe = as.character(genes), OrgDb = 'org.Hs.eg.db', ont="BP",qvalueCutoff = 0.01)
    common.go.terms<-intersect(all.go.enrich.terms,resultq@result$ID)
    r<-resultq@result
    go.es<-mean(-log(r[common.go.terms,6]))
    go.es.median<-median(-log(r[common.go.terms,6]))
    go.es.only.in.module<-mean(-log(r[,6]))
    go.es.median.only.in.module<-median(-log(r[,6]))
    return (list(mean=go.es,median=go.es.median,onlyinmodulemean=go.es,onlyinmodulemedian=go.es.median))
  }
  
  selected.go.enrich.module<-lapply(1:10,function(i){
    module<-dtc.final.target.list[[i]]
    module = unlist(unname(getBM(attributes = c('entrezgene_id'), filters = 'hgnc_symbol', values=module, mart=ensembl)))
    return (Selected_GO_enrichment_Module(module, universe))
    
  })

module.go.es.mean.list<-lapply(1:10, function(x) {
  return (selected.go.enrich.module[[x]][["mean"]])
})

module.go.es.mean<-mean(unlist(module.go.es.mean.list))

module.go.es.median.list<-lapply(1:10, function(x) {
  return (selected.go.enrich.module[[x]][["median"]])
})
module.go.es.mean.median<-mean(unlist(module.go.es.median.list))


module.onlyingo.es.mean.list<-lapply(1:10, function(x) {
  return (selected.go.enrich.module[[x]][["onlyinmodulemean"]])
})

module.onlyingo.es.mean<-mean(unlist(module.onlyingo.es.mean.list))

module.onlyingo.es.median.list<-lapply(1:10, function(x) {
  return (selected.go.enrich.module[[x]][["onlyinmodulemedian"]])
})
module.onlyingo.es.mean.median<-mean(unlist(module.onlyingo.es.median.list))

##################################KEGG############################################################

# inputs should be entrez gene id
All_KEGG_enrichment_Module = function(module, genes){
  
  resulth<- enrichKEGG(module, universe = as.character(genes), pvalueCutoff = 0.99, organism = "hsa", pAdjustMethod = "BH",qvalueCutoff = 0.99)
  return (resulth@result$ID)
}


kegg.enrich.module.all<-lapply(1:10,function(i){
  module<-dtc.final.target.list[[i]]
  module = unlist(unname(getBM(attributes = c('entrezgene_id'), filters = 'hgnc_symbol', values=module, mart=ensembl)))
  return (All_KEGG_enrichment_Module(module, universe))
  
})

all.kegg.enrich.terms<-unique(unlist(kegg.enrich.module.all))

Selected_KEGG_enrichment_Module = function(module, genes){
  resultq<- enrichKEGG(module, universe = as.character(genes),organism = "hsa", pAdjustMethod = "BH",qvalueCutoff = 0.01)
  common.kegg.terms<-intersect(all.kegg.enrich.terms,resultq@result$ID)
  r<-resultq@result
  kegg.es<-mean(-log(r[common.kegg.terms,6]))
  kegg.es.median<-median(-log(r[common.kegg.terms,6]))
  kegg.es.only.in.module<-mean(-log(r[,6]))
  kegg.es.median.only.in.module<-median(-log(r[,6]))
  return (list(mean=kegg.es,median=kegg.es.median,onlyinmodulemean=kegg.es,onlyinmodulemedian=kegg.es.median))
}


selected.kegg.enrich.module<-lapply(1:10,function(i){
  module<-dtc.final.target.list[[i]]
  module = unlist(unname(getBM(attributes = c('entrezgene_id'), filters = 'hgnc_symbol', values=module, mart=ensembl)))
  return (Selected_KEGG_enrichment_Module(module, universe))
  
})

module.kegg.es.mean.list<-lapply(1:10, function(x) {
  return (selected.kegg.enrich.module[[x]][["mean"]])
})

module.kegg.es.mean<-mean(unlist(module.kegg.es.mean.list))

module.kegg.es.median.list<-lapply(1:10, function(x) {
  return (selected.kegg.enrich.module[[x]][["median"]])
})
module.kegg.es.mean.median<-mean(unlist(module.kegg.es.median.list))


module.onlyinkegg.es.mean.list<-lapply(1:10, function(x) {
  return (selected.kegg.enrich.module[[x]][["onlyinmodulemean"]])
})

module.onlyinkegg.es.mean<-mean(unlist(module.onlyinkegg.es.mean.list))

module.onlyinkegg.es.median.list<-lapply(1:10, function(x) {
  return (selected.kegg.enrich.module[[x]][["onlyinmodulemedian"]])
})


##############################CH###################################################


# inputs should be entrez gene id
All_CH_enrichment_Module = function(module, genes){
  
  t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, entrez_gene)
  emh <- enricher(gene=module, TERM2GENE=t2g, universe = as.character(genes),pvalueCutoff = 0.99, qvalueCutoff = 0.99)
  return (emh@result$ID)
}


ch.enrich.module.all<-lapply(1:10,function(i){
  module<-dtc.final.target.list[[i]]
  module = unlist(unname(getBM(attributes = c('entrezgene_id'), filters = 'hgnc_symbol', values=module, mart=ensembl)))
  return (All_CH_enrichment_Module(module, universe))
  
})

all.ch.enrich.terms<-unique(unlist(ch.enrich.module.all))

Selected_CH_enrichment_Module = function(module, genes){
  t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, entrez_gene)
  emq <- enricher(gene=module, TERM2GENE=t2g, universe = as.character(genes),qvalueCutoff = 0.01)
  common.ch.terms<-intersect(all.ch.enrich.terms,emq@result$ID)
  r<-emq@result
  ch.es<-mean(-log(r[common.ch.terms,6]))
  ch.es.median<-median(-log(r[common.ch.terms,6]))
  ch.es.only.in.module<-mean(-log(r[,6]))
  ch.es.median.only.in.module<-median(-log(r[,6]))
  return (list(mean=ch.es,median=ch.es.median,onlyinmodulemean=ch.es,onlyinmodulemedian=ch.es.median))
}


selected.ch.enrich.module<-lapply(1:10,function(i){
  module<-dtc.final.target.list[[i]]
  module = unlist(unname(getBM(attributes = c('entrezgene_id'), filters = 'hgnc_symbol', values=module, mart=ensembl)))
  return (Selected_CH_enrichment_Module(module, universe))
  
})

module.ch.es.mean.list<-lapply(1:10, function(x) {
  return (selected.ch.enrich.module[[x]][["mean"]])
})

module.ch.es.mean<-mean(unlist(module.ch.es.mean.list))

module.ch.es.median.list<-lapply(1:10, function(x) {
  return (selected.ch.enrich.module[[x]][["median"]])
})
module.ch.es.mean.median<-mean(unlist(module.ch.es.median.list))


module.onlyinch.es.mean.list<-lapply(1:10, function(x) {
  return (selected.ch.enrich.module[[x]][["onlyinmodulemean"]])
})

module.onlyinch.es.mean<-mean(unlist(module.onlyinch.es.mean.list))

module.onlyinch.es.median.list<-lapply(1:10, function(x) {
  return (selected.ch.enrich.module[[x]][["onlyinmodulemedian"]])
})
module.onlyinch.es.mean.median<-mean(unlist(module.onlyinch.es.median.list))
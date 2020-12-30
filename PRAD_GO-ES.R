install.packages("msigdbr")
library("biomaRt")
library("clusterProfiler")
library("org.Hs.eg.db")
library("msigdbr")

final.target.module.prad<-lapply(1:length(PRAD_modules), function(i){
  
  target.module<-PRAD_modules[[i]][["targets"]]
  return(target.module)
  
})

universe.gene = rownames(regression.data$mRNA)
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
universe = unlist(unname(getBM(attributes = c('entrezgene_id'), filters = 'hgnc_symbol', values=universe.gene, mart=ensembl, useCache = FALSE)))

GO.ES.Score<-function(module.list, gene.names){
  
  GO_enrichment_ofModule = function(module, genes){
    result <- enrichGO(module, universe = as.character(genes), OrgDb = 'org.Hs.eg.db', ont="BP", qvalueCutoff = 0.99, pvalueCutoff = 0.99)
    return(result$ID)
  }
  
  go.enrich.module.all<-lapply(1:length(module.list),function(i){
    module<-module.list[[i]]
    module = unlist(unname(getBM(attributes = c('entrezgene_id'), filters = 'hgnc_symbol', values=module, mart=ensembl, useCache = FALSE)))
    return (GO_enrichment_ofModule(module, universe))
    
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
  
  selected.go.enrich.module<-lapply(1:length(module.list),function(i){
    module<-module.list[[i]]
    module = unlist(unname(getBM(attributes = c('entrezgene_id'), filters = 'hgnc_symbol', values=module, mart=ensembl)))
    return (Selected_GO_enrichment_Module(module, universe))
    
  })
  
  module.go.es.mean.list<-lapply(1:length(module.list), function(x) {
    return (selected.go.enrich.module[[x]][["mean"]])
  })
  
  module.go.es.mean<-mean(unlist(module.go.es.mean.list))
  cat("GO ES score is :",module.go.es.mean)
  
  module.go.es.median.list<-lapply(1:length(module.list), function(x) {
    return (selected.go.enrich.module[[x]][["median"]])
  })
  module.go.es.mean.median<-mean(unlist(module.go.es.median.list))
  
  cat("median GO ES score is :",module.go.es.mean.median)
  module.onlyingo.es.mean.list<-lapply(1:length(module.list), function(x) {
    return (selected.go.enrich.module[[x]][["onlyinmodulemean"]])
  })
  
  module.onlyingo.es.mean<-mean(unlist(module.onlyingo.es.mean.list))
  cat("GO ES score is :",module.onlyingo.es.mean)
  module.onlyingo.es.median.list<-lapply(1:length(module.list), function(x) {
    return (selected.go.enrich.module[[x]][["onlyinmodulemedian"]])
  })
  module.onlyingo.es.mean.median<-mean(unlist(module.onlyingo.es.median.list))
  cat("median GO ES score is :",module.onlyingo.es.mean.median)
  return(list(module.go.es.mean, module.go.es.mean.median,module.onlyingo.es.mean,module.onlyingo.es.mean.median))
  
}

cluster = parallel::makeCluster(detectCores()-5)
parallel::clusterExport(cl = cluster,
                        varlist = c("final.target.module.prad","de.genes"),
                        envir = environment())
invisible(parallel::clusterEvalQ(cl = cluster,expr = library(org.Hs.eg.db)))
go.es.prad= pbapply::pblapply(cl = cluster,FUN = GO.ES.Score(final.target.module.prad,de.genes))

stopCluster(cluster)

save(go.es.prad, file="go.es.prad.rda")

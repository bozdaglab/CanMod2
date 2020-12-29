require(clusterProfiler)
require(biomaRt)

######################PAIRWISE CORRELATION OF TARGET GENES OF MODULES######################

# inputs should be gene modules as list, gene names and gene expression of de genes
pairwise.corr.target.genes.modules<-function(module.list, gene.names, expression.cor){
  sig<-0
  pairwise.cor.target.gene<-lapply(1:length(module.list), function(i){
    
    final.target.gene<-module.list[[i]]
    cor_mat<-abs(expression.cor[final.target.gene,final.target.gene])
    diag(cor_mat) <- 0
    low_cor_mat <- cor_mat
    low_cor_mat[upper.tri(low_cor_mat)] <- 0
    cor.final.target.gene<-mean(low_cor_mat[lower.tri(low_cor_mat)])
    f<-0
    cor.random.target.gene<-vector()
    for (j in 1:1000){
      random.target.gene<-sample(x=gene.names,size=length(final.target.gene),replace= FALSE, prob= NULL)
      cor_mat<-abs(expression.cor[random.target.gene,random.target.gene])
      diag(cor_mat) <- 0
      low_cor_mat <- cor_mat
      low_cor_mat[upper.tri(low_cor_mat)] <- 0
      cor.random.target.gene[j]<-mean(low_cor_mat[lower.tri(low_cor_mat)])
      if(cor.final.target.gene<cor.random.target.gene[j]){
        f<-f+1
      }
    }
    p.value<-f/1000
    if(p.value<.05){
      sig<-sig+1
    }
    
    return(list(cor.final.target.gene,sig,cor.random.target.gene))
    
  })
  
  #x<-sum(sapply(pairwise.cor.target.gene, function(module){module[[2]] == 1}))
  
  num.significant.modules<-sum(sapply(pairwise.cor.target.gene, function(module){module[[2]] == 1}))
  num.not.significant.modules<-length(dtc.final.module.list)-num.significant.modules
  per.significant.modules<-(num.significant.modules/length(dtc.final.module.list))*100
  per.not.significant.modules<-(num.not.significant.modules/length(dtc.final.module.list))*100
  
  cat("Percentage of significant modules is: ", per.significant.modules)
  cat("Percentage of not significant modules is: ",per.not.significant.modules)
  
}


#############################ENRICHMENT OF MODULES#######################################

# In order to calculate GO.ES.score, KEGG.ES.score, CH.ES.score, inputs should be gene modules as list, gene names
install.packages("msigdbr")
library("biomaRt")
library("clusterProfiler")
library("org.Hs.eg.db")
library("msigdbr")


GO_enrichment_ofModule = function(module, genes){
  result <- enrichGO(module, universe = as.character(genes), OrgDb = 'org.Hs.eg.db', ont="BP", qvalueCutoff = 0.99, pvalueCutoff = 0.99)
  return(result$ID)
}
GO.ES.Score<-function(module.list, gene.names){
  # Generate universe as entrez ONCE using following three lines
  universe.gene = gene.names
  ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  universe = unlist(unname(getBM(attributes = c('entrezgene_id'), filters = 'hgnc_symbol', values=universe.gene, mart=ensembl)))
  
  # All_GO_enrichment_Module = function(module, genes){
  #   resulth <- enrichGO(module, universe = as.character(genes), OrgDb = 'org.Hs.eg.db', ont="BP", pvalueCutoff = 0.99, qvalueCutoff = 0.99)
  #   return (resulth@result$ID)
  # }
  
  go.enrich.module.all<-lapply(1:length(module.list),function(i){
    module<-module.list[[i]]
    module = unlist(unname(getBM(attributes = c('entrezgene_id'), filters = 'hgnc_symbol', values=module, mart=ensembl)))
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
}


##################################KEGG############################################################
KEGG_enrichment_ofModule = function(module, genes){
  result <- enrichKEGG(module, universe = as.character(genes), qvalueCutoff = 0.99, organism = "hsa", pAdjustMethod = "BH",pvalueCutoff = 0.99)
  return(result$ID)
}

KEGG.ES.Score<- function(module.list, gene.names){
  universe.gene = gene.names
  ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  universe = unlist(unname(getBM(attributes = c('entrezgene_id'), filters = 'hgnc_symbol', values=universe.gene, mart=ensembl)))
  
  # inputs should be entrez gene id
  # All_KEGG_enrichment_Module = function(module, genes){
  #   
  #   resulth<- enrichKEGG(module, universe = as.character(genes), pvalueCutoff = 0.99, organism = "hsa", pAdjustMethod = "BH",qvalueCutoff = 0.99)
  #   return (resulth@result$ID)
  # }
  
  
  kegg.enrich.module.all<-lapply(1:length(module.list),function(i){
    module<-module.list[[i]]
    module = unlist(unname(getBM(attributes = c('entrezgene_id'), filters = 'hgnc_symbol', values=module, mart=ensembl)))
    return (KEGG_enrichment_ofModule(module, universe))
    
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
  
  
  selected.kegg.enrich.module<-lapply(1:length(module.list),function(i){
    module<-module.list[[i]]
    module = unlist(unname(getBM(attributes = c('entrezgene_id'), filters = 'hgnc_symbol', values=module, mart=ensembl)))
    return (Selected_KEGG_enrichment_Module(module, universe))
    
  })
  
  module.kegg.es.mean.list<-lapply(1:length(module.list), function(x) {
    return (selected.kegg.enrich.module[[x]][["mean"]])
  })
  
  module.kegg.es.mean<-mean(unlist(module.kegg.es.mean.list))
  cat("KEGG ES score is :",module.kegg.es.mean)
  module.kegg.es.median.list<-lapply(1:length(module.list), function(x) {
    return (selected.kegg.enrich.module[[x]][["median"]])
  })
  module.kegg.es.mean.median<-mean(unlist(module.kegg.es.median.list))
  
  cat("median KEGG ES score is :",module.kegg.es.mean.median)
  module.onlyinkegg.es.mean.list<-lapply(1:length(module.list), function(x) {
    return (selected.kegg.enrich.module[[x]][["onlyinmodulemean"]])
  })
  
  module.onlyinkegg.es.mean<-mean(unlist(module.onlyinkegg.es.mean.list))
  cat("KEGG ES score is :",module.onlyinkegg.es.mean)
  module.onlyinkegg.es.median.list<-lapply(1:length(module.list), function(x) {
    return (selected.kegg.enrich.module[[x]][["onlyinmodulemedian"]])
  })
  module.onlyinkegg.es.mean.median<-mean(unlist(module.onlyinkegg.es.mean.median.list))
  cat("median KEGG ES score is :",module.onlyinkegg.es.mean.median)
}



##############################CH###################################################
CH_enrichment_ofModule = function(module, genes){
  t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
    dplyr::select(gs_name, entrez_gene)
  em <- enricher(gene = module, TERM2GENE=t2g, universe = genes, qvalueCutoff = 0.99, pvalueCutoff= 0.99)
  return(em$ID)
}

CH.ES.Score<-function(module.list, gene.names){
  universe.gene = gene.names
  ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
  universe = unlist(unname(getBM(attributes = c('entrezgene_id'), filters = 'hgnc_symbol', values=universe.gene, mart=ensembl)))

  # inputs should be entrez gene id
  # All_CH_enrichment_Module = function(module, genes){
  #   
  #   t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  #     dplyr::select(gs_name, entrez_gene)
  #   emh <- enricher(gene=module, TERM2GENE=t2g, universe = as.character(genes),pvalueCutoff = 0.99, qvalueCutoff = 0.99)
  #   return (emh@result$ID)
  # }
  
  
  ch.enrich.module.all<-lapply(1:length(module.list),function(i){
    module<-module.list[[i]]
    module = unlist(unname(getBM(attributes = c('entrezgene_id'), filters = 'hgnc_symbol', values=module, mart=ensembl)))
    return (CH_enrichment_ofModule(module, universe))
    
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
  
  
  selected.ch.enrich.module<-lapply(1:length(module.list),function(i){
    module<-module.list[[i]]
    module = unlist(unname(getBM(attributes = c('entrezgene_id'), filters = 'hgnc_symbol', values=module, mart=ensembl)))
    return (Selected_CH_enrichment_Module(module, universe))
    
  })
  
  module.ch.es.mean.list<-lapply(1:length(module.list), function(x) {
    return (selected.ch.enrich.module[[x]][["mean"]])
  })
  
  module.ch.es.mean<-mean(unlist(module.ch.es.mean.list))
  cat("CH ES score is :",module.ch.es.mean)
  module.ch.es.median.list<-lapply(1:length(module.list), function(x) {
    return (selected.ch.enrich.module[[x]][["median"]])
  })
  module.ch.es.mean.median<-mean(unlist(module.ch.es.median.list))
  
  cat("median CH ES score is :",module.ch.es.mean.median)
  module.onlyinch.es.mean.list<-lapply(1:length(module.list), function(x) {
    return (selected.ch.enrich.module[[x]][["onlyinmodulemean"]])
  })
  
  module.onlyinch.es.mean<-mean(unlist(module.onlyinch.es.mean.list))
  cat("CH ES score is :",module.onlyinch.es.mean)
  module.onlyinch.es.median.list<-lapply(1:length(module.list), function(x) {
    return (selected.ch.enrich.module[[x]][["onlyinmodulemedian"]])
  })
  module.onlyinch.es.mean.median<-mean(unlist(module.onlyinch.es.median.list))
  cat("median CH ES score is :",module.onlyinch.es.mean.median)
  
}

########################################################################################

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
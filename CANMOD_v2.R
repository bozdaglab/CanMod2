# Load utility functions to run CanMod ----------------------------------------------------
source("data/helper_functions.R")

### --- Global variables ---- #####
args = commandArgs(trailingOnly = T)
if (!exists("cancer.type")){
  cancer.type = ifelse(length(args) !=0 , args[1], "sample")
}

# Load required input for CanMod  ---------------------------------------------------------
load("data/miRNA_mRNA_interactions.rda")
colnames(putative.miRNA.mRNA.interactions) = c("miRNA","target")
load("data/TF_target_interactions.rda")
colnames(TF.target) = c("tf", "target")
load(paste0("data/", cancer.type, ".rda"))
regression.data = list(mRNA = mRNA, miRNA = miRNA, methyl = methyl, cna = cnv); rm(mRNA, miRNA, methyl, cnv)
{
  samples = colnames(regression.data$mRNA)[substr(colnames(regression.data$mRNA), 14, 15) %in% "01"]
  regression.data$mRNA = regression.data$mRNA[,samples]
  regression.data$miRNA = regression.data$miRNA[,samples]
  regression.data$methyl = regression.data$methyl[,samples]
  regression.data$cna = regression.data$cna[,samples]
  regression.data$miRNA.target.interactions = putative.miRNA.mRNA.interactions
  regression.data$miRNA.target.interactions = regression.data$miRNA.target.interactions[regression.data$miRNA.target.interactions$miRNA %in% rownames(regression.data$miRNA),]
  regression.data$miRNA.target.interactions = regression.data$miRNA.target.interactions[regression.data$miRNA.target.interactions$target %in% rownames(regression.data$mRNA),]
  regression.data$tf.target.interactions = TF.target
  regression.data$tf.target.interactions = regression.data$tf.target.interactions[regression.data$tf.target.interactions$tf %in% rownames(regression.data$mRNA),]
  regression.data$tf.target.interactions = regression.data$tf.target.interactions[regression.data$tf.target.interactions$target %in% rownames(regression.data$mRNA),]
  rownames(regression.data$miRNA) = gsub(rownames(regression.data$miRNA), pattern = "[-]", replacement = ".")
  regression.data$miRNA.target.interactions$miRNA = gsub(regression.data$miRNA.target.interactions$miRNA, pattern = "[-]", replacement = ".")
}

# Adjust hyperparameters ---------------------------------------------------------
drop_thr = 0.05
GC_thr = 0.45
RC_thr = 0.10

# STEP 1: Get GO-based cluster --------------------------------------------------------------------------
de.genes = rownames(regression.data$mRNA)
save(de.genes, file =  paste0("results/", cancer.type, "_Step1_de.genes.rda"))

if(file.exists(paste0(cancer.type, "_user.defined.gene.cluster.list.rda")) && file.exists(paste0(cancer.type, "_user.defined.similarities.rda"))){
  cat("CANMOD is using user-defined Gene Clusters and similarities.")
  load(paste0(cancer.type, "_user.defined.gene.cluster.list.rda"))
  load(paste0(cancer.type, "_user.defined.similarities.rda"))
}else{
  hsGO2 = godata('org.Hs.eg.db', keytype = "SYMBOL", ont="BP", computeIC=FALSE)
  
  cat("Start computing GO-based similarity matrix \n")
  time = proc.time()
  de.gene.bp.sim = mgeneSim(de.genes, semData=hsGO2, measure="Wang", combine="BMA", verbose=FALSE, drop = NULL)
  time = proc.time() - time
  cat(paste0("running time for computing gene similarity based on GO:", time[3], "\n"))
  
  cat(paste0("Done with obtaining similarity matrix in cluster_de_genes for", cancer.type, "\n")) 
  save(de.gene.bp.sim, file =  paste0("results/", cancer.type, "_Step1_de.gene.bp.sim.rda"))
  
  sm = reshape2::melt(de.gene.bp.sim)
  sm = sm[sm$Var1 != sm$Var2,]
  sm = sm[sm$value > GC_thr,]
  graph = graph_from_data_frame(sm[,c(1, 2)], directed = F) %>% set_edge_attr("weight", value = as.numeric(sm$value))
  gc.list = as.list(cluster_walktrap(graph))
  gene.cluster.list = lapply(1:length(gc.list), function(index){
    gc.list[[index]]
  })
  gene.cluster.df = list2df(gene.cluster.list)
  names(gene.cluster.df) = c("gene","cluster")
  save(gene.cluster.list, file =  paste0("results/", cancer.type, "_Step1_gene.cluster.list.rda"))
}
common.de.genes <- unique(unlist(gene.cluster.list))

# Update regression data
{
  regression.data$exp = regression.data$mRNA
  regression.data$mRNA = regression.data$mRNA[rownames(regression.data$mRNA) %in% common.de.genes,]
  regression.data$methyl = regression.data$methyl[rownames(regression.data$methyl) %in% common.de.genes,]
  regression.data$cna = regression.data$cna[rownames(regression.data$cna) %in% common.de.genes,]
  regression.data$miRNA.target.interactions = regression.data$miRNA.target.interactions[regression.data$miRNA.target.interactions$target %in% rownames(regression.data$mRNA),]
  regression.data$tf.target.interactions = regression.data$tf.target.interactions[regression.data$tf.target.interactions$target %in% rownames(regression.data$mRNA),]
}

target.sim = de.gene.bp.sim
target.cluster.list = gene.cluster.list
names(target.cluster.list) = as.character(1:length(target.cluster.list))
target.cluster.df = list2df(target.cluster.list); names(target.cluster.df) = c("target","cluster")
save(target.cluster.list, file =  paste0("results/", cancer.type, "_Step1_target.cluster.list.rda"))
save(target.cluster.df, file =  paste0("results/", cancer.type, "_Step1_target.cluster.df.rda"))


# STEP 2: Variable selection -----------------------------------
mRNA.targets = rownames(regression.data$mRNA)
cluster = parallel::makeCluster(parallel::detectCores()-2)
parallel::clusterExport(cl = cluster, 
                        varlist = c("mRNA.targets","regression.data", "getDF"),
                        envir = environment())
invisible(parallel::clusterEvalQ(cl = cluster,expr = library(data.table, HDCI)))
coefs = pbapply::pblapply(cl = cluster, X = 1:length(mRNA.targets), FUN = function(mRNA.index){
  mRNA = mRNA.targets[mRNA.index]
  cat("Start", mRNA, "\n")
  include.lncRNA = grepl(mRNA, pattern = "lncRNA")
  regression.df = getDF(regression.data = regression.data, mRNA = mRNA, include.lncRNA = include.lncRNA)
  if (is.null(regression.df)){
    return(NULL)
  }
  y.matrix = as.matrix(regression.df$mRNA)
  x.matrix = as.matrix(regression.df[,-1])
  colnames(x.matrix) = colnames(regression.df)[-1]
  candidate.regulators = colnames(x.matrix)
  if(ncol(x.matrix) == 1){
    if(cor(x.matrix, y.matrix) > 0){
      return(NULL)
    }else{
      if (cor.test(x.matrix, y.matrix)$p.value >= 0.05){
        return(NULL)
      }else{
        normal.lasso.dt = data.frame(regulator=candidate.regulators, target=mRNA, 
                                     coef = cor(x.matrix, y.matrix), 
                                     stringsAsFactors = F)
        normal.lasso.dt = as.data.table(normal.lasso.dt)
        normal.lasso.dt$pair = paste(normal.lasso.dt$regulator, normal.lasso.dt$target,sep="~")
        return(normal.lasso.dt)
      }
    }
  }
  
  selected.coefs = lapply(1:100, function(iter){
    y.matrix = as.matrix(regression.df$mRNA)
    x.matrix = as.matrix(regression.df[,-1])
    colnames(x.matrix) = colnames(regression.df)[-1]
    
    lasso.result = tryCatch({
      HDCI::Lasso(x = x.matrix, y = y.matrix,
                  fix.lambda = F,nfolds = 10,
                  cv.method = "cv1se",
                  standardize = T, parallel = F)
    }, warning = function(w) {
      NULL
    }, error = function(e) {
      NULL
    })
    
    if (is.null(lasso.result)){
      return(NULL)
    }
    
    coefs = lasso.result$beta
    names(coefs) =  colnames(x.matrix) 
    coefs = coefs[which(coefs != 0)]
    return(coefs)
  })
  
  normal.lasso.dt = NULL
  dt.list = lapply(1:length(selected.coefs), function(target.index){
    regulators = selected.coefs[[target.index]]
    if (is.null(regulators)){
      return(NULL)
    }
    target = mRNA
    regulator.names = names(regulators)
    regulator.coef =unname(regulators)
    if (length(regulator.names) > 0){
      df = data.frame(regulator=regulator.names,
                      target = rep(target,length(regulator.names)),
                      coef=as.vector(regulator.coef),
                      stringsAsFactors = F)
      df$pair = paste(df$regulator,df$target,sep="~")
      dt = as.data.table(df)
      return(dt)
    }else{
      return(NULL)
    }
  })
  if (is.null(dt.list)){
    normal.lasso.dt = NULL
    return(normal.lasso.dt)
  }
  
  normal.lasso.dt = tryCatch({
    rbindlist(dt.list)
  }, warning = function(w) {
    NULL
  }, error = function(e) {
    NULL
  })
  
  cat("Done", mRNA, "\n")
  return(normal.lasso.dt)
})
names(coefs) = mRNA.targets

bt.interval.list = pbapply::pblapply(1:length(mRNA.targets),FUN = function(mRNA.index){
  if (mRNA.index %% 100 == 0) cat(mRNA.index, "\n")
  mRNA = mRNA.targets[mRNA.index]
  cat(mRNA,"\n")
  include.lncRNA = grepl(mRNA, pattern = "lncRNA")
  regression.df = getDF(regression.data = regression.data, mRNA = mRNA, include.lncRNA = F)
  
  if (is.null(regression.df)){
    return(NULL)
  }
  
  y.matrix = as.matrix(regression.df$mRNA)
  x.matrix = as.matrix(regression.df[,-1])
  colnames(x.matrix) = colnames(regression.df)[-1]
  
  if(ncol(x.matrix) == 1){
    if(cor(x.matrix, y.matrix) > 0){
      return(NULL)
    }else{
      if (cor.test(x.matrix, y.matrix)$p.value >= 0.05){
        return(NULL)
      }else{
        bt.interval = matrix(c(-1000,1000), ncol = 1)
        colnames(bt.interval) = colnames(x.matrix)
        return(bt.interval)
      }
    }
  }
  
  num.bt.replications = 100
  btlasso.result = tryCatch({
    bootLasso(x = x.matrix, y = y.matrix,
              B = num.bt.replications,
              cv.method = "cv1se",
              type.boot = "paired",
              standardize = T)
  }, warning = function(w) {
    NULL
  }, error = function(e) {
    NULL
  })
  
  if (is.null(btlasso.result)){
    return(NULL)
  }
  
  bt.interval = btlasso.result$interval
  colnames(bt.interval) = colnames(x.matrix)
  return(bt.interval)
})
names(bt.interval.list) = mRNA.targets
stopCluster(cluster)

regulator.list = lapply(1:length(coefs), function(index){
  gene.name = mRNA.targets[index]
  coef.dt = coefs[[gene.name]]
  bt.interval.dt = bt.interval.list[[gene.name]]
  list(coef.dt=coef.dt, bt.interval.dt=bt.interval.dt)
})
names(regulator.list) = mRNA.targets

save(coefs, file =  paste0("results/", cancer.type, "_Step2_coefs.rda"))
save(bt.interval.list, file =  paste0("results/", cancer.type, "_Step2_bt.interval.list.wt.rda"))
save(regulator.list, file =  paste0("results/", cancer.type, "_Step2_regulator.list.wt.rda"))

# STEP 3: Cluster regulators based on shared targets similarity -------------------------------------------------------------------
regulator.target.pair.list = lapply(1:length(regulator.list), function(gene.index) {
  gene.pair.info.list = regulator.list[[gene.index]]
  gene.coef.dt = gene.pair.info.list$coef.dt
  if (nrow(gene.coef.dt) == 0 || is.null(gene.coef.dt)) {
    return(NULL)
  }
  gene.coef.stat.dt = gene.coef.dt[, list(count = .N, 
                                          median.coef = median(.SD$coef)),
                                   by = "pair"]
  gene.coef.stat.dt$regulator = sapply(strsplit(gene.coef.stat.dt$pair, split = "[~]"), function(v)
    v[1])
  gene.coef.stat.dt$target = sapply(strsplit(gene.coef.stat.dt$pair, split = "[~]"), function(v)
    v[2])
  gene.coef.stat.dt = gene.coef.stat.dt[, c(4, 5, 1, 2, 3), ]
  gene.bt.interval.dt = tryCatch({
    gene.pair.info.list$bt.interval.dt
  }, warning = function(w) {
    NULL
  }, error = function(e) {
    NULL
  })
  
  if (is.null(gene.bt.interval.dt) || nrow(gene.bt.interval.dt) == 0){
    return(NULL)
  }
  
  lower.percentile = sapply(1:nrow(gene.coef.stat.dt), function(row.index) {
    regulator = gene.coef.stat.dt$regulator[row.index]
    return(gene.bt.interval.dt[1, regulator])
  })
  gene.coef.stat.dt$lower.percentile = lower.percentile
  
  upper.percentile = sapply(1:nrow(gene.coef.stat.dt), function(row.index) {
    regulator = gene.coef.stat.dt$regulator[row.index]
    return(gene.bt.interval.dt[2, regulator])
  })
  gene.coef.stat.dt$upper.percentile = upper.percentile
  
  gene.coef.stat.dt$confidence =  (gene.coef.stat.dt$median.coef >= lower.percentile) &
    (gene.coef.stat.dt$median.coef <= upper.percentile)
  
  return(gene.coef.stat.dt)
})
regulator.target.pair.dt = rbindlist(regulator.target.pair.list)
regulator.target.pair.dt = regulator.target.pair.dt[which(regulator.target.pair.dt$confidence == T & regulator.target.pair.dt$count >=75)]
regulator.target.pair.dt = regulator.target.pair.dt[-which(regulator.target.pair.dt$regulator %in% c("CNA","Methyl"))]
to.removed.indicies = which(grepl(regulator.target.pair.dt$regulator, pattern = "hsa") & regulator.target.pair.dt$median.coef >=0 )

if(length(to.removed.indicies) != 0){
  regulator.target.pair.dt = regulator.target.pair.dt[-to.removed.indicies]
}

lasso.df = regulator.target.pair.dt[,c(1,2)]
save(lasso.df, file =  paste0("results/", cancer.type, "_Step3_lasso.df.rda"))

regulation.target.df = matrix(data = 0, 
                              nrow = length(unique(lasso.df$target)),
                              ncol = length(unique(lasso.df$regulator)),
                              dimnames = list(unique(lasso.df$target),
                                              unique(lasso.df$regulator)))
regulation.target.df = as.data.frame(regulation.target.df)
for (row.index in 1:nrow(regulation.target.df)){
  if (row.index %% 100 == 0) print(row.index)
  target = rownames(regulation.target.df)[row.index]
  regulators = lasso.df$regulator[which(lasso.df$target == target)]
  regulation.target.df[target,regulators] = 1
}
dim(regulation.target.df)

# Cluster regulators
{
  sm = reshape2::melt(1-as.matrix(dist(t(regulation.target.df), method = "binary")))
  sm = sm[sm$Var1 != sm$Var2,]
  sm = sm[sm$value > RC_thr,]
  graph = graph_from_data_frame(sm[,c(1, 2)], directed = F) %>% set_edge_attr("weight", value = as.numeric(sm$value))
  rc.list = as.list(cluster_walktrap(graph, weights =  E(graph)$weight))
  regulator.cluster.list = lapply(1:length(rc.list), function(index){
    rc.list[[index]]
  })
  regulator.cluster.df = list2df(regulator.cluster.list)
  names(regulator.cluster.df) = c("regulator","cluster")
}
save(regulator.cluster.list, file =  paste0("results/", cancer.type, "_Step3_regulator.cluster.list.rda"))
save(regulator.cluster.df, file =  paste0("results/", cancer.type, "_Step3_regulator.cluster.df.rda"))


# STEP 4: Generate candidate modules ------------------------------------------------------------------------
regulator.target.cluster.list = lapply(1:length(regulator.cluster.list), function(index){
  regulators = regulator.cluster.list[[index]]
  regulator.target.dt = lasso.df[lasso.df$regulator %in% regulators]
  regulator.target.dt = regulator.target.dt[regulator.target.dt$target %in% unique(target.cluster.df$target)]
  regulator.target.dt$cluster = as.numeric(target.cluster.df$cluster[match(regulator.target.dt$target, target.cluster.df$target)])
  return(regulator.target.dt)
}) 
save(regulator.target.cluster.list, file =  paste0("results/", cancer.type, "_Step4_regulator.target.cluster.list.rda"))

seed.target.list = lapply(1:length(regulator.target.cluster.list), function(k){
  regulator.cluster =  regulator.target.cluster.list[[k]]
  target.clusters = unique(regulator.cluster$cluster)
  target.clusters = lapply(1:length(target.clusters), function(i){
    target.cluster.index = target.clusters[[i]]
    seed.targets = unique(regulator.cluster$target[regulator.cluster$cluster == target.cluster.index] )
    all.targets.in.cluster = unique(target.cluster.df$target[target.cluster.df$cluster == target.cluster.index]  )
    all.targets.in.cluster = setdiff(all.targets.in.cluster, seed.targets)
    return(list(seed.targets= seed.targets, all.targets.in.cluster= all.targets.in.cluster))
  })
  return(target.clusters)
})
save(seed.target.list, file =  paste0("results/", cancer.type, "_Step4_seed.target.list.rda"))

expression.df<-rbind(regression.data$mRNA,regression.data$miRNA)
expression.cor <- cor(as.data.frame(t(regression.data$mRNA)))
cor.threshold = quantile(abs(expression.cor), 0.9)
expression.cor <- cor(as.data.frame(t(expression.df)))
cor.threshold.reg = quantile(abs(expression.cor), 0.9)

selected.seed.target.list =  lapply(1:length(seed.target.list), function(index){
  overall.seed.list =  seed.target.list[[index]]
  seed.regulators = regulator.cluster.list[[index]]
  seed.partner.list = lapply(1:length(overall.seed.list), function(i){
    seed.list = overall.seed.list[[i]]
    seed.targets = seed.list$seed.targets
    seed.partners = seed.list$all.targets.in.cluster
    if(length(seed.partners) == 0){
      return(unique(seed.targets))  
    }
    seed.target.partner.list = lapply(1:length(seed.partners), function(k){
      seed.partner = seed.partners[k]
      seed.target.cor = expression.cor[seed.partner, seed.targets]
      seed.reg.cor = expression.cor[seed.partner, seed.regulators]
      if(length(seed.targets) == 1 && ((abs(seed.target.cor) > cor.threshold) ||  (abs(seed.reg.cor) > cor.threshold.reg)) ){
        return(seed.partner)
      }
      else if(length(seed.targets) > 1 && ((length(names(which(abs(seed.target.cor) > cor.threshold))) >= 1) || 
                                           (length(names(which(abs(seed.reg.cor) > cor.threshold.reg))) >= 1)) ){
        return(seed.partner)
      }
      return(NULL)
    })
    seed.target.partner.list = unique(c(seed.targets, unlist(seed.target.partner.list)))
    return(seed.target.partner.list)
  })
})
save(selected.seed.target.list, file =  paste0("results/", cancer.type, "_Step4_selected.seed.target.list.rda"))

regulator.target.list =  lapply(1:length(regulator.cluster.list), function(index){
  regulators = regulator.cluster.list[[index]]
  target.list = selected.seed.target.list[[index]]
  regulator.target.list = lapply(target.list, function(targets){
    list(regulators=regulators, targets=targets)
  })
  return(regulator.target.list)
})
save(regulator.target.list, file =  paste0("results/", cancer.type, "_Step4_regulator.target.list.rda"))

simplified.regulator.target.list = list()
for (i in 1:length(regulator.target.list)){
  for (j in 1:length(regulator.target.list[[i]])){
    simplified.regulator.target.list = rlist::list.append(simplified.regulator.target.list, 
                                                          regulator.target.list[[i]][[j]])
  }
}
save(simplified.regulator.target.list, file =  paste0("results/", cancer.type, "_Step4_simplified.regulator.target.list.rda"))


# STEPS 5 and 6 iteratively -------------------------------------------------------------------
it.no = 1; drop = 1
while(drop >=drop_thr || it.no <= 2){
  cat(paste0('------------------- Iteration ',it.no,' -----------------\n'))
  # Step 5: Refine module using biclustering ----------------------------------------------------------------------------------------------
  if(it.no > 1){
    temp.module.list = final.module.list
  }else{
    temp.module.list = simplified.regulator.target.list
  }
  x = sapply (temp.module.list, function(r) r$regulators)
  prev_reg_size = length(unique(unlist(x)))
  y = sapply (temp.module.list, function(r) r$targets)
  prev_target_size = length(unique(unlist(y))) 
  
  module.expression.list = lapply(1:length(temp.module.list), function(index){
    module = temp.module.list[[index]]
    target.df = expression.df[module$targets,, drop = F]
    regulator.df = expression.df[module$regulators,, drop = F]
    module = list(target.df=target.df,
                  regulator.df=regulator.df)
  })
  
  {
    options(warn = -1)
    {
      means = vector()
      for(index in 1:length(module.expression.list)){
        module <- module.expression.list[[index]]
        exp.values = expression.cor[rownames(module$target.df), rownames(module$regulator.df), drop = F]
        exp = (as.matrix(abs(exp.values)))*100
        my.delta=summary(as.vector(exp))[5]
        my.alpha=1
        my.number=sum(ncol(exp),nrow(exp))
        if(ncol(exp)>1 && nrow(exp)>1){
          testCluster <-biclust(x = exp, method=BCCC(), delta = my.delta, alpha=my.alpha, number=my.number) 
          for(x in 1:ncol(testCluster@RowxNumber)){
            if(is.na(summary(as.vector(exp[rownames(exp)[testCluster@RowxNumber[,x]], colnames(exp)[testCluster@NumberxCol[x,]]]))[3])){
              means = c(means, summary(as.vector(exp[rownames(exp), colnames(exp)]))[4])
            }else{
              means = c(means, summary(as.vector(exp[rownames(exp)[testCluster@RowxNumber[,x]], colnames(exp)[testCluster@NumberxCol[x,]]]))[4])
            }
          }
        }else{
          means = c(means, summary(as.vector(exp[rownames(exp), colnames(exp)]))[4])
        }
      }
      thr = unname(summary(means)[2])
      my.modules = temp.module.list
      for(index in 1:length(module.expression.list)){
        module <- module.expression.list[[index]]
        exp.values = expression.cor[rownames(module$target.df), rownames(module$regulator.df)]
        exp = (as.matrix(abs(exp.values)))*100
        my.delta=summary(as.vector(exp))[5]
        my.alpha=1
        my.number=sum(ncol(exp),nrow(exp))
        if(ncol(exp)>1 && nrow(exp)>1){
          testCluster <- biclust(x = exp, method=BCCC(), delta = my.delta, alpha=my.alpha, number=my.number)
          cluster.toremove = vector()
          for(x in 1:ncol(testCluster@RowxNumber)){
            if(is.na(summary(as.vector(exp[rownames(exp)[testCluster@RowxNumber[,x]], colnames(exp)[testCluster@NumberxCol[x,]]]))[3])){
            }else{
              if(unname(summary(as.vector(exp[rownames(exp)[testCluster@RowxNumber[,x]], colnames(exp)[testCluster@NumberxCol[x,]]]))[4]) < thr){
                cluster.toremove = c(cluster.toremove, x)
              }
            }
          }
          
          if(length(cluster.toremove) == ncol(testCluster@RowxNumber)){
            my.modules[[index]]$regulators = NA
            my.modules[[index]]$targets = NA
          }else{
            cluster.tokeep = setdiff(1:ncol(testCluster@RowxNumber), cluster.toremove)
            if(ncol(testCluster@RowxNumber) == 1){
              selection.tar = t(testCluster@RowxNumber)
              selection.reg = t(testCluster@NumberxCol)
            }else{
              selection.tar = t(testCluster@RowxNumber)
              selection.reg = testCluster@NumberxCol
            }
            my.modules[[index]]$regulators = colnames(exp)[colSums(selection.reg[cluster.tokeep,,drop = F])>0]
            my.modules[[index]]$targets = rownames(exp)[colSums(selection.tar[cluster.tokeep,,drop = F])>0]
          }
        }else{
        }
      }
    }
    my.modules1 = my.modules; rm(my.modules)
    
    {
      means = vector()
      for(index in 1:length(module.expression.list)){
        module <- module.expression.list[[index]]
        exp.values = expression.cor[rownames(module$target.df), rownames(module$regulator.df), drop = F]
        exp = (t(as.matrix(abs(exp.values))))*100
        my.delta=summary(as.vector(exp))[5]
        my.alpha=1
        my.number=sum(ncol(exp),nrow(exp))
        if(ncol(exp)>1 && nrow(exp)>1){
          testCluster <- biclust(x = exp, method=BCCC(), delta = my.delta, alpha=my.alpha, number=my.number)
          for(x in 1:ncol(testCluster@RowxNumber)){
            if(is.na(summary(as.vector(exp[rownames(exp)[testCluster@RowxNumber[,x]], colnames(exp)[testCluster@NumberxCol[x,]]]))[3])){
              means = c(means, summary(as.vector(exp[rownames(exp), colnames(exp)]))[4])
            }else{
              means = c(means, summary(as.vector(exp[rownames(exp)[testCluster@RowxNumber[,x]], colnames(exp)[testCluster@NumberxCol[x,]]]))[4])
            }
          }
        }else{
          means = c(means, summary(as.vector(exp[rownames(exp), colnames(exp)]))[4])
        }
      }
      
      thr = unname(summary(means)[2])
      my.modules = temp.module.list
      for(index in 1:length(module.expression.list)){
        module <- module.expression.list[[index]]
        exp.values = expression.cor[rownames(module$target.df), rownames(module$regulator.df)]
        exp = (t(as.matrix(abs(exp.values))))*100
        my.delta=summary(as.vector(exp))[5]
        my.alpha=1
        my.number=sum(ncol(exp),nrow(exp))
        
        if(ncol(exp)>1 && nrow(exp)>1){
          testCluster <- biclust(x = exp, method=BCCC(), delta = my.delta, alpha=my.alpha, number=my.number) 
          
          cluster.toremove = vector()
          
          for(x in 1:ncol(testCluster@RowxNumber)){
            if(is.na(summary(as.vector(exp[rownames(exp)[testCluster@RowxNumber[,x]], colnames(exp)[testCluster@NumberxCol[x,]]]))[3])){
            }else{
              if(unname(summary(as.vector(exp[rownames(exp)[testCluster@RowxNumber[,x]], colnames(exp)[testCluster@NumberxCol[x,]]]))[4]) < thr){
                cluster.toremove = c(cluster.toremove, x)
              }
            }
          }
          
          if(length(cluster.toremove) == ncol(testCluster@RowxNumber)){
            my.modules[[index]]$regulators = NA
            my.modules[[index]]$targets = NA
          }else{
            cluster.tokeep = setdiff(1:ncol(testCluster@RowxNumber), cluster.toremove)
            if(ncol(testCluster@RowxNumber) == 1){
              selection.reg = t(testCluster@RowxNumber)
              selection.tar = t(testCluster@NumberxCol)
            }else{
              selection.reg = t(testCluster@RowxNumber)
              selection.tar = testCluster@NumberxCol
            }
            my.modules[[index]]$targets = colnames(exp)[colSums(selection.tar[cluster.tokeep,,drop = F])>0]
            my.modules[[index]]$regulators = rownames(exp)[colSums(selection.reg[cluster.tokeep,,drop = F])>0]
            
          } 
        }else{
        }
      }
    }
    my.modules2 = my.modules; rm(my.modules)
    options(warn = 0)
    
    my.modules = mapply(function(x, y) {list(list(selected.regulators = setdiff(union(x$regulators, y$regulators), NA), 
                                                  selected.targets = setdiff(union(x$targets, y$targets), NA)))}, my.modules1, my.modules2)
    
    
    toremove = unlist(sapply(my.modules, function(m) (length(m$selected.regulators) == 0 || length(m$selected.targets) == 0)||is.na(m$selected.regulators)[1] || is.na(m$selected.targets)[1]))
    temp.module.list = my.modules[!toremove]
    
    module.expression.list = lapply(1:length(temp.module.list), function(index){
      module = temp.module.list[[index]]
      target.df = expression.df[module$selected.targets,, drop = F]
      regulator.df = expression.df[module$selected.regulators,, drop = F]
      module = list(target.df=target.df,
                    regulator.df=regulator.df)
    })
    
    filtered.bic.module.list = temp.module.list
    filtered.bic.target.list = lapply(filtered.bic.module.list, function(module){
      module$selected.targets
    })
    filtered.bic.target.df = qdapTools::list2df(filtered.bic.target.list)
    colnames(filtered.bic.target.df) = c("target","module")
    
    filtered.bic.regulator.list = lapply(filtered.bic.module.list, function(module){
      module$selected.regulators
    })
    filtered.bic.regulator.df = qdapTools::list2df(filtered.bic.regulator.list)
    colnames(filtered.bic.regulator.df) = c("regulator","module")
    
    target.module.df = matrix(data = 0, 
                              nrow = length(unique(filtered.bic.target.df$target))+length(unique(filtered.bic.regulator.df$regulator)),
                              ncol = length(unique(filtered.bic.target.df$module)),
                              dimnames = list(c(unique(filtered.bic.target.df$target),unique(filtered.bic.regulator.df$regulator)),
                                              unique(filtered.bic.target.df$module)))
    
    target.module.df = as.data.frame(target.module.df)
    for (row.index in 1:nrow(target.module.df)){
      x = rownames(target.module.df)[row.index]
      modules = as.numeric(unique(c(filtered.bic.target.df$module[which(filtered.bic.target.df$target == x)],
                                    filtered.bic.regulator.df$module [which(filtered.bic.regulator.df$regulator == x)])))
      target.module.df[x,modules]= 1
    }
  }
  save(filtered.bic.module.list, file = paste0("results/", cancer.type, "_Step5_filtered.bic.module.list.", it.no, ".rda"))
  
  # STEP 6 -------------------------------------------------------------------
  {
    sm = reshape2::melt(1-as.matrix(dist(t(target.module.df), method = "binary")))
    sm = sm[sm$value >=0.8,]
    graph = graph_from_data_frame(sm[,c(1, 2)], directed = F) %>%
      set_edge_attr("weight", value = as.numeric(sm$value))
    module.cluster.list = as.list(cluster_walktrap(graph, weights =  E(graph)$weight))
  }
  final.module.list = lapply(1:length(module.cluster.list), function(index){
    modules = as.numeric(module.cluster.list[[index]])
    detailed.modules = filtered.bic.module.list[modules]
    combined.regulators = lapply(detailed.modules, function(module){
      module$selected.regulators
    })
    combined.regulators = unique(unlist(combined.regulators))
    combined.targets = lapply(detailed.modules, function(module){
      module$selected.targets
    }) 
    combined.targets = unique(unlist(combined.targets))
    return(list(regulators = combined.regulators, targets = combined.targets))
  })
  save(final.module.list, file = paste0("results/", cancer.type, "_Step6_final.module.list.", it.no, ".rda"))
  
  x = sapply (final.module.list, function(r) r$regulators)
  this_reg_size = length(unique(unlist(x)))
  this_module_size = length(final.module.list)
  y = sapply (final.module.list, function(r) r$targets)
  this_target_size = length(unique(unlist(y)))
  drop = ((prev_reg_size+prev_target_size)-(this_reg_size+this_target_size))/(prev_reg_size+prev_target_size)
  if(it.no != 0){
    cat(paste0('Drop is ', round(drop, digits=3), '\n'))
  }
  cat(paste0('It_',it.no, ': Unique ', this_reg_size, ' regulators and ', this_target_size, ' targets - ', this_module_size,' modules\n'))
  
  
  if(drop >=drop_thr || it.no == 1){
    it.no = it.no + 1
  }else{
    load(paste0("results/", cancer.type, "_Step6_final.module.list.", (it.no-1), ".rda"))
    save(final.module.list, file = paste0("results/", cancer.type, '_final.module.list.rda'))
    cat(paste0('--------- Not sufficient drop between Iterations ', (it.no-1), ' and ', it.no, ' ---------\n'))
    cat(paste0((it.no-1), '.iteration is selected!\n'))
    cat("CANMOD is done!\n")
    cat(paste0("Final module list is saved as \'", cancer.type, "_final.module.list.rda\' under \'results\' folder!\n"))
    break
  }
}

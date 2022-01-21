#write the output to a txt file
for(i in 1:17){
  if(i==1){
    load("~/CanMod_BRCA/Step4_simplified.regulator.target.list_0.1.rda") # load after step4 simplified.regulator.target.list
  }else{
    simplified.regulator.target.list<-final.module.list
  }
  
  module.expression.list = lapply(1:length(simplified.regulator.target.list), function(index){
    module = simplified.regulator.target.list[[index]]
    target.df = expression.df[module$targets,, drop = F]
    regulator.df = expression.df[module$regulators,, drop = F]
    module = list(target.df=target.df,
                  regulator.df=regulator.df)
  })
  
  {
    {
      means = vector()
      for(index in 1:length(module.expression.list)){
        module <- module.expression.list[[index]]
        exp.values = expression.cor[rownames(module$target.df), rownames(module$regulator.df), drop = F]
        exp = (as.matrix(abs(exp.values)))*100
        my.delta.max= summary(as.vector(exp))[6]
        my.delta.min= summary(as.vector(exp))[1]
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
      my.modules = simplified.regulator.target.list
      for(index in 1:length(module.expression.list)){
        module <- module.expression.list[[index]]
        exp.values = expression.cor[rownames(module$target.df), rownames(module$regulator.df)]
        exp = (as.matrix(abs(exp.values)))*100
        my.delta.max= summary(as.vector(exp))[6]
        my.delta.min= summary(as.vector(exp))[1]
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
            selection.tar = t(testCluster@RowxNumber)
            selection.reg = testCluster@NumberxCol
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
        my.delta.max= summary(as.vector(exp))[6]
        my.delta.min= summary(as.vector(exp))[1]
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
      my.modules = simplified.regulator.target.list
      for(index in 1:length(module.expression.list)){
        module <- module.expression.list[[index]]
        exp.values = expression.cor[rownames(module$target.df), rownames(module$regulator.df)]
        exp = (t(as.matrix(abs(exp.values))))*100
        my.delta.max= summary(as.vector(exp))[6]
        my.delta.min= summary(as.vector(exp))[1]
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
            
            selection.reg = t(testCluster@RowxNumber)
            selection.tar = testCluster@NumberxCol
            
            my.modules[[index]]$targets = colnames(exp)[colSums(selection.tar[cluster.tokeep,,drop = F])>0]
            my.modules[[index]]$regulators = rownames(exp)[colSums(selection.reg[cluster.tokeep,,drop = F])>0]
          }
        }else{
        }
      }
    }
    my.modules2 = my.modules; rm(my.modules)
    
    my.modules = mapply(function(x, y) {list(list(selected.regulators = setdiff(union(x$regulators, y$regulators), NA), 
                                                  selected.targets = setdiff(union(x$targets, y$targets), NA)))}, my.modules1, my.modules2)
    
    #my.modules= my.modules2
    toremove = unlist(sapply(my.modules, function(m) (length(m$selected.regulators) == 0 || length(m$selected.targets) == 0)||is.na(m$selected.regulators)[1] || is.na(m$selected.targets)[1]))
    simplified.regulator.target.list = my.modules[!toremove]
    #simplified.regulator.target.list = my.modules
    
    module.expression.list = lapply(1:length(simplified.regulator.target.list), function(index){
      module = simplified.regulator.target.list[[index]]
      target.df = expression.df[module$selected.targets,, drop = F]
      regulator.df = expression.df[module$selected.regulators,, drop = F]
      module = list(target.df=target.df,
                    regulator.df=regulator.df)
    })
    
    filtered.bic.module.list = simplified.regulator.target.list
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
      #if (row.index %% 100 == 0) #print(row.index)
      #row.index=3
      x = rownames(target.module.df)[row.index]
      #cat("target: ",x,"\n")
      modules = as.numeric(unique(c(filtered.bic.target.df$module[which(filtered.bic.target.df$target == x)],
                                    filtered.bic.regulator.df$module [which(filtered.bic.regulator.df$regulator == x)])))
      #cat("module: ",modules,"\n")
      target.module.df[x,modules]= 1
      #cat("row.index: ",row.index,"\n")
    }
  }
  
  {
    sm = reshape2::melt(1-as.matrix(dist(t(target.module.df), method = "binary")))
    sm = sm[sm$value >=0.8,]
    graph = graph_from_data_frame(sm[,c(1, 2)], directed = F) %>%
      set_edge_attr("weight", value = as.numeric(sm$value))
    module.cluster.list = as.list(cluster_walktrap(graph, weights =  E(graph)$weight)); 
    # rm(sm, d, df, graph)
    # module.cluster.list = as.list(cluster_walktrap(graph))
  }
  # for each cluster modules, recruit back targets 
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
  
  save(filtered.bic.module.list, file = paste0(cancer.type, "_Step5.gc45.rc10.it",i,".rda"))
  save(final.module.list, file = paste0(cancer.type, "_Step6.gc45.rc10.it",i,".rda"))
  
  
  filtered.bic.module.list.unique.gene<-lapply(1:length(filtered.bic.module.list),function(i){
    return (filtered.bic.module.list[[i]][["selected.targets"]])
  })
  
  filtered.bic.module.list.unique.gene.len<-lapply(1:length(filtered.bic.module.list),function(i){
    return (length(filtered.bic.module.list[[i]][["selected.targets"]]))
  })
  filtered.bic.module.list.unique.gene.len.df<-do.call(rbind,filtered.bic.module.list.unique.gene.len)
  summary(filtered.bic.module.list.unique.gene.len.df)
  plot(density(filtered.bic.module.list.unique.gene.len.df))
  table(filtered.bic.module.list.unique.gene.len.df)
  
  filtered.bic.module.list.unique.regulators<-lapply(1:length(filtered.bic.module.list),function(i){
    return (filtered.bic.module.list[[i]][["selected.regulators"]])
  })
  
  
  final.module.list.unique.gene<-lapply(1:length(final.module.list),function(i){
    return (final.module.list[[i]][["targets"]])
  })
  
  
  final.module.list.unique.regulator<-lapply(1:length(final.module.list),function(i){
    return (final.module.list[[i]][["regulators"]])
  })
  
  
  cat(length(filtered.bic.module.list))
  cat(" ", length(unique(unlist(filtered.bic.module.list.unique.gene))))
  cat(" ",length(unique(unlist(filtered.bic.module.list.unique.regulators))),"\n")
  
  cat(length(final.module.list))
  cat(" ",length(unique(unlist(final.module.list.unique.gene))))
  cat(" ",length(unique(unlist(final.module.list.unique.regulator))),"\n")
  
}




a<-read.delim2("3.txt",sep = "",header = F) # read a plain txt file of 34 lines.
colnames(a)<-c("#Module","#Gene","#Regulator")
i=1
module.begin<-393 #after step4 number of modules
gene.begin<-1326 #after step4 number of genes
regulator.begin<-133 #after step4 number of regulator
while(i<35){
  if(i==1){
    a$step5.module.drop[i]<-(module.begin-a$`#Module`[i])/module.begin
    a$step5.gene.drop[i]<-(gene.begin-a$`#Gene`)/gene.begin
    a$step5.regulator.drop[i]<-(regulator.begin-a$`#Regulator`)/regulator.begin
  }else{
    if(i%%2==1){
      a$step5.module.drop[i]<-(a$`#Module`[i-2]-a$`#Module`[i])/a$`#Module`[i-2]
      a$step5.gene.drop[i]<-(a$`#Gene`[i-2]-a$`#Gene`[i])/a$`#Gene`[i-2]
      a$step5.regulator.drop[i]<-(a$`#Regulator`[i-2]-a$`#Regulator`[i])/a$`#Regulator`[i-2]
    }
    if(i%%2==0){
      a$step5.module.drop[i]<-NA
      a$step5.gene.drop[i]<-NA
      a$step5.regulator.drop[i]<-NA
    }
  }
  i=i+1
  cat(i,"\n")
}

write.csv(a,file="BRCA45RC20.csv") # write it to a csv file
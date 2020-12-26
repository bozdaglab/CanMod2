#len.final.target.gene<-list()
cor.final.target.gene<-vector()
sig<-0
pairwise.cor.target.gene2<-lapply(1:20, function(i){
  
  final.target.gene<-dtc.final.module.list[[i]][["targets"]]
  #len.final.target.gene<-rlist::list.append(len.final.target.gene,length(final.target.gene))
  cor.final.target.gene[i]<-mean(abs(expression.cor[final.target.gene,final.target.gene]))
  f<-0
  cor.random.target.gene<-vector()
  for (j in 1:1000){
    random.target.gene<-sample(x=rownames(variable.selection.data$mRNA),size=length(final.target.gene),replace= FALSE, prob= NULL)
    cor.random.target.gene[j]<-mean(abs(expression.cor[random.target.gene,random.target.gene]))
    if(cor.final.target.gene[i]<cor.random.target.gene[j]){
      f<-f+1
    }
  }
  p.value<-f/1000
  if(p.value<.05){
    sig<-sig+1
  }
  
  return(list(cor.final.target.gene,sig,cor.random.target.gene))
  
})

sum(sapply(pairwise.cor.target.gene2, function(module){module[[2]] == 1}))
save(pairwise.cor.target.gene, file = "codes/test/test_final/ev/pairwise.cor.target.gene.rda")

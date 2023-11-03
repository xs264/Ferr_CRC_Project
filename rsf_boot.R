rsf_boot.os5 <- function(seed) {
  set.seed(seed)
  boot.indices <- sample(1:nrow(ml), replace = TRUE)
  ml.boot <- ml[boot.indices, ]
  set.seed(seed)
  best.mod<-rfsrc( Surv(time = os5.months, event = os5.status) ~., 
                   data =ml.boot, ntree = ntree, nsplit = nsplit,  nodesize=nodesize, forest = T)
  
  partial_nocov1 <- plot.variable(best.mod, xvar.names = gene_list, partial=TRUE, 
                                  show.plot=F,surv.type = "surv",time=36,smooth.lines = F)
  
  partial_nocov2 <- plot.variable(best.mod, xvar.names = gene_list, partial=TRUE, 
                                  show.plot=F,surv.type = "surv",time=60,smooth.lines = F)
  
  partial.data<-data.frame()
  for(i in 1:length(gene_list)){
    name <- gene_list[i] 
    partial.dt <- partial_nocov1$pData[[i]]
    group <- rep("Time = 3-year",length(partial.dt$x.uniq))
    gene.name <- rep(name,length(partial.dt$x.uniq))
    Gene.val <- partial.dt$x.uniq
    yhat <- partial.dt$yhat
    result<-cbind(gene.name=gene.name, Gene.exp=Gene.val,
                  survprob=yhat,
                  Time_point=group)
    partial.data<-rbind(partial.data,result)
  }
  for(i in 1:length(gene_list)){
    name <- gene_list[i] 
    partial.dt <- partial_nocov2$pData[[i]]
    group <- rep("Time = 5-year",length(partial.dt$x.uniq))
    gene.name <- rep(name,length(partial.dt$x.uniq))
    Gene.val <- partial.dt$x.uniq
    yhat <- partial.dt$yhat
    result<-cbind(gene.name=gene.name, Gene.exp=Gene.val,
                  survprob=yhat,
                  Time_point=group)
    partial.data<-rbind(partial.data,result)
  }
  
  output <- list("gene.name" = partial.data$gene.name, 
                 "Gene.exp" = partial.data$Gene.exp,
                 "survprob"=partial.data$survprob,
                 "Time_point"=partial.data$Time_point)
  return(output)
  
}
rsf_boot.os <- function(seed) {
  set.seed(seed)
  boot.indices <- sample(1:nrow(ml), replace = TRUE)
  ml.boot <- ml[boot.indices, ]
  set.seed(seed)
  best.mod<-rfsrc( Surv(time = os.months, event = os.status) ~., 
                   data =ml.boot, ntree = ntree, nsplit = nsplit,  nodesize=nodesize, forest = T)
  
  partial_nocov1 <- plot.variable(best.mod, xvar.names = gene_list, partial=TRUE, 
                                  show.plot=F,surv.type = "surv",time=time1,smooth.lines = smooth.lines,npts =npts)
  
  partial_nocov2 <- plot.variable(best.mod, xvar.names = gene_list, partial=TRUE, 
                                  show.plot=F,surv.type = "surv",time=time2,smooth.lines = smooth.lines,npts =npts)
  
  partial.data<-data.frame()
  for(i in 1:length(gene_list)){
    name <- gene_list[i] 
    partial.dt <- partial_nocov1$pData[[i]]
    group <- rep("Time = 3-year",length(partial.dt$x.uniq))
    gene.name <- rep(name,length(partial.dt$x.uniq))
    Gene.val <- partial.dt$x.uniq
    yhat <- partial.dt$yhat
    result<-cbind(gene.name=gene.name, Gene.exp=Gene.val,
                  survprob=yhat,
                  Time_point=group)
    partial.data<-rbind(partial.data,result)
  }
  for(i in 1:length(gene_list)){
    name <- gene_list[i] 
    partial.dt <- partial_nocov2$pData[[i]]
    group <- rep("Time = 5-year",length(partial.dt$x.uniq))
    gene.name <- rep(name,length(partial.dt$x.uniq))
    Gene.val <- partial.dt$x.uniq
    yhat <- partial.dt$yhat
    result<-cbind(gene.name=gene.name, Gene.exp=Gene.val,
                  survprob=yhat,
                  Time_point=group)
    partial.data<-rbind(partial.data,result)
  }
  
  output <- list("gene.name" = partial.data$gene.name, 
                 "Gene.exp" = partial.data$Gene.exp,
                 "survprob"=partial.data$survprob,
                 "Time_point"=partial.data$Time_point)
  return(output)
  
}
rsf_boot.rfs <- function(seed) {
  set.seed(seed)
  boot.indices <- sample(1:nrow(ml), replace = TRUE)
  ml.boot <- ml[boot.indices, ]
  set.seed(seed)
  best.mod<-rfsrc( Surv(time = rfs.months, event = rfs.status) ~., 
                   data =ml.boot, ntree = ntree, nsplit = nsplit,  nodesize=nodesize, forest = T)
  
  partial_nocov1 <- plot.variable(best.mod, xvar.names = gene_list, partial=TRUE, 
                                  show.plot=F,surv.type = "surv",time=36,smooth.lines = F)
  
  partial_nocov2 <- plot.variable(best.mod, xvar.names = gene_list, partial=TRUE, 
                                  show.plot=F,surv.type = "surv",time=60,smooth.lines = F)
  
  partial.data<-data.frame()
  for(i in 1:length(gene_list)){
    name <- gene_list[i] 
    partial.dt <- partial_nocov1$pData[[i]]
    group <- rep("Time = 3-year",length(partial.dt$x.uniq))
    gene.name <- rep(name,length(partial.dt$x.uniq))
    Gene.val <- partial.dt$x.uniq
    yhat <- partial.dt$yhat
    result<-cbind(gene.name=gene.name, Gene.exp=Gene.val,
                  survprob=yhat,
                  Time_point=group)
    partial.data<-rbind(partial.data,result)
  }
  for(i in 1:length(gene_list)){
    name <- gene_list[i] 
    partial.dt <- partial_nocov2$pData[[i]]
    group <- rep("Time = 5-year",length(partial.dt$x.uniq))
    gene.name <- rep(name,length(partial.dt$x.uniq))
    Gene.val <- partial.dt$x.uniq
    yhat <- partial.dt$yhat
    result<-cbind(gene.name=gene.name, Gene.exp=Gene.val,
                  survprob=yhat,
                  Time_point=group)
    partial.data<-rbind(partial.data,result)
  }
  
  output <- list("gene.name" = partial.data$gene.name, 
                 "Gene.exp" = partial.data$Gene.exp,
                 "survprob"=partial.data$survprob,
                 "Time_point"=partial.data$Time_point)
  return(output)
  
}
rsf_boot.pfs5 <- function(seed) {
  set.seed(seed)
  boot.indices <- sample(1:nrow(ml), replace = TRUE)
  ml.boot <- ml[boot.indices, ]
  set.seed(seed)
  best.mod<-rfsrc( Surv(time = pfs5.months, event = pfs5.status) ~., data =ml.boot, ntree = ntree, nsplit = nsplit,  nodesize=nodesize, forest = T)
  
  partial_nocov1 <- plot.variable(best.mod, xvar.names = gene_list, partial=TRUE, 
                                  show.plot=F,surv.type = "surv",time=36,smooth.lines = F)
  
  partial_nocov2 <- plot.variable(best.mod, xvar.names = gene_list, partial=TRUE, 
                                  show.plot=F,surv.type = "surv",time=60,smooth.lines = F)
  
  partial.data<-data.frame()
  for(i in 1:length(gene_list)){
    name <- gene_list[i] 
    partial.dt <- partial_nocov1$pData[[i]]
    group <- rep("Time = 3-year",length(partial.dt$x.uniq))
    gene.name <- rep(name,length(partial.dt$x.uniq))
    Gene.val <- partial.dt$x.uniq
    yhat <- partial.dt$yhat
    result<-cbind(gene.name=gene.name, Gene.exp=Gene.val,
                  survprob=yhat,
                  Time_point=group)
    partial.data<-rbind(partial.data,result)
  }
  for(i in 1:length(gene_list)){
    name <- gene_list[i] 
    partial.dt <- partial_nocov2$pData[[i]]
    group <- rep("Time = 5-year",length(partial.dt$x.uniq))
    gene.name <- rep(name,length(partial.dt$x.uniq))
    Gene.val <- partial.dt$x.uniq
    yhat <- partial.dt$yhat
    result<-cbind(gene.name=gene.name, Gene.exp=Gene.val,
                  survprob=yhat,
                  Time_point=group)
    partial.data<-rbind(partial.data,result)
  }
  
  output <- list("gene.name" = partial.data$gene.name, 
                 "Gene.exp" = partial.data$Gene.exp,
                 "survprob"=partial.data$survprob,
                 "Time_point"=partial.data$Time_point)
  return(output)
  
}
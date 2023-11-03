# Survival analysis one-by-one
FML <-function(surv,varb,covariate) {
  as.formula(paste0(surv,paste( append(varb,covariate),collapse = '+')))
}
#example: call FML("1") and get "Surv() ~ metabolite + SexGroup + SideGroup + stageGroup + age"
#Get all the results function
getall<-function(i,data,num_pred){
  sum<-summary(coxph(i,data = data))
  HR <- round(sum$coefficients[,2],2)
  P_Value <- round(sum$coefficients[,5],4)
  LCI <- format(round(sum$conf.int[,3],2),2)
  UCI <- format(round(sum$conf.int[,4],2),2)
  CI95 <- paste0(LCI,'-',UCI)
  SE <- sum$coefficients[,3]
  # show statistical significance
  sig<- rep(NA,num_pred)
  for (j in 1:num_pred){
    if (P_Value[j]<0.05) {sig[j]="*"}
    else if (P_Value[j] <0.1) {sig[j]="."}
    else {sig[j]=""}
  }
  result<-cbind(HR=HR,CI95=CI95,SE=SE,LCI=LCI,UCI=UCI,P_Value=P_Value,sig=sig)
} 

getall.exact<-function(i,data,num_pred){
  sum<-summary(coxph(i,data = data))
  HR <- sum$coefficients[,2]
  P_Value <- sum$coefficients[,5]
  LCI <- sum$conf.int[,3]
  UCI <- sum$conf.int[,4]
  CI95 <- paste0(LCI,'-',UCI)
  SE <- sum$coefficients[,3]
  # show statistical significance
  sig<- rep(NA,num_pred)
  for (j in 1:num_pred){
    if (P_Value[j]<0.05) {sig[j]="*"}
    else if (P_Value[j] <0.1) {sig[j]="."}
    else {sig[j]=""}
  }
  result<-cbind(HR=HR,CI95=CI95,SE=SE,LCI=LCI,UCI=UCI,P_Value=P_Value,sig=sig)
} 

coxphtable<-function(met.list,Surv,covar,data){
  out_multi=data.frame()
  num_pred=1+length(covar)
  num_met <- length(met.list)
  for (i in met.list){
    out<-as.data.frame(getall.exact(FML(Surv,i,covar),data,num_pred))
    out_multi<-rbind(out_multi,out)
  }
  #only extract rows that contain outputs of biomarker :
  out_multi_mtlog=data.frame()
  for (i in 1:length(met.list)){
    out_multi_mtlog<-rbind(out_multi[num_pred*(num_met-i)+1,],out_multi_mtlog)
  }
  # output<-cbind(met_name,out_multi_mtlog)
  P_Value <- out_multi_mtlog$P_Value
  out_multi_mtlog$FDR.P <- round(p.adjust(as.numeric(P_Value), method = "fdr", n = num_met),3)
  output<-out_multi_mtlog
  return(output)
}

FML.int <-function(surv,varb,covar) {
  as.formula(paste0(surv,paste(varb,"+","sex*",varb,"+",paste(covar,collapse = '+'))))
}

coxphtable.int<-function(name.os.sex,Surv,covar,data){
  out_multi=data.frame()
  num_pred=2+length(covar)
  num_met <- length(name.os.sex)
  for (i in name.os.sex){
    out<-as.data.frame(getall.exact(FML.int(Surv,i,covar),data,num_pred))
    out_multi<-rbind(out_multi,out)
  }
  #only extract rows that contain outputs of metabolites :
  out_multi_mtlog=data.frame()
  for (i in 1:length(name.os.sex)){
    out_multi_mtlog<-rbind(out_multi[num_pred*(num_met-i)+num_pred,],out_multi_mtlog)
  }
  # out_multi_mtlog
  # output<-cbind(met_name,out_multi_mtlog)
  # P_Value <- out_multi_mtlog$P_Value
  output<-out_multi_mtlog[c(6,7)]
  return(output)
}

#' Computes AUC of AROC and ROC along with Specificity and Sensitivity of the Youden Index
#'
#' @param data Healthy control as 0 and Disease control as 1. csv file with header as false. First column defining Sample information. Second column = Group, Third Column = Covariate ( Must be included even with random data to process the code). Rest of the columns are to be used as metabolites
#' @param cov T,F argument specifying whether covariate adjsuted ROC will be computed
#'
#' @return Either a list (adjusted + not adjusted) or a data matrix
#' @export
#'
#' @examples
#'
#' D<-as.data.frame(cbind(paste0("sample_",1:26),c(rep(1,13),rep(0,13)),c(rep("2",7),rep("3",12),rep("4",7)),runif(26),runif(26)))
#' Datas<-rbind(c("Sample","Group","Covariate","some metabolite1","some metabolite2"),D)
#' roc_all(Datas)

roc_all<-function(data,cov=T){
  if (cov ==T){
    colnames(data)<-data[1,];data<-data[-1,]
    rownames(data)<-data[,1];data<-data[,-1]
    data_Nadj<-data[,-2]
    for(n in 1:ncol(data_Nadj)){data_Nadj[,n]<-as.numeric(data_Nadj[,n])}

    Result_unadj<-matrix(data = NA,nrow = (ncol(data_Nadj)-1),ncol = 3)
    rownames(Result_unadj)<-colnames(data_Nadj)[2:ncol(data_Nadj)]
    colnames(Result_unadj)<-c("AUC","Specificity","Sensitivity")
    for (i in 2:ncol(data_Nadj)){
      roc_res <- pROC::roc(response = data_Nadj$Group, predictor = data_Nadj[,i], percent=TRUE,quiet = T)
      Result_unadj[(i-1),1]<-roc_res[["auc"]]
      for (n in 2:3){
        Result_unadj[(i-1),n]<-pROC::coords(roc_res, "best",  transpose = FALSE,best.method="youden")[1,n]}
    }
    data_adj<-data
    for(n in 1:ncol(data_adj)){data_adj[,n]<-as.numeric(data_adj[,n])}

    ##### adjusted #####
    Result_adj <-matrix(data = NA,nrow = (ncol(data_adj)-2),ncol = 3)
    rownames(Result_adj)<-colnames(data_adj)[3:ncol(data_adj)];colnames(Result_adj)<-c("AUC","Specificity","Sensitivity")
    LETTERS702 <- c(sapply(LETTERS, function(x) paste0(x, LETTERS)))
    LETTERS37232 <- c(LETTERS,LETTERS702, sapply(LETTERS, function(x) paste0(x,
                                                                             LETTERS702)))
    LETTERS37232<-LETTERS37232[-365]
    colnames(data_adj)[3:ncol(data_adj)]<-LETTERS37232[1:(ncol(data_adj)-2)]
    for (i in 3:ncol(data_adj)){
      AROC <-ROCnReg::AROC.sp(paste0(colnames(data_adj)[i],"~Covariate") ,group = "Group",tag.h = "0",data = data_adj)
      Tem<-aroc_compute(AROC)
      Result_adj_YI<-matrix(ncol = 3, nrow = 500)
      colnames(Result_adj_YI)<-c("YI","Specificity","Sensitivity")
      Result_adj_YI[,1]<-(Tem$specificity+Tem$sensitivity-1);Result_adj_YI[,2]<-Tem$specificity;Result_adj_YI[,3]<-Tem$sensitivity
      Result_adj_YI<-as.data.frame(Result_adj_YI)
      Result_adj[(i-2),1]<-AROC[["AUC"]][["est"]]
      for (n in 2:3){
        Result_adj[(i-2),n]<-Result_adj_YI[order(Result_adj_YI$YI,decreasing = T),][1,][,n]}
    }
    Final<-list()
    Final$adj<-Result_adj*100
    Final$unadj<-Result_unadj
    Final
  }
  else if (cov ==F){
    colnames(data)<-data[1,];data<-data[-1,]
    rownames(data)<-data[,1];data<-data[,-1]
    data[,1]<-as.factor(data[,1])
    data_Nadj<-data[,-2]
    for(n in 2:ncol(data_Nadj)){data_Nadj[,n]<-as.numeric(data_Nadj[,n])}

    Result_unadj<-matrix(data = NA,nrow = (ncol(data_Nadj)-1),ncol = 3)
    rownames(Result_unadj)<-colnames(data_Nadj)[2:ncol(data_Nadj)]
    colnames(Result_unadj)<-c("AUC","Specificity","Sensitivity")
    for (i in 2:ncol(data_Nadj)){
      roc_res <- pROC::roc(response = data_Nadj$Group, predictor = data_Nadj[,i], percent=TRUE,quiet = T)
      Result_unadj[(i-1),1]<-roc_res[["auc"]]
      for (n in 2:3){
        Result_unadj[(i-1),n]<-pROC::coords(roc_res, "best",  transpose = FALSE,best.method="youden")[1,n]}
    }
    Result_unadj
  }
  else (print("need to specifcy the cov argument"))
}

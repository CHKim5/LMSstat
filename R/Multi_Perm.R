#' Computes PERMANOVA with vegan using multiple group information
#'
#' @param Data csv file with Header as False First column with Sample Second column with Multilevel(Mixomics) so that it can be compatible with other multivariate statistics Third column with Group information. Rest of the columns are the metabolites to be tested.
#' @param Classification csv file with Header as False ;First column with Sample; every other column can be used to indicate multiple classifiers
#' @param method Dissimilarity index c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq",chord")
#'
#' @return vegan::adonis2(Dist_Met~Group,method = "bray",by = NULL,data = x_y_coord_G)
#' @export
#'
#' @examples D<-as.data.frame(cbind(LETTERS,rep(1,26),c(rep("A",7),rep("B",6),rep("C",13)),runif(26),runif(26)))
#' Datas<-rbind(c("Sample","Multilevel","Group","some metabolite1","some metabolite2"),D)
#' Datas_G<-rbind(c("Sample","some class0","some class1","some class2"),D)
#' Multi_Perm(Datas,Datas_G)
#' Note that the code recognizes each class as a factor.
Multi_Perm<-function(Data,Classification,method="bray"){
  colnames(Data)<-Data[1,];Data<-Data[-1,];Data<-Data[,-2]
  rownames(Data)<-Data[,1];Data<-Data[,-1]
  colnames(Classification)<-Classification[1,];Classification<-Classification[-1,]
  rownames(Classification)<-Classification[,1];Classification<-Classification[,-1]
  Data<-Data[,-1]
  for (n in 1: ncol(Classification)){Classification[,n] <-as.factor(Classification[,n])}
  for (n in 1:ncol(Data)){Data[,n]<-as.numeric(Data[,n])}
  Dist_Met<-as.matrix(vegan::vegdist(Data,method = method))
  NMDS=vegan::metaMDS(Dist_Met,k = 2,trymax = 2000)#Bray curtis
  x_y_coord<-as.data.frame(vegan::scores(NMDS,display = "sites"))
  result<-list()
for (n in 1:ncol(Classification)){
  x_y_coord_G <- cbind(x_y_coord, Classification[,n])
  colnames(x_y_coord_G)[3]<-"Group"
  M_p_value<-vegan::adonis2(Dist_Met~Group,method = method,by = NULL,data = x_y_coord_G)
  eval(parse(text = paste0("result$",colnames(Classification)[n],"<-M_p_value")))
}
  result
}

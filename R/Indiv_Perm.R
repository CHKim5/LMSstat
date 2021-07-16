#' Computes PERMANOVA with vegan using specific group information
#'
#' @param Data csv file with Header as False First column with Sample Second column with Multilevel(Mixomics) so that it can be compatible with other multivariate statistics Third column with Group information. Rest of the columns are the metabolites to be tested.
#' @param method Dissimilarity index c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq",chord")
#'
#' @return vegan::adonis2(Dist_Met~Group,method = "bray",by = NULL,data = x_y_coord_G)
#' @export
#'
#' @examples Df<-data(Data)
#' Indiv_Perm(Df)
Indiv_Perm<-function(Data,method="bray"){
  colnames(Data)<-Data[1,];Data<-Data[-1,];Data<-Data[,-2]
  rownames(Data)<-Data[,1];Data<-Data[,-1]
  Data[,1]<-as.factor(Data[,1]) # Group as factors
  for (n in 2:ncol(Data)){Data[,n]<-as.numeric(Data[,n])}
  Dist_Met<-as.matrix(vegan::vegdist(Data[,2:ncol(Data)],method = method))
  NMDS=vegan::metaMDS(Dist_Met,k = 2,trymax = 2000)#Bray curtis
  x_y_coord<-as.data.frame(vegan::scores(NMDS,display = "sites"))
  x_y_coord_G <- cbind(x_y_coord, Data$Group)
  colnames(x_y_coord_G)[3]<-"Group"
  M_p_value<-vegan::adonis2(Dist_Met~Group,method = method,by = NULL,data = x_y_coord_G)
  M_p_value
}

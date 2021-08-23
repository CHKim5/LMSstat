#' Computes Scheirer-Ray-Hare Test
#'
#' @param data csv file with Header as False First column with Sample Second column with Multilevel(Mixomics) so that it can be compatible with other multivariate statistics Third column with Group information. Rest of the columns are the metabolites to be tested.
#' @param Adjust_p_value Set True if FDR adjustments are to be made. If not set False
#' @param Adjust_method adjustment methods frequently used. "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"
#'
#' @return csv files including significant variables (Multilevel, Group, interaction) and a Venn diagram
#' @export
#'
#' @examples data(Data)
#' Result<-SRH(Data)
SRH<-function(data,Adjust_p_value = T,Adjust_method = "BH"){
  LETTERS702 <- c(sapply(LETTERS, function(x) paste0(x, LETTERS)))
  LETTERS37232 <- c(LETTERS, LETTERS702, sapply(LETTERS, function(x) paste0(x,
                                                                            LETTERS702)))
  LETTERS210729 <- c(LETTERS, LETTERS702, LETTERS37232, sapply(LETTERS,
                                                               function(x) paste0(x, LETTERS37232)))
  LETTERS210729 <- LETTERS210729[-365];LETTERS210729 <- LETTERS210729[-10205];LETTERS210729 <- LETTERS210729[-267101]
  colnames(data)<-data[1,]; data<-data[-1,]
  df<-as.data.frame(matrix(data = NA,ncol = 3,nrow=(ncol(data)-3)))
  df_name<-as.data.frame(matrix(data = NA,ncol = 2,nrow=(ncol(data)-3)))
  colnames(df)<-c("Multilevel","Group","Interaction");colnames(df_name)<-c("Ori","LETTERS")
  df_name$Ori<-colnames(data)[4:ncol(data)]
  df_name$LETTERS<-LETTERS210729[1:nrow(df_name)]
  colnames(data)[4:ncol(data)]<-df_name$LETTERS
  print("Started SRH")
  for(a in 1:nrow(df_name)){
    eval(parse(text = paste0("df[",a,",]<-rcompanion::scheirerRayHare(",colnames(data)[a+3],"~Multilevel*Group,data = data,verbose = F)$p.value[1:3]")))
  }
  print("Finished SRH")
  rownames(df)<-df_name$Ori

  {if (Adjust_p_value == T) {
    print("###########################################")
    print(paste0("adjusted according to the ",
                 Adjust_method, " method"))
    print("###########################################")
    df <- apply(df, 2, function(x) {
      p.adjust(x, method = Adjust_method)
    })
    df<-as.data.frame(df)

  }
    else {
      print("###########################################")
      print("p_value not adjusted")
      print("###########################################")
    }}



  ML_sig<-df %>% dplyr::filter(Multilevel<0.05)
  G_sig<-df %>% dplyr::filter(Group<0.05)
  Int_sig<-df %>% dplyr::filter(Interaction<0.05)
  Result<-list()
  Result[["Multilevel"]]<-ML_sig
  Result[["Group"]]<-G_sig
  Result[["Interaction"]]<-Int_sig
  Result[["Result"]]<-df
  print("Started to print Venn_diagram")
  VennDiagram::venn.diagram(list(rownames(ML_sig),rownames(G_sig),rownames(Int_sig)),
                            category.names = c("Multilevel","Group","Interaction"),
                            filename = "Significant_variable_venndiagram.jpg",
                            resolution = 1200,
                            col=c("#440154ff", '#21908dff', '#fde725ff'),
                            fill = c(ggplot2::alpha("#440154ff",0.3), ggplot2::alpha('#21908dff',0.3), ggplot2::alpha('#fde725ff',0.3)),
                            cex = 4,
                            lwd = 1,
                            height = 16 ,
                            units = "in",
                            width = 16,
                            cat.cex = 4,
                            cat.pos = c(-27, 27, 135),
                            cat.dist = c(0.055, 0.055, 0.085) )
  print("Finished printing Venn_diagram")
  Result
  write.csv(G_sig,"Significant for variable_Group.csv")
  write.csv(ML_sig,"Significant for variable_Multilevel.csv")
  write.csv(Int_sig,"Significant for variable_Interaction.csv")
  write.csv(df,"Result.csv")
}

#' Automatically save PCA plot for the selected dimension metrics
#'
#' @param Data uploaded data
#' @param color colors used for ggplots.color=c("#FF3300","#FF6600","#FFCC00","#99CC00","#0066CC","#660099")
#' @param legend_position legend position "none","left","right","bottom","top"
#' @param fig_width figure size
#' @param fig_height figure size
#' @param components PCA components c(1,3)
#' @param names whether the names are to be indicated T,F
#' @param labsize label size
#' @param dotsize dot size
#'
#' @return PCA plot
#' @export
#'
#' @examples data(Data)
#' PCA(Data,components = c(1,2),legend_position = "none")
PCA<-function(Data,
               color = c("#FF3300","#660099","#FFCC00","#99CC00","#0066CC","#FF6600"),
               legend_position = "none",
               fig_width = 24,
               fig_height = 20,
               components = c(1,2),
              names = F,
              dotsize = 3,
              labsize = 3
               ){
  ifelse(!dir.exists(file.path(getwd(), "PCA")), dir.create(file.path(getwd(), "PCA")), FALSE)
  colnames(Data) <- Data[1, ]
  Data <- Data[-1, -2]
  Data<-Data %>% dplyr::arrange(Group)
  rownames(Data)<-Data[,1];Data<-Data[,-1]
  Data[,2:ncol(Data)]<-sapply(Data[,2:ncol(Data)],function(x) as.numeric(x))
  PCA_Data_input<-Data#subset samples of interest, columns: Group(Male, Female etc)+metabolite, rownames : sample names
  PCA_Data_input[,1]<-as.factor(PCA_Data_input[,1]) # Group as factors
  Bef_merge<-prcomp(PCA_Data_input[,2:ncol(PCA_Data_input)],scale = T)
  PCA_res<-cbind(as.data.frame(Bef_merge$x),Data$Group)
  colnames(PCA_res)[ncol(PCA_res)]<-"Group"
  explained_variance <- Bef_merge$sdev^2/sum(Bef_merge$sdev^2)
    sem_res<-eval(parse(text=paste0("ggplot2::ggplot(data = PCA_res, ggplot2::aes(x = PC",components[1],", y =PC",components[2],",
                                                       color = Group,label = rownames(PCA_res)))")))+
    ggplot2::geom_point(size = dotsize)+
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(values = color)+
    ggplot2::labs(x = paste0("PC",components[1],": ",round(explained_variance[components[1]]*100,1),"%"),
                  y = paste0("PC",components[2],": ",round(explained_variance[components[2]]*100,1),"%"))+
    ggplot2::theme(plot.title=ggplot2::element_text(size=10,
                                                    face="bold",
                                                    color="Black",
                                                    hjust=0.5,
                                                    lineheight=1.2),  # title
                   plot.subtitle=ggplot2::element_text(size=15,
                                                       hjust=0.5),  # subtitle
                   plot.caption=ggplot2::element_text(size=15),  # caption
                   axis.title.x=ggplot2::element_text(vjust=0,
                                                      size=10),  # X axis title
                   axis.title.y=ggplot2::element_text(size=10),  # Y axis title
                   legend.position = legend_position)+
    eval(parse(text=paste0("ggplot2::stat_ellipse(ggplot2::aes(x = PC",components[1],", y =PC",components[2],",
                                                       color = Group), type = 'norm')")))+
    ggplot2::ggtitle(paste("Component", components[1]," & Component", components[2]," PCA result"))
  if (names == T){
    sem_res+ggrepel::geom_text_repel(size=labsize)}
  else {sem_res}
  ggplot2::ggsave(paste(components[1],components[2],"label",names,"PCA.png",sep = "_"),
                  path=paste0(getwd(),"/PCA"),
                  height = fig_height,
                  width = fig_width,units = "cm")
  Result<-list()
  Result$explained_variance<-explained_variance
  Result$coordinates<-PCA_res
  Result
  }

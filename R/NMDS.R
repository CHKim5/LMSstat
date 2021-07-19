#' Automatically save NMDS plots for the selected dimension metrics
#'
#' @param Data uploaded data
#' @param methods c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn",  "raup", "binomial", "mahalanobis", "chisq","chord")
#' @param color colors used for ggplots.color=c("#FF3300","#FF6600","#FFCC00","#99CC00","#0066CC","#660099")
#' @param legend_position legend position "none","left","right","bottom","top"
#' @param fig_width figure size
#' @param fig_height figure size
#' @param names whether the names are to be indicated T,F
#' @param labsize label size
#' @param dotsize dot size
#'
#' @return NMDS plot
#' @export
#'
#' @examples data(Data)
#' NMDS(Data,methods = c("euclidian","manhattan","bray"))

NMDS<-function(Data,
               methods=c("euclidean","bray"),
               color = c("#FF3300","#660099","#FFCC00","#99CC00","#0066CC","#FF6600"),
               legend_position = "bottom",
               fig_width = 24,
               fig_height = 20,
               dotsize = 3,
               names = F,
               labsize = 3
               ){
  ifelse(!dir.exists(file.path(getwd(), "NMDS")), dir.create(file.path(getwd(), "NMDS")), FALSE)
  colnames(Data) <- Data[1, ]
  Data <- Data[-1, -2]
  Data<-Data %>% dplyr::arrange(Group)
  rownames(Data)<-Data[,1];Data<-Data[,-1]
  Data[,2:ncol(Data)]<-sapply(Data[,2:ncol(Data)],function(x) as.numeric(x))
  NMDS_Data_input<-Data#subset samples of interest, columns: Group(Male, Female etc)+metabolite, rownames : sample names
  NMDS_Data_input[,1]<-as.factor(NMDS_Data_input[,1]) # Group as factors
  for(method in methods){
  Dist_Met<-as.matrix(vegan::vegdist(NMDS_Data_input[,2:ncol(NMDS_Data_input)],method = method))
  NMDS=vegan::metaMDS(Dist_Met,k = 2,trymax = 2000)#Bray curtis
  x_y_coord<-as.data.frame(vegan::scores(NMDS,display = "sites"))
  x_y_coord_G <- cbind(x_y_coord, NMDS_Data_input$Group)
  colnames(x_y_coord_G)[3]<-"Group"
  M_p_value<-vegan::adonis2(Dist_Met~Group,
                     method = method,
                     by = NULL,
                     data = x_y_coord_G)[["Pr(>F)"]][1] #Model P_Value
  sem_res<-ggplot2::ggplot(data = x_y_coord_G, ggplot2::aes(x = NMDS1, y = NMDS2,
                                                   color = Group,label=rownames(x_y_coord_G))) +
    ggplot2::geom_point(size = dotsize)+
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(values = color)+
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
    ggplot2::stat_ellipse(ggplot2::aes(x=NMDS1,
                                       y=NMDS2,color=Group),
                          type = "norm")+
    ggplot2::ggtitle(paste(method,"index/distance, ","Stress Level : ",NMDS$stress,"Model_p_value",M_p_value))
  if (names == T){
    sem_res+ggrepel::geom_text_repel(size=labsize)}
  else {sem_res}
  ggplot2::ggsave(paste(method,"NMDS.png",sep = "_"),
                  path=paste0(getwd(),"/NMDS"),
                  height = fig_height,
                  width = fig_width,units = "cm")
  }}

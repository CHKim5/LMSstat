#' Volcano plot for every group comparison
#'
#' @param data data inheriting from Allstats
#' @param asterisk statistics inheriting from Allstats  # "t-test", "u-test", "Scheffe", "Dunn"
#' @param reverse reverse the direction of fold change
#' @param fig_width figure size
#' @param FC_log Fold change log transformation value
#' @param pval_log p_value log transformation value
#' @param fig_height figure size
#' @param dotsize dotsize
#' @param x_limit x axis limt c( -3,3)
#' @param y_limit y axis limit c(0, 100)
#' @param pval_intercept intercept for identification 0.05
#' @param sig_label label significant variables
#' @param color colors used for ggplots.color=c("#FF3300","#FF6600","#FFCC00","#99CC00","#0066CC","#660099")
#' @param fixed_limit whether the limit should be fixed or not T, F
#' @param max_overlap maximum overlap for labels
#' @param FC_range significant fold change range
#'
#' @return Volcano Plot
#' @export
#'
#' @examples
#' data(Data)
#' Test<-Allstats(Data)
#' Volcano(Test,asterisk = "t-test")
Volcano<-function(data,
                  asterisk = "t-test",
                  reverse = F,
                  fig_width = NA,
                  FC_log = 2,
                  pval_log = 10,
                  fig_height = NA,
                  dotsize=2,
                  fixed_limit = F,
                  x_limit = c(-15,15),
                  y_limit = c(-2,6),
                  pval_intercept = 0.05,
                  max_overlap = 20,
                  FC_range = c(-1.5,1.5),
                  sig_label = T,
                  color = c("blue","black", "red")){
  ifelse(!dir.exists(file.path(getwd(), "Volcano")), dir.create(file.path(getwd(), "Volcano")), FALSE)
  data[["Data_renamed"]]<-data[["Data_renamed"]] %>% plyr::mutate(
    ZZZZ =data[["Data_renamed"]][,2] )
  data[["Data_renamed"]]<-data[["Data_renamed"]][,c(-1,-2)]
  data[["Data_renamed"]][,1:(ncol(data[["Data_renamed"]])-1)]<-sapply(data[["Data_renamed"]][,1:(ncol(data[["Data_renamed"]])-1)],
                                                                      function(x) as.numeric(x))
  NAMES<-colnames(data[["Data"]][3:ncol(data[["Data"]])])
  data[["Data_renamed"]]$ZZZZ<-as.factor(data[["Data_renamed"]]$ZZZZ)
  Comb<-gtools::combinations(length(unique(data[["Data"]]$Group)),2,unique(data[["Data"]]$Group))

  {if (asterisk =="Scheffe"){
    p_val_data<-data[["Anova_PostHoc"]]}
    else if (asterisk =="t-test"){
      p_val_data<-data[["t_test"]]}
    else if (asterisk =="u-test"){
      p_val_data<-data[["u_test"]]}
    else if (asterisk =="Dunn"){
      p_val_data<-data[["Dunn"]]}
    else {
      print("Wrong asterisk input must be one of (Dunn,Scheffe,u_test,t_test)")}}

  {if(length(unique(data[["Data"]]$Group)) !=2){
    for(a in 1:nrow(p_val_data)){
      for(b in 1:ncol(p_val_data)){
        if(is.nan(p_val_data[a,b])==T){
          p_val_data[a,b]<-1
        }
      }
    }}
    else if(length(unique(data[["Data"]]$Group))==2){
      for(b in 1:length(p_val_data)){
        if(is.nan(p_val_data[b])==T){
          p_val_data[b]<-1
        }
      }
    }}

  for (a in 1: nrow(Comb)){
    {if(reverse ==T){
      G1<- Comb[a,2]
      G2<- Comb[a,1]
    }
      else if (reverse ==F){
        G1<- Comb[a,1]
        G2<- Comb[a,2]
      }
      data[["Data"]][,3:ncol(data[["Data"]])]<-apply(data[["Data"]][,3:ncol(data[["Data"]])],2,function(x)as.numeric(x))
      Part_dat_1<-data[["Data"]] %>% dplyr::filter(Group == G1)
      Part_dat_2<-data[["Data"]] %>% dplyr::filter(Group == G2)
      Volc_Dat<-matrix(nrow = (ncol(data[["Data"]])-2),ncol = 3)
    }
    for (i in 3:ncol(Part_dat_1)){
      Volc_Dat[(i-2),1]<- log(mean(Part_dat_1[,i])/mean(Part_dat_2[,i]),base = FC_log) }
    if(length(unique(data[["Data"]]$Group))>2)
    {if(asterisk == "Scheffe"){
      {if (paste0(Comb[a,1],"-",Comb[a,2],"___","ANO_posthoc")%in% colnames(p_val_data)==T){
        part_p_val<-p_val_data[,paste0(Comb[a,1],"-",Comb[a,2],"___","ANO_posthoc")]
      }
        else {part_p_val<-p_val_data[,paste0(Comb[a,2],"-",Comb[a,1],"___","ANO_posthoc")]}}
    }
      else if(asterisk == "Dunn"){
        {if (paste0(Comb[a,1]," - ",Comb[a,2],"___","Kru_posthoc(Dunn)")%in% colnames(p_val_data)==T){
          part_p_val<-p_val_data[,paste0(Comb[a,1]," - ",Comb[a,2],"___","Kru_posthoc(Dunn)")]
        }
          else {part_p_val<-p_val_data[,paste0(Comb[a,2]," - ",Comb[a,1],"___","Kru_posthoc(Dunn)")]}}
      }
      else{
        {if (paste0(Comb[a,1],"-",Comb[a,2],"___",asterisk)%in% colnames(p_val_data)==T){
          part_p_val<-p_val_data[,paste0(Comb[a,1],"-",Comb[a,2],"___",asterisk)]
        }
          else {part_p_val<-p_val_data[,paste0(Comb[a,2],"-",Comb[a,1],"___",asterisk)]}}

      }}
    else  if(length(unique(data[["Data"]]$Group))==2){
      {part_p_val<-p_val_data}
    }

    Volc_Dat[,2]<-part_p_val
    rownames(Volc_Dat)<-colnames(data[["Data"]])[3:ncol(data[["Data"]])]
    Volc_Dat<-as.data.frame(Volc_Dat)
    colnames(Volc_Dat)<-c("log2FoldChange","pvalue","Direction")
    Volc_Dat$Direction <- "Not Significant"
    Volc_Dat$Direction[Volc_Dat$log2FoldChange > FC_range[2] & Volc_Dat$pvalue < pval_intercept] <- "UP"
    Volc_Dat$Direction[Volc_Dat$log2FoldChange < FC_range[1] & Volc_Dat$pvalue < pval_intercept] <- "DOWN"
    KK<-Volc_Dat %>% dplyr::filter(Direction!="Not Significant")
    sem_res<-ggplot2::ggplot(data=Volc_Dat, ggplot2::aes(x=log2FoldChange, y=-log(pvalue,base = pval_log), col=Direction)) +
      ggplot2::geom_point(size= dotsize) +
      ggplot2::labs(title=paste0(G1," divided by ",G2," Volcano plot"),
                    x =eval(parse(text = paste0("expression(Log[",FC_log,"]~(Fold~Change))"))),
                    y = eval(parse(text = paste0("expression(Log[",pval_log,"]~(P-value))"))))+
      ggplot2::theme_classic()+ ggplot2::theme(plot.title=ggplot2::element_text(size=25,
                                                                                color="Black",
                                                                                hjust=0.5,
                                                                                lineheight=1.2),  # title
                                               plot.subtitle=ggplot2::element_text(size=15,
                                                                                   hjust=0.5),  # subtitle
                                               plot.caption=ggplot2::element_text(size=15),  # caption
                                               axis.title.x=ggplot2::element_text(vjust=0,
                                                                                  size=15),  # X axis title
                                               axis.title.y=ggplot2::element_text(size=15),  # Y axis title
                                               axis.text.x=ggplot2::element_text(size=10,
                                                                                 vjust=.5),  # X axis text
                                               axis.text.y=ggplot2::element_text(size=10))+
      ggplot2::scale_color_manual(values=color[c("DOWN","Not Significant","UP")%in% unique(Volc_Dat$Direction)]) +
      ggplot2::geom_vline(xintercept=c(0), col=c("black")) +
      ggplot2::geom_vline(xintercept=FC_range[1],linetype = 'dashed', col=c("grey"))+
      ggplot2::geom_vline(xintercept=FC_range[2],linetype = 'dashed', col=c("grey"))+
      ggplot2::geom_hline(yintercept=-log(pval_intercept,base = pval_log), col="black")
    {if(fixed_limit ==T){
      sem_res<-sem_res+ggplot2::xlim(x_limit[1], x_limit[2])+
        ggplot2::ylim(y_limit[1],y_limit[2])}
      else {sem_res<-sem_res}}

    {if (sig_label== T){
      sem_res<- sem_res+ggrepel::geom_text_repel(size=3,show.legend = F,
                                                 ggplot2::aes(label =ifelse(Volc_Dat$Direction!="Not Significant",
                                                                            rownames(Volc_Dat),"")),
                                                 max.overlaps = max_overlap)
      sem_res
    }
      else {sem_res}}

    ggplot2::ggsave(filename = paste(G1, "divided by", G2, "Volcano plot.png",collapse = ""),
                    path=paste0(getwd(),"/Volcano"),
                    width = fig_width,
                    height = fig_height)

  }
}

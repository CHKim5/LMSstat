#' Automatically save boxplots for all metabolites including asterisks for significance.
#'
#' @param data List inheriting from Allstats
#' @param asterisk Choose asterisk to plot ("Dunn","Scheffe","u_test","t_test")
#' @param significant_variable_only show significant variable only (T,F)
#' @param color colors used for ggplots.
#' @param legend_position legend position "none","left","right","bottom","top"
#' @param order order of the groups; order = c("A","B","C")
#' @param tip_length significance tip length
#' @param label_size significance label size
#' @param step_increase significance step increase
#' @param width box width
#' @param fig_width figure size
#' @param fig_height figure size
#'
#' @return ggboxplot
#' @export
#'
#' @examples D<-as.data.frame(cbind(LETTERS,rep(1,26),c(rep("A",7),rep("B",6),rep("C",13)),runif(26),runif(26)))
#' Datas<-rbind(c("Sample","Multilevel","Group","some metabolite1","some metabolite2"),D)
#' Test<-Allstats(Datas)
#' AS_boxplot(Test,asterisk = "Dunn",significant_variable_only = F,color = c("#FF3300","#FF6600","#FFCC00","#99CC00","#0066CC","#660099") )
AS_boxplot<-function(data,
                     asterisk="t_test",
                     significant_variable_only=F,
                     color=c("#FF3300","#FF6600","#FFCC00","#99CC00","#0066CC","#660099"),
                     legend_position = "none",
                     order = NULL,
                     tip_length = 0.01,
                     label_size = 2.88,
                     step_increase = 0.05,
                     width = 0.3,
                     fig_width = NA,
                     fig_height = NA){
  {#Summary
    ###Plot_data_prep###
    ifelse(!dir.exists(file.path(getwd(), "barplot")), dir.create(file.path(getwd(), "boxplot")), FALSE)
    data[["Data_renamed"]]<-data[["Data_renamed"]] %>% plyr::mutate(
      ZZZZ =data[["Data_renamed"]][,2] )
    data[["Data_renamed"]]<-data[["Data_renamed"]][,c(-1,-2)]
    data[["Data_renamed"]][,1:(ncol(data[["Data_renamed"]])-1)]<-sapply(data[["Data_renamed"]][,1:(ncol(data[["Data_renamed"]])-1)],
                                                    function(x) as.numeric(x))
    NAMES<-colnames(data[["Data"]][3:ncol(data[["Data"]])])
    data[["Data_renamed"]]$ZZZZ<-as.factor(data[["Data_renamed"]]$ZZZZ)
    colnames(data[["Data_renamed"]])[ncol(data[["Data_renamed"]])]<-"Group"
    Comb<-gtools::combinations(length(unique(data[["Data"]]$Group)),2,unique(data[["Data"]]$Group))
    ckey<-color[1:length(unique(data[["Data"]]$Group))]

    if (asterisk =="Scheffe"){
      p_val_data<-data[["Anova_PostHoc"]]}
    else if (asterisk =="t_test"){
      p_val_data<-data[["t_test"]]}
    else if (asterisk =="u_test"){
      p_val_data<-data[["u_test"]]}
    else if (asterisk =="Dunn"){
      p_val_data<-data[["Dunn"]]}
    else {
      print("Wrong asterisk input must be one of (Dunn,Scheffe,u_test,t_test)")}}
  ###Plots###
  suppressWarnings(
    if (length(unique(data[["Data_renamed"]][["Group"]]))!=2) {
    for (number in 1:nrow(p_val_data)){
      if (number == round(nrow(p_val_data)/4)){print("25% complete")}
      if (number == round(nrow(p_val_data)/3)){print("33% complete")}
      if (number == round(nrow(p_val_data)/2)){print("50% complete")}
      if (number == round(2*nrow(p_val_data)/3)){print("66% complete")}
      if (number == round(3*nrow(p_val_data)/4)){print("75% complete")}
      if (number == nrow(p_val_data)){print("100% complete")}
      stat.test <- as.data.frame(matrix(data = NA, nrow = nrow(Comb),
                                        ncol = 2))
      stat.test[, 1:2] <- Comb
      if (asterisk == "Dunn"){
        rownames(stat.test)<-paste0(stat.test[,1]," - ",stat.test[,2])}
      else{
        rownames(stat.test)<-paste0(stat.test[,1],"-",stat.test[,2])}
      colnames(p_val_data)<-as.data.frame(strsplit(colnames(p_val_data),"___"))[1,]
      t1<-as.data.frame(t(as.data.frame(p_val_data)[number,]))
      t2<-stat.test
      t3<-merge(t2,t1,by=0)
      stat.test<-t3[,2:4]
      colnames(stat.test) <- c("group1", "group2",
                               "p")
      stat.test<-stat.test %>% plyr::mutate(
        p.adj.signif = case_when(
          p > 0.05                   ~ 'NS',
          p <= 0.05  &p > 0.01                  ~ '*',
          p <= 0.01                   ~ '**'
        )
      )
      stat.test<-stat.test[stat.test$p.adj.signif != "NS",]
      if (length(stat.test>4)){
        ggpubr::ggboxplot(data[["Data_renamed"]],
                          x = "Group",
                          y = colnames(data[["Data_renamed"]])[number],
                          color = "Group",
                          palette = ckey,
                          add = "jitter",
                          order = order,
                          width = width)+
          ggplot2::scale_y_continuous(label=ecoflux::scientific_10x)+
          ggplot2::labs(title=NAMES[number],x =NULL, y = "Intensity")+
          ggplot2::theme_classic()+ ggplot2::theme(plot.title=element_text(size=15,
                                                         face="bold",
                                                         color="Black",
                                                         hjust=0.5,
                                                         lineheight=1.2),  # title
                                 plot.subtitle=element_text(size=15,
                                                            hjust=0.5),  # subtitle
                                 plot.caption=element_text(size=15),  # caption
                                 axis.title.x=element_text(vjust=10,
                                                           size=15),  # X axis title
                                 axis.title.y=element_text(size=12),  # Y axis title
                                 axis.text.x=element_text(size=10,
                                                          vjust=.5),  # X axis text
                                 axis.text.y=element_text(size=10),
                                 legend.position = legend_position)+
          ggpubr::stat_pvalue_manual(stat.test,
                                     y.position = 1.05 * max(data[["Data_renamed"]][,number]),
                                     step.increase = step_increase,
                                     label.size = label_size,
                                     tip.length = tip_length,
                                     label = "p.adj.signif",
                                     size = 3.5,
                                     vjust = 0.05)
        ggplot2::ggsave(filename = paste(NAMES[number],"boxplot.png",collapse = ""),
               path=paste0(getwd(),"/boxplot"),
               width = fig_width,
               height = fig_height)}
      else if(significant_variable_only == F) {
        ggpubr::ggboxplot(data[["Data_renamed"]],
                          x = "Group",
                          y = colnames(data[["Data_renamed"]])[number],
                          color = "Group",
                          palette = ckey,
                          add = "jitter",
                          order = order,
                          width = width)+
          ggplot2::scale_y_continuous(label=ecoflux::scientific_10x)+
          ggplot2::labs(title=NAMES[number],x =NULL, y = "Intensity")+
          ggplot2::theme_classic()+ ggplot2::theme(plot.title=element_text(size=15,
                                                                           face="bold",
                                                                           color="Black",
                                                                           hjust=0.5,
                                                                           lineheight=1.2),  # title
                                                   plot.subtitle=element_text(size=15,
                                                                              hjust=0.5),  # subtitle
                                                   plot.caption=element_text(size=15),  # caption
                                                   axis.title.x=element_text(vjust=10,
                                                                             size=15),  # X axis title
                                                   axis.title.y=element_text(size=12),  # Y axis title
                                                   axis.text.x=element_text(size=10,
                                                                            vjust=.5),  # X axis text
                                          axis.text.y=element_text(size=10),
                                          legend.position = legend_position)

        ggplot2::ggsave(filename = paste(NAMES[number],"boxplot.png",collapse = ""),
               path=paste0(getwd(),"/boxplot"),
               width = fig_width,
               height = fig_height)}}}
    else if (length(unique(data[["Data_renamed"]][["Group"]]))==2) {
      p_val_data<-as.data.frame(p_val_data)
      if (asterisk=="t_test"){
        colnames(p_val_data)<-colnames(data[["Result"]])[1]
      }
      else if(asterisk == "u_test"){
        colnames(p_val_data)<-colnames(data[["Result"]])[2]
      }

      for (number in 1:nrow(p_val_data)){
        if (number == round(nrow(p_val_data)/4)){print("25% complete")}
        if (number == round(nrow(p_val_data)/3)){print("33% complete")}
        if (number == round(nrow(p_val_data)/2)){print("50% complete")}
        if (number == round(2*nrow(p_val_data)/3)){print("66% complete")}
        if (number == round(3*nrow(p_val_data)/4)){print("75% complete")}
        if (number == nrow(p_val_data)){print("100% complete")}
        df <- data_summary(data[["Data_renamed"]],
                           varname = colnames(data[["Data_renamed"]])[number],
                           groupnames = c("Group"))
        df <- df %>% plyr::mutate(ymax = len + SEM)
        colnames(df)[1] <- "Group"
        stat.test <- as.data.frame(matrix(data = NA, nrow = nrow(Comb),
                                          ncol = 2))
        stat.test[, 1:2] <- Comb
        if (asterisk == "Dunn") {
          rownames(stat.test) <- paste0(stat.test[, 1],
                                        " - ", stat.test[, 2])
        }
        else {
          rownames(stat.test) <- paste0(stat.test[, 1],
                                        "-", stat.test[, 2])
        }
        colnames(p_val_data) <- as.data.frame(strsplit(colnames(p_val_data),
                                                       "___"))[1, ]
        t1 <- as.data.frame(t(as.data.frame(p_val_data[number,
        ])))
        rownames(t1)<-colnames(p_val_data)

        t2 <- stat.test
        t3 <- merge(t2, t1, by = 0)
        stat.test <- t3[, 2:4]
        colnames(stat.test) <- c("group1", "group2",
                                 "p")
        stat.test <- stat.test %>% plyr::mutate(p.adj.signif = case_when(p >
                                                                           0.05 ~ "NS", p <= 0.05 & p > 0.01 ~ "*",
                                                                         p <= 0.01 ~ "**"))
        stat.test <- stat.test[stat.test$p.adj.signif !=
                                 "NS", ]
        if (length(stat.test>4)){
          ggpubr::ggboxplot(data[["Data_renamed"]],
                            x = "Group",
                            y = colnames(data[["Data_renamed"]])[number],
                            color = "Group",
                            palette = ckey,
                            add = "jitter",
                            order = order,
                            width = width)+
            ggplot2::scale_y_continuous(label=ecoflux::scientific_10x)+
            ggplot2::labs(title=NAMES[number],x =NULL, y = "Intensity")+
            ggplot2::theme_classic()+ ggplot2::theme(plot.title=element_text(size=15,
                                                                    face="bold",
                                                                    color="Black",
                                                                    hjust=0.5,
                                                                    lineheight=1.2),  # title
                                            plot.subtitle=element_text(size=15,
                                                                       hjust=0.5),  # subtitle
                                            plot.caption=element_text(size=15),  # caption
                                            axis.title.x=element_text(vjust=10,
                                                                      size=15),  # X axis title
                                            axis.title.y=element_text(size=12),  # Y axis title
                                            axis.text.x=element_text(size=10,
                                                                     vjust=.5),  # X axis text
                                            axis.text.y=element_text(size=10),
                                            legend.position = legend_position,
                                            step.increase = step_increase,
                                            label.size = label_size,
                                            tip.length = tip_length)+
            ggpubr::stat_pvalue_manual(stat.test, y.position = 1.05 * max(data[["Data_renamed"]][,number]), step.increase = 0.1,label = "p.adj.signif",size = 3.5,vjust = 0.05)
          ggplot2::ggsave(filename = paste(NAMES[number],"boxplot.png",collapse = ""),
                 path=paste0(getwd(),"/boxplot"),
                 width = fig_width,
                 height = fig_height)}
        else if(significant_variable_only == F) {
          ggpubr::ggboxplot(data[["Data_renamed"]],
                            x = "Group",
                            y = colnames(data[["Data_renamed"]])[number],
                            color = "Group",
                            palette = ckey,
                            add = "jitter",
                            order = order,
                            width = width)+
            ggplot2::scale_y_continuous(label=ecoflux::scientific_10x)+
            ggplot2::labs(title=NAMES[number],x =NULL, y = "Intensity")+
            ggplot2::theme_classic()+ ggplot2::theme(plot.title=element_text(size=15,
                                                           face="bold",
                                                           color="Black",
                                                           hjust=0.5,
                                                           lineheight=1.2),  # title
                                   plot.subtitle=element_text(size=15,
                                                              hjust=0.5),  # subtitle
                                   plot.caption=element_text(size=15),  # caption
                                   axis.title.x=element_text(vjust=10,
                                                             size=15),  # X axis title
                                   axis.title.y=element_text(size=12),  # Y axis title
                                   axis.text.x=element_text(size=10,
                                                            vjust=.5),  # X axis text
                                   axis.text.y=element_text(size=10),
                                   legend.position = legend_position)

          ggplot2::ggsave(filename = paste(NAMES[number],"boxplot.png",collapse = ""),
                 path=paste0(getwd(),"/boxplot"),
                 width = fig_width,
                 height = fig_height)}}}

    )
}

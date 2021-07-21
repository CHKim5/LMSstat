#' Heatmap with top annotation bar
#'
#' @param data data inheriting from D_tran
#' @param col colors for heatmap
#' @param col_lim color boundaries c(-3, 0, 3)
#' @param reverse Reverse column and rows
#' @param distance "pearson", "manhattan","euclidean","spearman","kendall"
#' @param rownames rownames
#' @param colnames colnames
#' @param Hsize Width & Height c(a,b)
#' @param g_legend Annotation legend title
#' @param h_legend Heatmap legend title
#' @param T_size Title text size
#' @param R_size row text size
#' @param C_size column text size
#' @param Gcol Color for top_annotation bar c("ASD" = "black","HGH"="red","LAC"="blue","LUE" ="grey","SDF" = "yellow","WEI"="green")
#' @param Title title
#' @param dend_h dendrite height
#'
#' @return Heatmap
#' @export
#'
#' @examples data(Data)
#' data<-D_tran(Data,param="Auto",save=F)
#' AS_heatmap(data)
#' dev.off()
AS_heatmap<-function(data,col = c("green", "white", "red"),
                     col_lim = c(-3,0,3),
                     reverse = F,
                     distance = "euclidean",
                     rownames =TRUE,
                     colnames = FALSE,
                     Hsize = c(12,6),
                     g_legend = "Group",
                     h_legend = "Intensity",
                     T_size = 10,
                     R_size = 7,
                     C_size = 5,
                     Gcol = c("ASD" = "black","HGH"="red","LAC"="blue","LUE" ="grey","SDF" = "yellow","WEI"="green"),
                     Title = "Sample Heatmap",
                     dend_h = 0.5){
m_for_heatmap<-as.matrix(data[,2:ncol(data)])
colors = circlize::colorRamp2(col_lim, col)
if (reverse == T){
  m_for_heatmap<-t(m_for_heatmap)
  row_log =TRUE
  col_log =FALSE} else if (reverse ==F){
    row_log =FALSE
    col_log =TRUE}

Gcol_D<-as.data.frame(Gcol)
ha = ComplexHeatmap::HeatmapAnnotation(
  Group = data$Group,col = list(Group = Gcol),
  annotation_legend_param = F,show_annotation_name = F,show_legend = F)
G_I = ComplexHeatmap::Legend(labels = rownames(Gcol_D), title = g_legend, legend_gp = grid::gpar(fill = Gcol_D[,1]))
M_I = ComplexHeatmap::Legend(col_fun = colors, title = h_legend)
pd = ComplexHeatmap::packLegend(list = list(G_I,M_I))
kt<-ComplexHeatmap::Heatmap(t(m_for_heatmap),
                        col=colors,
                        column_title = Title,
                        column_title_gp = grid::gpar(fontsize = T_size),
                        name = "legend",
                        row_names_gp = grid::gpar(fontsize = R_size),
                        column_names_gp = grid::gpar(fontsize = C_size),
                        cluster_rows = row_log,
                        cluster_columns = col_log,
                        border_gp = grid::gpar(col = "black", lty = 1),
                        column_dend_height = grid::unit(dend_h, "cm"),
                        clustering_distance_rows = distance,
                        clustering_distance_columns = distance,
                        show_column_names = colnames,
                        show_row_names = rownames,
                        width = grid::unit(Hsize[1], "cm"), height = grid::unit(Hsize[2], "cm"),
                        top_annotation = ha,
                        show_heatmap_legend = F)
pdf(paste0(Title,"Heatmap.pdf"))
ComplexHeatmap::draw(kt, annotation_legend_list = pd)
}

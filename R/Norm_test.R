#' Normality test with Shaprio_Wilk_test
#'
#' @param Data input data file
#'
#' @return dataframe with shapiro wilk test
#' @export
#'
#' @examples Result<-Norm_test(Data)
#' write.csv(Result,"Normality_test_Result.csv")
#' 
Norm_test<-function(Data){
colnames(Data)<-Data[1,]; Data<-Data[-1,]
Data[,4:ncol(Data)]<-sapply(Data[,4:ncol(Data)],function(x) as.numeric(x))
Data<-as.data.frame(sapply(Data[,4:ncol(Data)], function(x)shapiro.test(x)[["p.value"]]))
colnames(Data)<-"Shapiro_Wilk_p_val"
Data<-Data %>% mutate(Normality = case_when(Shapiro_Wilk_p_val<=0.05~"Non_parametric",
                                            Shapiro_Wilk_p_val>0.05~"Parametric"))
Data}

#' Automatically processes T-test, U-test, Anova, Scheffe(Anova Post-Hoc), Krukal Wallis, Dunn-test(BH adjusted,(Kurkal Wallis Post-Hoc)) while allowing adjustment of FDR
#'
#' @param Data csv file with Header as False First column with Sample Second column with Multilevel(Mixomics) so that it can be compatible with other multivariate statistics Third column with Group information. Rest of the columns are the metabolites to be tested.
#' @param Adjust_p_value Set True if FDR adjustments are to be made. If not set False
#' @param Adjust_method adjustment methods frequently used. "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"
#'
#' @return List including Result Matrix of p-values, converted datas.
#' @export
#'
#' @examples data(Data)
#' Result<-Allstats(Data)
#'
#'
Allstats<-function (Data, Adjust_p_value = T, Adjust_method = "BH")
{

  LETTERS702 <- c(sapply(LETTERS, function(x) paste0(x, LETTERS)))
  LETTERS37232 <- c(LETTERS, LETTERS702, sapply(LETTERS, function(x) paste0(x,
                                                                            LETTERS702)))
  LETTERS210729 <- c(LETTERS, LETTERS702,LETTERS37232, sapply(LETTERS, function(x) paste0(x,
                                                                                          LETTERS37232)))
  LETTERS210729 <- LETTERS210729[-365]
  LETTERS210729<-LETTERS210729[-10205]
  LETTERS210729<-LETTERS210729[-267101]
  colnames(Data) <- Data[1, ]
  Data <- Data[-1, -2]
  Data<-Data %>% dplyr::arrange(Group)
  Data_renamed <- Data
  nmet <- ncol(Data) - 2
  colnames(Data_renamed) <- c(colnames(Data[1:2]), LETTERS210729[1:nmet])
  rownames(Data_renamed) <- Data[, 1]
  Data_renamed_raw <- Data_renamed[, -c(1, 2)]
  Data_renamed_raw <- apply(Data_renamed_raw, 2, as.numeric)
  Data_final <- cbind(Data[, 1:2], Data_renamed_raw)
  if (length(unique(Data_final$Group)) == 2) {
    suppressWarnings({
      for (Gnum in 1:length(unique(Data$Group))) {
        assign(paste("Data", unique(Data$Group)[Gnum],
                     sep = "_"), dplyr::filter(Data_final, Group ==
                                          unique(Data_final$Group)[Gnum]))
      }
      Result <- matrix(data = NA, nrow = (ncol(Data_final) -
                                            2), ncol = (2 * choose(length(unique(Data$Group)),
                                                                   2)))
      Rounder <- 1
      for (Gnum_1 in 1:length(unique(Data$Group))) {
        for (Gnum_2 in 1:length(unique(Data$Group))) {
          if (Gnum_2 <= Gnum_1) {
            next
          }
          else {
            for (met in 1:(ncol(Data) - 2)) {
              assign(paste(as.character(LETTERS210729[met]),
                           unique(Data$Group)[Gnum_1], unique(Data$Group)[Gnum_2],
                           "T_test", sep = "_"), t.test(x = (eval(parse(text = paste("Data",
                                                                                     unique(Data$Group)[Gnum_1], sep = "_")))[,
                                                                                                                              eval(parse(text = "LETTERS210729[met]"))]),
                                                        y = (eval(parse(text = paste("Data",
                                                                                     unique(Data$Group)[Gnum_2], sep = "_")))[,
                                                                                                                              eval(parse(text = "LETTERS210729[met]"))]))[["p.value"]])
            }
          }
        }
      }
      for (C in 1:length(unique(Data$Group))) {
        for (H in 1:length(unique(Data$Group))) {
          if (H <= C) {
            next
          }
          else {
            if (is.na(Result[1, Rounder]) == T) {
              for (met_1 in 1:(ncol(Data_final) - 2)) {
                Result[met_1, Rounder] <- eval(parse(text = paste(LETTERS210729[met_1],
                                                                  unique(Data$Group)[C], unique(Data$Group)[H],
                                                                  "T_test", sep = "_")))
              }
              Rounder <- Rounder + 1
            }
          }
        }
      }
      print("T-test has finished")
      for (Gnum_1 in 1:length(unique(Data$Group))) {
        for (Gnum_2 in 1:length(unique(Data$Group))) {
          if (Gnum_2 <= Gnum_1) {
            next
          }
          else {
            for (met in 1:(ncol(Data) - 2)) {
              assign(paste(as.character(LETTERS210729[met]),
                           unique(Data$Group)[Gnum_1], unique(Data$Group)[Gnum_2],
                           "U_test", sep = "_"), wilcox.test(x = (eval(parse(text = paste("Data",
                                                                                          unique(Data$Group)[Gnum_1], sep = "_")))[,
                                                                                                                                   eval(parse(text = "LETTERS210729[met]"))]),
                                                             y = (eval(parse(text = paste("Data",
                                                                                          unique(Data$Group)[Gnum_2], sep = "_")))[,
                                                                                                                                   eval(parse(text = "LETTERS210729[met]"))]))[["p.value"]])
            }
          }
        }
      }
      for (C in 1:length(unique(Data$Group))) {
        for (H in 1:length(unique(Data$Group))) {
          if (H <= C) {
            next
          }
          else {
            if (is.na(Result[1, Rounder]) == T) {
              for (met_1 in 1:(ncol(Data_final) - 2)) {
                Result[met_1, Rounder] <- eval(parse(text = paste(LETTERS210729[met_1],
                                                                  unique(Data$Group)[C], unique(Data$Group)[H],
                                                                  "U_test", sep = "_")))
              }
              Rounder <- Rounder + 1
            }
          }
        }
      }
      print("U-test has finished")
      Names <- NULL
      for (C in 1:length(unique(Data$Group))) {
        for (H in 1:length(unique(Data$Group))) {
          if (H <= C) {
            next
          }
          else {
            Names <- rbind(Names, paste(unique(Data$Group)[C],
                                        unique(Data$Group)[H], sep = "-"))
          }
        }
      }
      rownames(Result) <- colnames(Data)[3:(ncol(Data_final))]
      colnames(Result) <- c(paste(Names[, 1], "t-test",
                                  sep = "___"), paste(Names[, 1], "u-test",
                                                      sep = "___"))
      if (Adjust_p_value == T) {
        print("###########################################")
        print(paste0("adjusted according to the ",
                     Adjust_method, " method"))
        print("###########################################")
        Result <- apply(Result, 2, function(x) {
          p.adjust(x, method = Adjust_method)
        })
      }
      else {
        print("###########################################")
        print("p_value not adjusted")
        print("###########################################")
      }
      rm(list = setdiff(ls(), c("Data_renamed", "Data",
                                "Result", "LETTERS210729", "event",
                                "P_hoc", "Colors", "significant_variable_only")))
      print("statistical test has finished")
      Result_T <- Result[, 1]
      Result_U <- Result[, 2]
      print("subsets have been made")
      Final <- list()
      Final$Data <- Data
      Final$Data_renamed <- Data_renamed
      Final$Result <- Result
      Final$t_test <- Result_T
      Final$u_test <- Result_U
      Final
    })
  }
  else if (length(unique(Data_final$Group)) > 2) {
    suppressWarnings({
      for (Gnum in 1:length(unique(Data$Group))) {
        assign(paste("Data", unique(Data$Group)[Gnum],
                     sep = "_"), filter(Data_final, Group ==
                                          unique(Data_final$Group)[Gnum]))
      }
      Result <- matrix(data = NA, nrow = (ncol(Data_final) -
                                            2), ncol = (2 + 4 * choose(length(unique(Data$Group)),
                                                                       2)))
      Rounder <- 1
      for (Gnum_1 in 1:length(unique(Data$Group))) {
        for (Gnum_2 in 1:length(unique(Data$Group))) {
          if (Gnum_2 <= Gnum_1) {
            next
          }
          else {
            for (met in 1:(ncol(Data) - 2)) {
              assign(paste(as.character(LETTERS210729[met]),
                           unique(Data$Group)[Gnum_1], unique(Data$Group)[Gnum_2],
                           "T_test", sep = "_"), t.test(x = (eval(parse(text = paste("Data",
                                                                                     unique(Data$Group)[Gnum_1], sep = "_")))[,
                                                                                                                              eval(parse(text = "LETTERS210729[met]"))]),
                                                        y = (eval(parse(text = paste("Data",
                                                                                     unique(Data$Group)[Gnum_2], sep = "_")))[,
                                                                                                                              eval(parse(text = "LETTERS210729[met]"))]))[["p.value"]])
            }
          }
        }
      }
      for (C in 1:length(unique(Data$Group))) {
        for (H in 1:length(unique(Data$Group))) {
          if (H <= C) {
            next
          }
          else {
            if (is.na(Result[1, Rounder]) == T) {
              for (met_1 in 1:(ncol(Data_final) - 2)) {
                Result[met_1, Rounder] <- eval(parse(text = paste(LETTERS210729[met_1],
                                                                  unique(Data$Group)[C], unique(Data$Group)[H],
                                                                  "T_test", sep = "_")))
              }
              Rounder <- Rounder + 1
            }
          }
        }
      }
      print("T-test has finished")
      for (Gnum_1 in 1:length(unique(Data$Group))) {
        for (Gnum_2 in 1:length(unique(Data$Group))) {
          if (Gnum_2 <= Gnum_1) {
            next
          }
          else {
            for (met in 1:(ncol(Data) - 2)) {
              assign(paste(as.character(LETTERS210729[met]),
                           unique(Data$Group)[Gnum_1], unique(Data$Group)[Gnum_2],
                           "U_test", sep = "_"), wilcox.test(x = (eval(parse(text = paste("Data",
                                                                                          unique(Data$Group)[Gnum_1], sep = "_")))[,
                                                                                                                                   eval(parse(text = "LETTERS210729[met]"))]),
                                                             y = (eval(parse(text = paste("Data",
                                                                                          unique(Data$Group)[Gnum_2], sep = "_")))[,
                                                                                                                                   eval(parse(text = "LETTERS210729[met]"))]))[["p.value"]])
            }
          }
        }
      }
      for (C in 1:length(unique(Data$Group))) {
        for (H in 1:length(unique(Data$Group))) {
          if (H <= C) {
            next
          }
          else {
            if (is.na(Result[1, Rounder]) == T) {
              for (met_1 in 1:(ncol(Data_final) - 2)) {
                Result[met_1, Rounder] <- eval(parse(text = paste(LETTERS210729[met_1],
                                                                  unique(Data$Group)[C], unique(Data$Group)[H],
                                                                  "U_test", sep = "_")))
              }
              Rounder <- Rounder + 1
            }
          }
        }
      }
      print("U-test has finished")
      for (met in 1:(ncol(Data) - 2)) {
        Ano_In <- aov(eval(parse(text = LETTERS210729[met])) ~
                        Group, data = Data_final)
        assign(paste(as.character(LETTERS210729[met]),
                     "Anova", sep = "_"), summary(Ano_In)[[1]][["Pr(>F)"]][1])
      }
      if (is.na(Result[1, Rounder]) == T) {
        for (met_1 in 1:(ncol(Data_final) - 2)) {
          Result[met_1, Rounder] <- eval(parse(text = paste(LETTERS210729[met_1],
                                                            "Anova", sep = "_")))
        }
        Rounder <- Rounder + 1
      }
      else {
        print("Anova went wrong")
      }
      for (met in 1:(ncol(Data) - 2)) {
        assign(paste(LETTERS210729[met], "Anova_Post_Hoc",
                     sep = "_"), DescTools::PostHocTest(aov(eval(parse(text = LETTERS210729[met])) ~
                                                              Group, data = Data_final), method = "scheffe"))
      }
      for (Ano_numb in 1:nrow(eval(parse(text = paste(LETTERS210729[met],
                                                      "Anova_Post_Hoc", sep = "_")))[["Group"]])) {
        if (is.na(Result[1, Rounder]) == T) {
          for (met_1 in 1:(ncol(Data_final) - 2)) {
            Result[met_1, Rounder] <- eval(parse(text = paste(LETTERS210729[met_1],
                                                              "Anova_Post_Hoc", sep = "_")))[["Group"]][Ano_numb,
                                                                                                        4]
          }
          Rounder <- Rounder + 1
        }
        else {
          print("Anova went wrong")
        }
      }
      print("Anova & PostHoc has finished")
      for (met in 1:(ncol(Data) - 2)) {
        Kru_In <- kruskal.test(eval(parse(text = LETTERS210729[met])) ~
                                 Group, data = Data_final)
        assign(paste(as.character(LETTERS210729[met]),
                     "Kruskal_Wallis", sep = "_"), Kru_In[["p.value"]])
      }
      if (is.na(Result[1, Rounder]) == T) {
        for (met_1 in 1:(ncol(Data_final) - 2)) {
          Result[met_1, Rounder] <- eval(parse(text = paste(LETTERS210729[met_1],
                                                            "Kruskal_Wallis", sep = "_")))
        }
        Rounder <- Rounder + 1
      }
      else {
        print("Kruskal went wrong")
      }
      for (met in 1:(ncol(Data) - 2)) {
        assign(paste(LETTERS210729[met], "Dunn_Post_Hoc",
                     sep = "_"), FSA::dunnTest(eval(parse(text = LETTERS210729[met])) ~
                                                 Group, data = Data_final, method = "bh"))
      }
      for (Kru_numb in 1:length(eval(parse(text = paste(LETTERS210729[met],
                                                        "Dunn_Post_Hoc", sep = "_")))[["res"]][["Comparison"]])) {
        if (is.na(Result[1, Rounder]) == T) {
          for (met_1 in 1:(ncol(Data_final) - 2)) {
            Result[met_1, Rounder] <- eval(parse(text = paste(LETTERS210729[met_1],
                                                              "Dunn_Post_Hoc", sep = "_")))[["res"]][["P.adj"]][Kru_numb]
          }
          Rounder <- Rounder + 1
        }
        else {
          print("Dunn went wrong")
        }
      }
      print("Kruskal Wallis & PostHoc has finished")
      Names <- NULL
      for (C in 1:length(unique(Data$Group))) {
        for (H in 1:length(unique(Data$Group))) {
          if (H <= C) {
            next
          }
          else {
            Names <- rbind(Names, paste(unique(Data$Group)[C],
                                        unique(Data$Group)[H], sep = "-"))
          }
        }
      }
      AN_Post_names <- NULL
      for (nrow in 1:nrow(eval(parse(text = paste(LETTERS210729[met],
                                                  "Anova_Post_Hoc", sep = "_")))[["Group"]])) {
        AN_Post_names <- rbind(AN_Post_names, rownames(eval(parse(text = paste(LETTERS210729[met],
                                                                               "Anova_Post_Hoc", sep = "_")))[["Group"]])[nrow])
      }
      DU_post_names <- NULL
      for (nrow in 1:length(eval(parse(text = paste(LETTERS210729[met],
                                                    "Dunn_Post_Hoc", sep = "_")))[["res"]][["Comparison"]])) {
        DU_post_names <- rbind(DU_post_names, (eval(parse(text = paste(LETTERS210729[met],
                                                                       "Dunn_Post_Hoc", sep = "_")))[["res"]][["Comparison"]])[nrow])
      }
      rownames(Result) <- colnames(Data)[3:(ncol(Data_final))]
      colnames(Result) <- c(paste(Names[, 1], "t-test",
                                  sep = "___"), paste(Names[, 1], "u-test",
                                                      sep = "___"), "Anova", paste(AN_Post_names,
                                                                                   "ANO_posthoc", sep = "___"), "Kruskal_Wallis",
                            paste(DU_post_names, "Kru_posthoc(Dunn)",
                                  sep = "___"))
      if (Adjust_p_value == T) {
        print("###########################################")
        print(paste0("adjusted according to the ",
                     Adjust_method, " method"))
        print("###########################################")
        Result <- apply(Result, 2, function(x) {
          p.adjust(x, method = Adjust_method)
        })
      }
      else {
        print("###########################################")
        print("p_value not adjusted")
        print("###########################################")
      }
      rm(list = setdiff(ls(), c("Data_renamed", "Data",
                                "Result", "LETTERS210729", "event",
                                "P_hoc", "Colors", "significant_variable_only")))
      print("statistical test has finished")
    })
    Result_T <- Result[, 1:choose(length(unique(Data$Group)),
                                  2)]
    Result_U <- Result[, (choose(length(unique(Data$Group)),
                                 2) + 1):(2 * choose(length(unique(Data$Group)), 2))]
    Result_Ano <- as.data.frame(Result[, (2 * choose(length(unique(Data$Group)),
                                                     2) + 1):(2 * choose(length(unique(Data$Group)), 2) +
                                                                1)])
    colnames(Result_Ano) <- colnames(Result)[(2 * choose(length(unique(Data$Group)),
                                                         2) + 1)]
    Result_Ano_P <- Result[, (2 * choose(length(unique(Data$Group)),
                                         2) + 2):(3 * choose(length(unique(Data$Group)), 2) +
                                                    1)]
    Result_Kru <- as.data.frame(Result[, (3 * choose(length(unique(Data$Group)),
                                                     2) + 2):(3 * choose(length(unique(Data$Group)), 2) +
                                                                2)])
    colnames(Result_Kru) <- colnames(Result)[(3 * choose(length(unique(Data$Group)),
                                                         2) + 2)]
    Result_Kru_P <- Result[, (3 * choose(length(unique(Data$Group)),
                                         2) + 3):(4 * choose(length(unique(Data$Group)), 2) +
                                                    2)]
    print("subsets have been made")
    Final <- list()
    Final$Data <- Data
    Final$Data_renamed <- Data_renamed
    Final$Result <- Result
    Final$Anova <- Result_Ano
    Final$Anova_PostHoc <- Result_Ano_P
    Final$KW <- Result_Kru
    Final$Dunn <- Result_Kru_P
    Final$t_test <- Result_T
    Final$u_test <- Result_U
    Final
  }
}

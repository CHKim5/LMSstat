# LMSstat
Automation of statistical test with an identical data input aiming to reduce arduous work searching for packages and changing data input.

The package includes

* Simple Statistics :u-test, t-test, post hocs of Anova and Kruskal Wallis with FDR adjusted values

* Bar, Box, Dot plots with significance (u-test, t-test, post hocs of Anova and Kruskal Wallis)

## Instructions

### Installation
install.packages("devtools")

devtools::install_github("CHKim5/LMSstat")

### Univariate statistics

data(Data) # Sample data 

Statfile<-Allstats(Data,Adjust_p_value = T, Adjust_method = "BH")

head(Statfile[["Result"]]) # includes all statistical results

write.csv(Statfile[["Result"]],"p_value_result.csv")  # Write csv with all the p-value included

### Plots

AS_boxplot(Statfile) # aesthetics can be adjusted   

AS_dotplot(Statfile)

AS_barplot(Statfile)

#### Adjustable parameters

* asterisk = "t_test"  
* significant_variable_only = F  # If set to TRUE, insignificant results will not be plotted
* color = c("#FF3300", "#FF6600", "#FFCC00", "#99CC00", "#0066CC", "#660099") # Colors for the plots
* legend_position = "none" #  "none","left","right","bottom","top"
* order = NULL # Order of the groups c("some_metabolite1","some_metabolite3","some_metabolite2")
* tip_length = 0.01 # significance tip length
* label_size = 2.88 # significance label size
* step_increase = 0.05 #significance step increase
* width = 0.3 # box width ; size = 3 # dot size
* fig_width = NA #figure size 
* fig_height = NA #figure size

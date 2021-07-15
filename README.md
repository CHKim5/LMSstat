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
#### Adjustable parameters
* Adjust_p_value = T # Set True if adjustment is needed
* Adjust_method = F # Adjustment methods frequently used. c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")

head(Statfile[["Result"]]) # includes all statistical results

write.csv(Statfile[["Result"]],"p_value_result.csv")  # Write csv with all the p-value included

### Plots

&emsp;&emsp;&emsp;&emsp;&emsp;**AS_boxplot(Statfile)**&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;**AS_dotplot(Statfile)**&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;**AS_barplot(Statfile)**
<p align="center">
<img src="https://user-images.githubusercontent.com/77651662/125714687-1908d6eb-b2bd-4e25-8ef0-62c24466c32a.png" width="225" height="250">
<img src="https://user-images.githubusercontent.com/77651662/125714704-d7dab67e-03c0-4e35-b86a-36723f7c63de.png" width="225" height="250">
<img src="https://user-images.githubusercontent.com/77651662/125715925-0878ec77-30bf-4859-8e56-316d98b6d520.jpg" width="350" height="250">
</p>

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

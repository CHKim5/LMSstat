# LMSstat
Automation of statistical test with an identical data input aiming to reduce arduous work searching for packages and changing data input.

The package includes

* Simple Statistics :u-test, t-test, post hocs of Anova and Kruskal Wallis with FDR adjusted values

* Bar, Box, Dot plots with significance (u-test, t-test, post hocs of Anova and Kruskal Wallis)

* PERMANOVA

## Instructions

### Installation

```
install.packages("devtools")
devtools::install_github("CHKim5/LMSstat")
```


### Basic structure of the Data
#### Used in

* Simple statistics
* Barplot, Boxplot, Dotplot  
* PERMANOVA

```
Data<-read.csv("statT.csv",header = F)
```
**provided within the package, can be called with data("Data")**

<p align="center">  
<img src="https://user-images.githubusercontent.com/77651662/125877205-e140d306-81d8-459f-8414-b8ef3bca63d7.PNG" width="750" height="400">
</p>
<p align="center">statT.csv</p>


#### Used in

* PERMANOVA

```
Classification<-read.csv("statT_G.csv",header = F)
```

 **provided within the package, can be called with data("Classification")**
<p align="center">
<img src="https://user-images.githubusercontent.com/77651662/125877154-d01ad8b0-25cd-448b-905d-749a6cc93552.PNG" width="500" height="400">
</p>
<p align="center">statT_G.csv</p>


### Univariate statistics

```
data(Data) # Sample data 

Statfile<-Allstats(Data,Adjust_p_value = T, Adjust_method = "BH")
```

#### Adjustable parameters
* Adjust_p_value = T # Set True if adjustment is needed
* Adjust_method = F # Adjustment methods frequently used. c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")

```
head(Statfile[["Result"]]) # includes all statistical results

write.csv(Statfile[["Result"]],"p_value_result.csv")  # Write csv with all the p-value included
```
### Plots

```
# Makes a subdirectory and saves boxplots for all the variables
AS_boxplot(Statfile,asterisk = "u_test") 

# Makes a subdirectory and saves dotplots for all the variables
AS_dotplot(Statfile,asterisk = "t_test") 

# Makes a subdirectory and saves barplots for all the variables
AS_barplot(Statfile,asterisk = "Scheffe") 
```


&emsp;&emsp;&emsp;&emsp;&emsp;**AS_boxplot(Statfile)**&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;**AS_dotplot(Statfile)**&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;**AS_barplot(Statfile)**
<p align="center">
<img src="https://user-images.githubusercontent.com/77651662/125714687-1908d6eb-b2bd-4e25-8ef0-62c24466c32a.png" width="225" height="250">
<img src="https://user-images.githubusercontent.com/77651662/125714704-d7dab67e-03c0-4e35-b86a-36723f7c63de.png" width="225" height="250">
<img src="https://user-images.githubusercontent.com/77651662/125715925-0878ec77-30bf-4859-8e56-316d98b6d520.jpg" width="350" height="250">
</p>

#### Adjustable parameters

* asterisk = "t_test"   #c("Dunn","Scheffe","u_test","t_test")
* significant_variable_only = F  # If set to TRUE, insignificant results will not be plotted
* color = c("#FF3300", "#FF6600", "#FFCC00", "#99CC00", "#0066CC", "#660099") # Colors for the plots
* legend_position = "none" #  "none","left","right","bottom","top"
* order = NULL # Order of the groups c("LAC","LUE","WEI","SDF","HGH","ASH")
* tip_length = 0.01 # significance tip length
* label_size = 2.88 # significance label size
* step_increase = 0.05 #significance step increase
* width = 0.3 # box width ; size = 3 # dot size
* fig_width = NA #figure size 
* fig_height = NA #figure size
### Multivariate statistics
#### PERMANOVA
```
data("Data")

data("Classification")
```

#### Single factor
PERMANOVA done with the Group column
```
Indiv_Perm(Data) # The group information is treated as a factor
```

<p align="center">
<img src="https://user-images.githubusercontent.com/77651662/125876315-c51448c9-eef1-4299-b85e-34e62a99bea8.PNG" width="750" height="225">
</p>

#### Multiple Factors

Loops PERMANOVA over different classes provided by Classification

```
Result<-Multi_Perm(Data,Classification) # The group information is treated as factors
```
<p align="center">
<img src="https://user-images.githubusercontent.com/77651662/125875797-89b696dd-f2a7-4ff3-aec9-74ac4b9075f9.PNG" width="750" height="400">
</p>



#### Adjustable parameters
* method = Dissimilarity index c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq",chord")

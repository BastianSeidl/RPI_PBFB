##### READ ME #####

#Output: 
#Shapiro_Result: p-values to check, whether the individual samples,are non-parametric distributed
#Boxplot for Length, Width and Geometric Mean per sample to show the overlap between samples
#Kruskal and Dune Results for samples of original size, to show their significant differences
#Sub_PBFB_result: p-values to check, whether test results are an artifact of sample size 

##### Loading of Data and required packages ######

#Loading the data into R with working directory
setwd()
PBFB <- read.csv("data/PBFB.csv")

#Loading of package for later use of Dunn Test
library(FSA)

#####Calculate geometric mean for body size#####

#Creating empty object to store results of Geometric Mean
Geometric_Mean <- NULL

#For-loop for calculating Geometric Mean and saving the result
for(i in length(PBFB)){
  result <-  sqrt(PBFB$Length*PBFB$Width)
  Geometric_Mean <- c(Geometric_Mean, result)
}

#Combing the existing PBFB with the Geometric Mean
PBFB <- cbind(PBFB, Geometric_Mean)

##### Checking whether Length and Width are normally distributed #####

#Creating vector with unique names of Samples
list_sample <- unique(PBFB$Sample)

#Creating empty data frame to store results of Shapiro Test
Shapiro_Result <- data.frame(`Sample` = character(),
                             `Length` = numeric(),
                             `Width` = numeric(),
                             `Geometric Mean` =  numeric())

#For-loop for calculating Geometric Mean and saving the result
for (i in 1:length(list_sample)) {
  cur_sample <- list_sample[i]
  cur_dat <- PBFB[cur_sample == PBFB$Sample,]
  Shapiro_Length <- shapiro.test(cur_dat$Length)
  Shapiro_Width <- shapiro.test(cur_dat$Width)
  Shapiro_Mean <- shapiro.test(cur_dat$Geometric_Mean)
  Shapiro_Result[i, "Sample"] <- list_sample[i] 
  Shapiro_Result[i, "Length"] <- Shapiro_Length["p.value"]
  Shapiro_Result[i, "Width"] <- Shapiro_Width["p.value"]
  Shapiro_Result[i, "Geometric.Mean"] <- Shapiro_Mean["p.value"]
}

##### Creating boxplots for comparing the sample Length, Width and Geometric Mean #####

windows(width=5, height=10)
par(xaxs="i", yaxs="i", mfrow=c(3,1))

#Creating the boxplot
boxplot(Length ~ Sample, data = PBFB, main = "Length", xlab = "", ylab = "Length in mm")
mtext("A", side = 3, line = -2.3, at = 0.7, cex = 2)
boxplot(Width ~ Sample, data = PBFB, main = "Width", xlab = "", ylab = "Width in mm")
mtext("B", side = 3, line = -2.3, at = 0.7, cex = 2)
boxplot(Geometric_Mean ~ Sample, data = PBFB, main = "Geometric Mean", xlab = "Sample Group", ylab = "Geometric Mean")
mtext("C", side = 3, line = -2.3, at = 0.7, cex = 2)

#Adding the number of measured specimen per sample 
for (i in 1:length(list_sample)) {
  cur_sample <- list_sample[i]
  cur_dat <- PBFB[cur_sample == PBFB$Sample,]
  mtext(paste0("n = ", nrow(cur_dat)), side = 1, line = 2.1, at = i, cex = 0.8)
}

par(xaxs="i", yaxs="i", mfrow = c(1, 1)) 

##### Test for significant differences between samples #####

kruskal_length <- kruskal.test(Length ~ Sample, data = PBFB)
kruskal_width <- kruskal.test(Width ~ Sample, data = PBFB)
kruskal_mean <- kruskal.test(Geometric_Mean ~ Sample, data = PBFB)

##### Dunn Test for comparing the sample one by one #####

PBFB$Sample <- factor(PBFB$Sample)
library(FSA)

dunn_length <- dunnTest(Length ~ Sample, data=PBFB, method="bonferroni")
dunn_width <- dunnTest(Width ~ Sample, data=PBFB, method="bonferroni")
dunn_mean <- dunnTest(Geometric_Mean ~ Sample, data=PBFB, method="bonferroni")

##### Sample standardization using 1000 subsamples of 60 measurements per sample #####
#Setting a seed to make it reproducable
set.seed(735)

#setting the number of subsamples
num_subPBFB <- 1000

# Create an empty dataframe to store the results
sub_PBFB_result <- data.frame(Repeat = numeric(),
                                   `Kuskal_p-value` = numeric(),
                                   `A001-A018` = numeric(),
                                   `A001-A033` = numeric(),
                                   `A018-A033` = numeric(),
                                   `A001-A045` = numeric(),
                                   `A018-A045` = numeric(),
                                   `A033-A045` = numeric())


#For-loop for sampling 60 measurements per sample and doing 1000 runs
for(i in 1:num_subPBFB){
  
  # Create an empty dataframe to store the temporary results
  sub_PBFB_temp <- data.frame(Sample = character(), Geometric_Mean = numeric())
  
  sub_PBFB_temp_result <- data.frame(Repeat = numeric(),
                                     `Kuskal_p-value` = numeric(),
                                     `A001-A018` = numeric(),
                                     `A001-A033` = numeric(),
                                     `A018-A033` = numeric(),
                                     `A001-A045` = numeric(),
                                     `A018-A045` = numeric(),
                                     `A033-A045` = numeric())
  
  
  for (j in 1:length(list_sample)) {
    cur_sample <- list_sample[j]
    cur_dat <- PBFB[cur_sample == PBFB$Sample, ]
    result_df <- sample(cur_dat$Geometric_Mean, 60, replace = FALSE)
    # Create a temporary dataframe with the current sample and corresponding geometric means
    temp_df <- data.frame(Sample = rep(cur_sample, 60), Geometric_Mean = result_df)
    # Append the temporary dataframe to the result dataframe
    sub_PBFB_temp <- rbind(sub_PBFB_temp, temp_df)
  }

  kruskal_mean <- kruskal.test(Geometric_Mean ~ Sample, data = sub_PBFB_temp)
  dunn_mean <- dunnTest(Geometric_Mean ~ Sample, data=sub_PBFB_temp, method="bonferroni")

  dunn_adj.p <- dunn_mean$res$P.adj
  
  sub_PBFB_temp_result <- data.frame(Repeat = i,
                                     `Kuskal_p-value` = kruskal_mean$p.value,
                                     `A001-A018` = dunn_adj.p[1],
                                     `A001-A033` = dunn_adj.p[2],
                                     `A018-A033` = dunn_adj.p[3],
                                     `A001-A045` = dunn_adj.p[4],
                                     `A018-A045` = dunn_adj.p[5],
                                     `A033-A045` = dunn_adj.p[6])
  
  sub_PBFB_result <- rbind(sub_PBFB_result, sub_PBFB_temp_result)

}


#Check

larger.05 <- sub_PBFB_result[sub_PBFB_result$Kuskal_p.value > 0.05, ]

test1.18 <- sub_PBFB_result[sub_PBFB_result$A001.A018 > 0.05, ]

test1.33 <- sub_PBFB_result[sub_PBFB_result$A001.A033 > 0.05, ]

test18.33 <- sub_PBFB_result[sub_PBFB_result$A018.A033 > 0.05, ]

test1.45 <- sub_PBFB_result[sub_PBFB_result$A001.A045 > 0.05, ]

test18.45 <- sub_PBFB_result[sub_PBFB_result$A018.A045 > 0.05, ]

test33.45 <- sub_PBFB_result[sub_PBFB_result$A033.A045 > 0.05, ]

#Combine results for use in boxplot 
result_counts <- c(abs(nrow(larger.05)-num_subPBFB), abs(nrow(test1.18)-num_subPBFB), abs(nrow(test1.33)-num_subPBFB), abs(nrow(test18.33)-num_subPBFB),abs(nrow(test1.45)-num_subPBFB), abs(nrow(test18.45)-num_subPBFB), abs(nrow(test33.45)-num_subPBFB))

# Create a bar plot for visualizing results of sample standartization
barplot(result_counts, names.arg = c("Kruskal", "1-18", "1-33", "18-33", "1-45", "18-45", "33-45"), 
        xlab = "Entry", ylab = "Count", main = "p values <0.05 (significant differnce in...)")


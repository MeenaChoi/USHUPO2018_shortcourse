##############################
##############################
## US HUPO 2018 Short course - Section 5 : Introduction to R and basic statistics
## Date: March 10, 2018
## Created by Meena Choi
##############################
##############################

# Summary

# Creating a new RStudio project
# Reading in data in R
# Data exploration, subsetting and replacement
# Visualizing data
# Select random sample and randomize MS run orders.
# Calculate simple statistics and visualize them using ggplot2.
# Statistical hypothesis testing by t-test.
# Saving your work

#############################################
## 1. Create a new Rstudio project
#############################################

# From the menu, select **File > New Project...**, 
# then select **Existing Directory** and choose the directory where you downloaded this script and the example datasets for this tutorial. 
# All the output files we'll be creating in this tutorial will be saved in the 'working directory' that now has been set by Rstudio.
  
# Let's verify the working directory path with the get working directory command.
getwd()


#############################################
## 2. Reading in data
#############################################

iprg <- read.csv("iPRG_example_runsummary.csv")


#############################################
## 3. Data exploration
#############################################

#`class` shows the type of a variable, in this case a 'data.frame`.
class(iprg)

# `dim` shows the dimension of a data.frame, which are the number of rows and the number of columns
dim(iprg)

#`colnames` is short for column names. 
colnames(iprg)
rownames(irpg)
#`head` shows the first 6 rows of data. Try `tail` to show the last 6 rows of data.
head(iprg)
tail(iprg)
# Let's explore the type of every column/variable and a summary for the value range for every column.
summary(iprg)

iprg[100, 2]
# Inspect the possible values for the `Conditions` and the `BioReplicate` (8th) column 
# using the named and numbered column selection syntax for data frames.
unique(iprg[, 'Condition'])
unique(iprg[, 4])

unique(iprg[, c('Condition', 'BioReplicate', 'Run')])
# Select subsets of rows from iprg dataset: 
# i.e we might be interested in working from here on only with the Condition1
# or all measurements on one particular MS run.

iprg.condition1 <- iprg[iprg$Condition == 'Condition1', ]
iprg.condition1.bio1 <- iprg[iprg$Condition == 'Condition1' & iprg$BioReplicate == '1', ]
nrow(iprg.condition1.bio1)

#############################################
## 4. Summarizing and Visualizing data
#############################################

### 4.1 Histogram
# Make a histogram of all the MS1 intensities, quantified by Skyline, for `iPRG_example`.

hist(iprg$Intensity)

# Our histogram looks quite skewed. How does this look on log-scale? 
# Do you recognize this distribution? The distribution for log2-transformed intensities looks very similar to the normal distribution. 

hist(iprg$Log2Intensity, 
     main='Histogram of iPRG data',
     xlab='log2 transformed intensities')

hist(log2(iprg$Intensity), 
     xlab="log2 transformed intensities", main="Histogram of iPRG data")

# We look at the summary for the log2-transformed values including the value for the mean.
summary(iprg$Log2Intensity)


### 4.2 Boxplot or box-and-whisker plot

# Boxplots are extremely useful becasue they allow us to quickly visualize the data distribution, 
# without making assumptions of the distribution type (non-parametric). 

# Let's make the boxplot with `ggplot2`, one of the most popular and powerful R packages 
# for making graphics. 

# Load ggplot2
library(ggplot2)

ggplot(aes_string(x='Run', y='Log2Intensity'), data=iprg)+
  geom_boxplot(aes_string(fill='Condition'))

# Let's rename all axis labels and title, and rotate the x-axis labels 90 degrees. 
# We can add those specifications using the `labs` and `theme` functions of the `ggplot2` package.

ggplot(aes_string(x='Run', y='Log2Intensity'), data=iprg)+
  geom_boxplot(aes_string(fill='Condition'))+
  labs(title='Log2 transformed intensity distribution per MS run', 
       y='Log2(Intensity)',
       x='MS run')+
  theme(axis.text.x=element_text(angle=90))

# And easily switch from a boxplot to a violin plot representation by changing the `geom` type. 

ggplot(aes_string(x='Run', y='Log2Intensity'), data=iprg)+
  geom_violin(aes_string(fill='Condition'))+
  labs(title='Log2 transformed intensity distribution per Subject', 
       y='Log2(Intensity)',
       x='MS run')+
  theme(axis.text.x=element_text(angle=90))


#############################################
# 5. Randomization
#############################################

## 5.1 Random selection of samples from a larger set

# This particular dataset contains a total of 10 subjects across conditions. 
# Suppose we label them from 1 to 10 
# and randomly would like to select 3 subjects we can do this using the `sample` function. 
# When we run `sample` another time, different subjects will be selected. Try this a couple times.

sample(1000, 10)
sample(10, 3)

# Now suppose we would like to select the same randomly selected samples every time, 
# then we can use a random seed number.

set.seed(3728)
sample(10, 3)

set.seed(3728)
sample(10, 3)


## 5.2 Completely randomized order of MS runs 

# We can also create a random order using all elements of iPRG dataset. 
# Again, we can achieve this using `sample`, asking for exactly the amount of samples in the subset. 
# This time, each repetition gives us a different order of the complete set.

msrun <- unique(iprg$Run)
msrun

# randomize order among all 12 MS runs
sample(msrun, length(msrun))

# different order will be shown.
sample(msrun, length(msrun))


## 5.3 Randomized block design

## Allow to remove known sources of variability that you are not interested in.
## Group conditions into blocks such that the conditions in a block are as similar as possible.
## Randomly assign samples with a block.

# This particular dataset contains a total of 12 MS runs across 4 conditions, 
# 3 technical replicates per condition. Using the `block.random` function in the `psych` package, 
# we can achieve randomized block designs!

# use 'psych' package. Load the package first.
library(psych)

msrun <- unique(iprg[, c('Condition','Run')])
msrun

# 4 Conditions of 12 MS runs randomly ordered
block.random(n=12, c(Group=4))


#############################################
# 6. Basic statistical summaries in R
#############################################

## 6.1 Calculate simple statistics

# Let's start with one protein as an example 
# and calculate the mean, standard deviation, standard error of the mean across all replicates per condition. 
# We then store all the computed statistics into a single summary data frame for easy access.

# We can use the **aggregate** function to compute summary statistics

# check what proteins are in dataset, show all protein names
unique(iprg$Protein)

# Let's start with one protein, named "sp|P44015|VAC2_YEAST"
oneproteindata <- iprg[iprg$Protein == "sp|P44015|VAC2_YEAST", ]

# If you want to see more details, 
?aggregate

### 6.1.1 Calculate mean per groups
# splits 'oneproteindata' into subsets by 'Condition', 
# then, compute 'FUN=mean' of 'log2Int'
sub.mean <- aggregate(Log2Intensity ~ Condition, data=oneproteindata, FUN=mean)
sub.mean

### 6.1.2 Calculate SD(standard deviation) per groups
# The same as mean calculation above. 'FUN' is changed to 'sd'.
sub.sd <- aggregate(Log2Intensity ~ Condition, data=oneproteindata, FUN=sd)
sub.sd

### 6.1.3 Count the number of observation per groups
# The same as mean calculation. 'FUN' is changed 'length'.
sub.len <- aggregate(Log2Intensity ~ Condition, data=oneproteindata, FUN=length)
sub.len

### 6.1.4 Calculate SE(standard error of mean) per groups
# SE = sqrt(s^2 / n)
sub.se <- sqrt(sub.sd$Log2Intensity^2/sub.len$Log2Intensity)
sub.se

# make the summary table including the results above (mean, sd, se and length).
summaryresult <- data.frame(Group=c("Condition1", "Condition2", "Condition3", "Condition4"),
                            mean=sub.mean$Log2Intensity,
                            sd=sub.sd$Log2Intensity, 
                            se=sub.se, 
                            length=sub.len$Log2Intensity)
summaryresult

## 6.2 Visualization with error bars for descriptive purpose
#‘error bars’ can have a variety of meanings or conclusions if what they represent is not precisely specified. 
# Below we provide some examples of which types of error bars are common. 
# We're using the summary of protein `sp|P44015|VAC2_YEAST` from the previous section and the `ggplot2` package as it provides a convenient way to make easily adaptable plots.

# Let's draw plots with mean and error bars

# means without any errorbar
ggplot(aes(x=Group, y=mean, colour=Group), data=summaryresult)+
  geom_point(size=3)

# Let's change a number of visual properties to make the plot more atttractive
# Let's change the labels of x-axis and y-axis and title: 
# add labs(title="Mean", x="Condition", y='Log2(Intensity)')
# Let's change background color for white : add theme_bw()
# Let's change size or color of labels of axes and title, text of x-axis : in theme
# Let's change the position of legend :'none' remove the legend
# Let's make the box for legend
# Let's remove the box for legend key.

ggplot(aes(x=Group, y=mean, colour=Group), data=summaryresult)+
  geom_point(size=3)+
  labs(title="Mean", x="Group", y='Log2(Intensity)')+
  theme_bw()+
  theme(plot.title = element_text(size=25, colour="darkblue"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=13),
        legend.position = 'bottom',
        legend.background = element_rect(colour='black'),
        legend.key = element_rect(colour='white'))

# Very similar but now as a bar plot.
ggplot(aes(x=Group, y=mean, fill=Group), data=summaryresult)+
  geom_bar(position=position_dodge(), stat='identity')+
  scale_x_discrete('Group')+
  labs(title="Mean", x="Group", y='Log2(Intensity)')+
  theme_bw()+
  theme(plot.title = element_text(size=25, colour="darkblue"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=13),
        legend.background = element_rect(colour='black'),
        legend.key = element_rect(colour='white'))

# For the sake of this tutorial, we'll continue adding error bars for different statistics with the point plots. 
# We'll leave it as an exercise to add error bars to the barplots. 
# Let's first add the standard deviation, then the standard error of the mean. Which one is smaller?

# mean with SD
ggplot(aes(x=Group, y=mean, colour=Group), data=summaryresult)+
  geom_point(size=3)+
  geom_errorbar(aes(ymax = mean + sd, ymin=mean - sd), width=0.5)+
  scale_x_discrete('Group')+
  labs(title="Mean with SD", x="Group", y='Log2(Intensity)')+
  theme_bw()+
  theme(plot.title = element_text(size=25, colour="darkblue"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=13),
        legend.position = 'bottom',
        legend.background = element_rect(colour='black'),
        legend.key = element_rect(colour='white'))

# mean with SE
ggplot(aes(x=Group, y=mean, colour=Group), data=summaryresult)+
  geom_point(size=3)+
  geom_errorbar(aes(ymax = mean + se, ymin=mean - se), width=0.1)+
  labs(title="Mean with SE", x="Condition", y='Log2(Intensity)')+
  theme_bw()+
  theme(plot.title = element_text(size=25, colour="darkblue"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=13),
        legend.position = 'bottom',
        legend.background = element_rect(colour='black'),
        legend.key = element_rect(colour='white'))

# The SE is narrow than the SD!

## 6.3 Calculate the confidence interval

# Now that we've covered the standard error of the mean and the standard deviation, 
# let's investigate how we can add custom confidence intervals (CI) for our measurement of the mean. 
# We'll add these CI's to the summary results we previously stored for protein `sp|P44015|VAC2_YEAST`

# Confidence interval : mean + or - (SE * alpha /2 { quantile of t distribution})$

# 95% confident interval
# Be careful for setting quantile for two-sided. need to divide by two for error.
# For example, 95% confidence interval, right tail is 2.5% and left tail is 2.5%.

summaryresult$ciw.lower.95 <- summaryresult$mean - qt(0.975,summaryresult$len-1)*summaryresult$se
summaryresult$ciw.upper.95 <- summaryresult$mean + qt(0.975,summaryresult$len-1)*summaryresult$se
summaryresult

# mean with 95% two-sided confidence interval
ggplot(aes(x=Group, y=mean, colour=Group), data=summaryresult)+
  geom_point(size=3)+
  geom_errorbar(aes(ymax = ciw.upper.95, ymin=ciw.lower.95), width=0.1)+
  labs(title="Mean with 95% confidence interval", x="Condition", y='Log2(Intensity)')+
  theme_bw()+
  theme(plot.title = element_text(size=25, colour="darkblue"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=13),
        legend.position = 'bottom',
        legend.background = element_rect(colour='black'),
        legend.key = element_rect(colour='white'))

# Let's repeat that one more time for the 99% two-sided confidence interval. 

# mean with 99% two-sided confidence interval
summaryresult$ciw.lower.99 <- summaryresult$mean - qt(0.995,summaryresult$len-1)*summaryresult$se
summaryresult$ciw.upper.99 <- summaryresult$mean + qt(0.995,summaryresult$len-1)*summaryresult$se
summaryresult

ggplot(aes(x=Group, y=mean, colour=Group), data=summaryresult)+
  geom_point(size=3)+
  geom_errorbar(aes(ymax = ciw.upper.99, ymin=ciw.lower.99), width=0.1)+
  labs(title="Mean with 99% confidence interval", x="Condition", y='Log2(Intensity)')+
  theme_bw()+
  theme(plot.title = element_text(size=25, colour="darkblue"),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=13),
        legend.position = 'bottom',
        legend.background = element_rect(colour='black'),
        legend.key = element_rect(colour='white'))


#############################################
# 7. Statistical hypothesis test in R
## Two sample t-test for one protein with one feature
#############################################

# Next, we'll perform a t-test whether protein `sp|P44015|VAC2_YEAST` has a change in abundance 
# between Condition 1 and Condition 2.

# If you want to see more details, 
?t.test

# First, get two conditions only, because t.test only works for two groups (conditions).
oneproteindata.condition12 <- oneproteindata[which(oneproteindata$Condition %in% 
                                                     c('Condition1', 'Condition2')), ]
unique(oneproteindata.condition12$Condition)
unique(oneproteindata$Condition)

# t test for different abundance (log2Int) between Groups (Condition)
result <- t.test(oneproteindata.condition12$Log2Intensity ~ oneproteindata.condition12$Condition,
                 var.equal=FALSE)
# show the summary of t-test including confidence interval level with 0.95
result

# We can redo the t-test and change the confidence interval level for the log2 fold change.

result.ci90 <- t.test(oneproteindata.condition12$Log2Intensity ~ oneproteindata.condition12$Condition, 
                      var.equal=FALSE, 
                      conf.level=0.9)
result.ci90

# Let's have a more detailed look at what information we can learn from the results our t-test. 

# name of output
names(result)

# mean for each group
result$estimate 

# log2 transformed fold change between groups : Disease-Healthy
result$estimate[1]-result$estimate[2]

# test statistic value, T value
result$statistic 

# standard error
(result$estimate[1]-result$estimate[2])/result$statistic

# degree of freedom
result$parameter 

# p value for two-sides testing
result$p.value 

# 95% confidence interval for log2 fold change
result$conf.int 

# p value calculation for one side
1-pt(result$statistic, result$parameter)

# p value for two sides, which is the same as pvalue from t test (result$p.value)
2*(1-pt(result$statistic, result$parameter))

# We can also manually compute our t-test statistic using the formulas we descibed above and 
# compare it with the `summaryresult` 

summaryresult

summaryresult12 <- summaryresult[c(1,4), ]

# test statistic, It is the same as 'result$statistic' above.
diff(summaryresult12$mean) # same as result$estimate[1]-result$estimate[2]
sqrt(sum(summaryresult12$sd^2/summaryresult12$length)) # same as stand error

diff(summaryresult12$mean)/sqrt(sum(summaryresult12$sd^2/summaryresult12$length))


#############################################
# 8. Saving your work 
#############################################

# You can save plots to a number of different file formats. 
# PDF is by far the most common format because it's lightweight, cross-platform 
# and scales up well but jpegs, pngs and a number of other file formats are also supported. 
# Let's redo the last barplot but save it to the file system this time. 

# Let's save the boxplot as pdf file.
pdf('boxplot_log2intensity_distribution_byMSrun.pdf', width=10, height=8)

ggplot(aes_string(x='Run', y='Log2Intensity'), data=iprg)+
  geom_boxplot(aes_string(fill='Condition'))+
  labs(title='Log2 transformed intensity distribution per MS run', 
       y='Log2(Intensity)',
       x='MS run')+
  theme(axis.text.x=element_text(angle=90))

dev.off() 

# Finally, we can save this whole session you worked so hard on! 
  
save.image(file='Day1.RData')

# let's give it a rest for today. Saving an .RData is the easiest way to pick up your work right where you left it!



rm(list=ls())
load(file = 'Day1.RData')




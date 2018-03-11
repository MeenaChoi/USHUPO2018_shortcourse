##############################
##############################
## US HUPO 2018 Short course - Section 7 : Protein siginificance analysis with MSstats
## Date: March 11, 2018
## Created by Meena Choi
##############################
##############################

##############################
## 0. Load MSstats
##############################

library(MSstats)
?MSstats

##############################
## 1. Read data
##############################
# set working directory

# read skyline output
raw <- read.csv(file="iPRG_10ppm_2rt_15cut_nosingle.csv")

# read annotation table
annot <- read.csv("iPRG_skyline_annotation.csv", header=TRUE)
annot

# reformating and pre-processing for Skyline output.
quant <- SkylinetoMSstatsFormat(raw, annotation=annot)
head(quant)

##############################
## 2. Processing data
##     Transformation = log2
##     Normalization = equalize median
##     Model-based run-level summarization (TMP) after imputation
##############################

quant.processed <- dataProcess(raw = quant, 
                               logTrans=2, 
                               normalization = 'equalizeMedians',
                               summaryMethod = 'TMP', 
                               MBimpute=TRUE,
                               censoredInt='0',
                               cutoffCensored='minFeature',
                               maxQuantileforCensored = 0.999)

# show the name of outputs
names(quant.processed)

# show reformated and normalized data.
# 'ABUNDANCE' column has normalized log2 transformed intensities.
head(quant.processed$ProcessedData)

# This table includes run-level summarized log2 intensities. (column : LogIntensities)
# Now one summarized log2 intensities per Protein and Run.
# NumMeasuredFeature : show how many features are used for run-level summarization.
#         If there is no missing value, it should be the number of features in certain protein.
# MissingPercentage : the number of missing features / the number of features in certain protein.
head(quant.processed$RunlevelData)

# show which summarization method is used.
head(quant.processed$SummaryMethod)

##############################
## 3. Data visualization
##############################

dataProcessPlots(data = quant.processed, 
                 type="QCplot", 
                 width=7, height=7,
                 which.Protein = 'allonly',
                 address='iPRG_skyline_equalizeNorm_')

# It will generate profile plots per protein. It will take a while
# Please run at home. It takes a while.
dataProcessPlots(data = quant.processed, 
                 type="Profileplot", 
                 featureName="NA",
                 width=7, height=7,
                 summaryPlot = TRUE,
                 originalPlot = FALSE,
                 address="iPRG_skyline_equalizeNorm_")

# Instead, make profile plot for only some.
dataProcessPlots(data = quant.processed, 
                 type="Profileplot", 
                 featureName="NA",
                 width=7, height=7,
                 which.Protein = 'sp|P44015|VAC2_YEAST',
                 address="iPRG_skyline_equalizeNorm_P44015")


# sp|P55752|ISCB_YEAST
# sp|P44374|SFG2_YEAST
# sp|P44983|UTR6_YEAST
# sp|P44683|PGA4_YEAST
# sp|P55249|ZRT4_YEAST

# and first few proteins


# Please run at home. It takes a while.
dataProcessPlots(data = quant.processed, 
                 type="conditionplot", 
                 width=7, height=7,
                 address="iPRG_skyline_equalizeNorm_")

# Instead, make profile plot for only some.

# sp|P44015|VAC2_YEAST
# sp|P55752|ISCB_YEAST
# sp|P44374|SFG2_YEAST
# sp|P44983|UTR6_YEAST
# sp|P44683|PGA4_YEAST
# sp|P55249|ZRT4_YEAST

# and first few proteins



##############################
## 4. Model-based comparison and adjustment for multiple testing
##############################

unique(quant.processed$ProcessedData$GROUP_ORIGINAL)

comparison1<-matrix(c(-1,1,0,0),nrow=1)
comparison2<-matrix(c(-1,0,1,0),nrow=1)
comparison3<-matrix(c(-1,0,0,1),nrow=1)
comparison4<-matrix(c(0,-1,1,0),nrow=1)
comparison5<-matrix(c(0,-1,0,1),nrow=1)
comparison6<-matrix(c(0,0,-1,1),nrow=1)
comparison<-rbind(comparison1, comparison2, comparison3, comparison4, comparison5, comparison6)
row.names(comparison)<-c("C2-C1","C3-C1","C4-C1","C3-C2","C4-C2","C4-C3")

test <- groupComparison(contrast.matrix=comparison, data=quant.processed)

names(test)

# Show test result
# Label : which comparison is reported.
# log2FC : estimated log2 fold change between conditions.
# adj.pvalue : adjusted p value by BH
# issue : detect whether this protein has any issue for comparison
#    such as, there is measurement in certain group, or no measurement at all.
# MissingPercentage : the number of missing intensities/total number of intensities 
#     in conditions your are interested in for comparison
# ImputationPercentage : the number of imputed intensities/total number of intensities 
#     in conditions your are interested in for comparison
head(test$ComparisonResult)

# After fitting linear model, residuals and fitted values can be shown.
head(test$ModelQC)

# Fitted model per protein
head(test$fittedmodel)

# save testing result as .csv file
Skyline.intensity.comparison.result <- test$ComparisonResult

write.csv(Skyline.intensity.comparison.result, 
          file='testResult.iprg.skyline.proteinlevel.csv')

head(Skyline.intensity.comparison.result)
SignificantProteins <- Skyline.intensity.comparison.result[Skyline.intensity.comparison.result$adj.pvalue < 0.05 ,]
nrow(SignificantProteins)

## practice : please find the result for spike-in protein

# sp|P44015|VAC2_YEAST
# sp|P55752|ISCB_YEAST
# sp|P44374|SFG2_YEAST
# sp|P44983|UTR6_YEAST
# sp|P44683|PGA4_YEAST
# sp|P55249|ZRT4_YEAST


##############################
## 5. Visualization for testing result
##############################

groupComparisonPlots(Skyline.intensity.comparison.result, 
                     type="VolcanoPlot", 
                     address="testResult_iprg_skyline_")

groupComparisonPlots(data = Skyline.intensity.comparison.result, 
                     type = 'VolcanoPlot',
                     sig = 0.05, 
                     FCcutoff = 2^2, 
                     address = 'testResult_iprg_skyline_FCcutoff4_')

groupComparisonPlots(Skyline.intensity.comparison.result, 
                     type="Heatmap", 
                     address="testResult_iprg_skyline_")

# Please run at home. It takes a while.
groupComparisonPlots(Skyline.intensity.comparison.result, 
                     type="ComparisonPlot", 
                     address="testResult_iprg_skyline_")


##############################
## 6. Verify the model assumption
##############################

# normal quantile-quantile plots
# Please run at home. It takes a while.

modelBasedQCPlots(data=test, type="QQPlots", 
                  width=5, height=5, 
                  address="iPRG_skyline_equalizeNorm_testResult_")

# residual plots
# Please run at home. It takes a while.

modelBasedQCPlots(data=test, type="ResidualPlots", 
                  width=5, height=5, 
                  address="iPRG_skyline_equalizeNorm_testResult_")



##############################
## 7. Power calculation
##############################
test.power <- designSampleSize(data = test$fittedmodel, 
                               desiredFC = c(1.1, 1.6), 
                               FDR = 0.05,
                               power = TRUE,
                               numSample = 3)
test.power

designSampleSizePlots(data = test.power)



##############################
## 8. Sample size calculation
##############################
samplesize <- designSampleSize(data = test$fittedmodel, 
                               desiredFC = c(1.1, 1.6), 
                               FDR = 0.05,
                               power = 0.9,
                               numSample = TRUE)
samplesize

designSampleSizePlots(data = samplesize)



##############################
## 9. sample quantification
##############################
sampleQuant <- quantification(quant.processed)
head(sampleQuant)


##############################
##############################
## Extra practice
##############################
##############################


##############################
## quantile normalization
##############################
quant.processed.quantile <- dataProcess(raw = quant, 
                               logTrans=2, 
                               normalization = 'quantile',
                               summaryMethod = 'TMP', 
                               MBimpute=TRUE,
                               censoredInt='0',
                               cutoffCensored='minFeature',
                               maxQuantileforCensored = 0.999)

dataProcessPlots(data = quant.processed.quantile, 
                 type="QCplot", 
                 width=7, height=7,
                 which.Protein = 1,
                 address='iPRG_skyline_quantile_')

dataProcessPlots(data = quant.processed.quantile, 
                 type="Profileplot", 
                 featureName="NA",
                 width=7, height=7,
                 which.Protein = 1,
                 address="iPRG_skyline_quantile_1_")


##############################
##############################
## Extra :  peptide-level analysis
##############################
##############################

head(quant)
quant.pep <- quant
head(quant.pep)

##############################
## 1. Replace 'ProteinName' with combination of 'ProteinName' and 'PeptideSequence'
##############################

quant.pep$ProteinName <- paste(quant.pep$ProteinName, 
                               quant.pep$PeptideSequence, 
                               sep = '_')

length(unique(quant.pep$ProteinName)) # 30500 peptides

save(quant.pep, file='quant.pep.RData')


##############################
## 2. Process the data
##############################

# Please run at home. It takes a while.
quant.pep.processed <- dataProcess(raw = quant.pep, 
                               logTrans=2, 
                               normalization = 'equalizeMedians',
                               summaryMethod = 'TMP', 
                               MBimpute=TRUE,
                               censoredInt='0',
                               cutoffCensored='minFeature',
                               maxQuantileforCensored = 0.999)

save(quant.pep.processed, file='quant.pep.processed.RData')

##############################
## 3. Data visualization
##############################

dataProcessPlots(data = quant.pep.processed, type="QCplot", 
                 width=7, height=7,
                 which.Protein = 1,
                 address='iPRG_skyline_equalizeNorm_Peptidelevel_')

# Please run at home. It takes a while.
dataProcessPlots(data = quant.pep.processed, type="Profileplot", 
                 featureName="NA",
                 width=7, height=7,
                 summaryPlot = FALSE,
                 address="iPRG_skyline_equalizeNorm_Peptidelevel_")

# Please run at home. It takes a while.
dataProcessPlots(data = quant.pep.processed, type="conditionplot", 
                 address="iPRG_skyline_equalizeNorm_Peptidelevel_")


##############################
## 4. Model-based comparison and adjustment for multiple testing
##############################

unique(quant.pep.processed$ProcessedData$GROUP_ORIGINAL)

comparison1<-matrix(c(-1,1,0,0),nrow=1)
comparison2<-matrix(c(-1,0,1,0),nrow=1)
comparison3<-matrix(c(-1,0,0,1),nrow=1)
comparison4<-matrix(c(0,-1,1,0),nrow=1)
comparison5<-matrix(c(0,-1,0,1),nrow=1)
comparison6<-matrix(c(0,0,-1,1),nrow=1)
comparison<-rbind(comparison1, comparison2, comparison3, comparison4, comparison5, comparison6)
row.names(comparison)<-c("C2-C1","C3-C1","C4-C1","C3-C2","C4-C2","C4-C3")

# Please run at home. It takes a while.
test.Peptidelevel <- groupComparison(contrast.matrix=comparison, 
                                     data=quant.pep.processed)

testResult.iprg.skyline.Peptidelevel <- test.Peptidelevel$ComparisonResult
save(testResult.iprg.skyline.Peptidelevel, file='testResult.iprg.skyline.Peptidelevel.RData')
write.csv(testResult.iprg.skyline.Peptidelevel , file='testResult.iprg.skyline.Peptidelevel.csv')

head(testResult.iprg.skyline.Peptidelevel)
SignificantProteins <- testResult.iprg.skyline.Peptidelevel[testResult.iprg.skyline.Peptidelevel$adj.pvalue < 0.05 , ]
nrow(testResult.iprg.skyline.Peptidelevel)

## let's check one peptide. Check profile plot (page 6693)
testResult.iprg.skyline.Peptidelevel[testResult.iprg.skyline.Peptidelevel$Protein=='sp|P22146|GAS1_YEAST_ALNDADIYVIADLAAPATSINR', ]

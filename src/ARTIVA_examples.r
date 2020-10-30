#### Set up some working directories ####
setwd("~/Dropbox/HILDA/Causal modelling")
library(ARTIVA)

#### Example of choosing priors from ARTIVA.pdf ####
data(priors)
kmax = 2
choosePriors(kmax,priors)

#### Example from Bayesian Networks in R p.80 ####
data(simulatedProfiles)

targets <- c("1", "10", "20", "TF3", "45", "50")
parents <- c("TF1", "TF2", "TF3", "TF4", "TF5")

results <- ARTIVAnet(
  targetData = simulatedProfiles[targets, ],
  parentData = simulatedProfiles[parents, ],
  targetNames = targets,
  parentNames = parents,
  niter = 50000,
  savePictures = FALSE)

head(results[, -7])
traceNetworks(results[1:10, ], edgesThreshold = .5)
#geneNetworkSummary(results, edgesThreshold = 0.5) # displays a summary table

#### ARTIVAnet example from ARTIVA.pdf p.7 ####
data(simulatedProfiles)

targets <- c("1", "10", "20", "TF3", "45", "50")
parents <- c("TF1", "TF2", "TF3", "TF4", "TF5")

results <- ARTIVAnet(
  targetData = simulatedProfiles[targets, ],
  parentData = simulatedProfiles[parents, ],
  targetNames = targets,
  parentNames = parents,
  niter = 20000,
  savePictures = TRUE)

print(results)

#### geneNetworkSummary example from ARTIVA.pdf p.23 ####
data(simulatedProfiles)

# Name of the target gene to be analyzed with ARTIVA
targetGene = 1
# Names of the parent genes (typically transcription factors)
parentGenes = c("TF1", "TF2", "TF3", "TF4", "TF5")
# Run the ARTIVAsubnet function
result = ARTIVAsubnet(targetData = simulatedProfiles[targetGene,],
                          parentData = simulatedProfiles[parentGenes,],
                          targetName = targetGene,
                          parentNames = parentGenes,
                          segMinLength = 2,
                          edgesThreshold = 0.6,
                          niter= 2000,
                          savePictures=TRUE)
# Print a summary of the obtained network
geneNetworkSummary(result$network, edgesThreshold = 0.6)
traceNetworks(results$network, edgesThreshold = 0.6)

#### Yeast example from ARTIVA.pdf p.33 ####
# Datasets related to the analysis of the genomic response of the yeast
# Saccharomyces cerevisiae to an environmental stress induced by
# benomyl (a toxic compound). Responses were measured in 5 different strains:
# 1 wildtype (WT) and 4 knockout (YAP1, PDR1, PDR3 and YRR1). For each group,
# the measured expression values for 5 time-points was taken. In this context, 
# regulatory associations between parent and target genes are proposed if the
# deletion of a parent gene significantly alters the expression measurements of
# the target genes.
# Analysis of the yeast data is presented in the original article of
# ARTIVA (Lebre et al. BMC Syst. Biol, 2010)
####
# Load the yeast dataset
data(yeast)
# This is a a list that comprises information for the 18 clusters of genes
# whose expression is identically modified in strains deleted for
# YAP1, PDR1, PDR3 and YRR1 transcription factors,
# compared to the wild type strain.
# As an illustration : analysis of one cluster
cluster=4
# Different genes in a cluster is considered as repeated measurements.
# Organisation of the different time point measurements is described in
# variable : yeast[[cluster]]$dataDescription
# Because of repeated measurements, the minimum segment length is set to
# segMinLength = 1.
# The parentdata is the experiment design (YAP1, PDR1, PDR3 and YRR1
# deletion) described in variable: yeast[[cluster]]$parentData
# Time delay between parent and target genes is fixed to dyn=0. (I think this
# is set to 0 because we are looking for associations revealed by differences
# between strains (as well as between time points))
## Not run:
ARTIVAtest = ARTIVAsubnet(targetData = yeast[[cluster]]$targetData,
                          targetName = yeast[[cluster]]$targetName,
                          parentData = yeast[[cluster]]$parentData,
                          parentNames = row.names(yeast[[cluster]]$parentData),
                          dataDescription = yeast[[cluster]]$dataDescription,
                          outputPath = paste("ARTIVA_Results_Cluster", cluster, sep = ""),
                          dyn = 0,
                          segMinLength = 1,
                          edgesThreshold = 0.7,
                          saveIterations = TRUE,
                          saveEstimations = TRUE,
                          savePictures = TRUE,
                          niter = 20000)

# Detailed results can be found in the folder named
# "ARTIVA_Results_Cluster4" (with the subfolders "Estimations" for
# detailed results of the estimated parameters and "Pictures" for
# graphical representations).

#### Example with ARTIVAnet ####
data(simulatedProfiles)

targets <- c("1")
parents <- c("TF1", "TF2", "1")

results <- ARTIVAnet(
  targetData = simulatedProfiles[targets, 1:14],
  parentData = simulatedProfiles[parents, 1:14],
  targetNames = targets,
  parentNames = parents,
  niter = 50000,
  saveEstimations = TRUE,
  saveIterations = TRUE,
  savePictures = TRUE)

print(results)

#### Trace plots and autoregressive plots ####
pathtoCoeff <- 'ARTIVAnet/IterationsCoeffSamples_target_1.txt'
pathtoVar <- 'ARTIVAnet/IterationsVarianceSamples_target_1.txt'

coeffsamples <- read.table(pathtoCoeff)
varsamples <- read.table(pathtoVar)
matplot(coeffsamples[20:nrow(coeffsamples), 1:4], type = "l", lty = 1)
matplot(varsamples[20:2000, 5], type = "l", lty = 1)

acf(coeffsamples$V1)

coeffsamples <- read.table('ARTIVA_Results_Cluster4/IterationsCoeffSamples_target_Cluster 4.txt')
varsample <- read.table('ARTIVA_Results_Cluster4/IterationsVarianceSamples_target_Cluster 4.txt')
matplot(coeffsamples[20:nrow(coeffsamples), 1:3], type = "l", lty = 1)
matplot(varsample[20:nrow(coeffsamples), 1:3], type = "l", lty = 1)
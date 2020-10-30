setwd("~/Dropbox/HILDA/Causal modelling")
library(EDISON)

#### Example from EDISON.pdf p.4 ####
# Generate random gene network and simulate data from it
dataset = simulateNetwork(l=25)
# Run MCMC simulation to infer networks and changepoint locations
result = EDISON.run(dataset$sim_data, num.iter=500)
# Calculate posterior probabilities of changepoints
cps = calculateCPProbabilities(result)
# Calculate marginal posterior probabilities of edges in the network
network = calculateEdgeProbabilities(result)

#### My example ####
dataset = simulateNetwork(l=14,         # length of time series
                          k_bar=1,      # number of changepoint
                          fixed = TRUE, # fix the location
                          cps = 7,      # location of changepoints
                          q = 3)        # number of parents
matplot(t(dataset$sim_data), type = 'l')

Options <- defaultOptions()             # Set default options
Options$lmax <- 3                       
Options$minPhase <- 5
Options$maxCP <- 1
Options$maxTF <- 3
Options$m <- 3

df <- dataset$sim_data
df3 <- rbind(df, df)
df3 <- rbind(df3, df)

result = EDISON.run(df3, 
                    num.iter = 500,
                    options = Options)
cps = calculateCPProbabilities(result)
matplot(t(cps$global.cps), type = 'l')
matplot(t(cps$node.cps), type = 'l')


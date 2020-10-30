# Install developer version
devtools::install_github(
  'FrankD/EDISON', ref='MultipleTimeSeries', subdir='Package/EDISON'
  )

# Test script for repeated sample time series inference
library(EDISON)

set.seed(10)

#### Brief test ####
# should take a few seconds to run
# 
# Generate a network with 5 nodes and no changepoints
test = generateNetwork(l=10, q=5, k_bar=0)

# Generate 100 repeated measurements for 10 time points
test.data = lapply(1:100, simulateNetwork, l=10, net=test)

# Make array and put dimensions in the right order (repeats, time points, variables)
test.data.array = sapply(test.data, function(x) x$sim_data, simplify='array')
test.data.array = aperm(test.data.array, c(3,2,1))

# Specify inference with no changepoints
edison.options = defaultOptions()
edison.options$cp.fixed = TRUE

# Run EDISON
edison.test = EDISON.run(test.data.array, num.iter=10000, options=edison.options)

# Calculate posterior probabilities of the edges in the network
calculateEdgeProbabilities(edison.test)$probs.segs

#### Longer test ####
#
# Generate a network with 9 nodes and no changepoints
test = generateNetwork(l=16, q=9, k_bar=0)

# Generate 1000 repeated measurements for 16 time points
test.data = lapply(1:1000, simulateNetwork, l=16, net=test)

# Make array and put dimensions in the right order (repeats, time points, variables)
test.data.array = sapply(test.data, function(x) x$sim_data, simplify='array')
test.data.array = aperm(test.data.array, c(3,2,1))

# Run and time EDISON (this originally took 1.2 days to run)
start.time <- Sys.time()
edison.test = EDISON.run(test.data.array, num.iter=10000, options=edison.options)
end.time <- Sys.time()
print(end.time - start.time)
# > Time difference of 1.444573 mins
# 
# Calculate posterior probabilities of the edges in the network
calculateEdgeProbabilities(edison.test)$probs.segs

# Test script for repeated sample time series inference
# 
# Install the developer version of EDISON (check forward slashes if error)
devtools::install_github(repo='FrankD/EDISON',ref='MultipleTimeSeries', 
                         subdir='/Package/EDISON')

# Setup
library(EDISON)
set.seed(10)

# Generate a timeseries network with 3 nodes, 6 timepoints and no changepoints
homogenous_DBN = generateNetwork(l=6, q=3, k_bar=0, fixed = TRUE)
homogenous_DBN$network

          [,1]      [,2]      [,3]
[1,] 0.0000000 0.4829785  0.000000
[2,] 0.0000000 0.0000000  0.000000
[3,] 0.9255213 0.0000000 -1.265198

# simulateNetwork(net = homogenous_DBN)

# Generate 1000 repeated measurements for 6 time points
test.data = lapply(1:1000, simulateNetwork, l=6, net=homogenous_DBN)

# Make array and put dimensions in the right order (repeats, time points, variables)
test.data.array = sapply(test.data, function(x) x$sim_data, simplify='array')
test.data.array = aperm(test.data.array, c(3,2,1))

# Specify inference with no changepoints
edison.options = defaultOptions()
edison.options$cp.fixed = FALSE   # default False
edison.options$m = 1000           # default 1
edison.options$minPhase = 1       # default 2
edison.options$maxCP = 5          # default 10
edison.options$lmax = 3           # default 5
edison.options$dyn = 1            # default 1 
edison.options$maxTF = 3



# Run EDISON (if returns error, check maxTF)
edison.test = EDISON.run(test.data.array, num.iter=5000, options=edison.options)
# Error in sample.int(length(x), size, replace, prob) : 
#   cannot take a sample larger than the population when 'replace = FALSE'

# Calculate posterior probabilities of the edges in the network at each time point
(posteriors <- calculateEdgeProbabilitiesTimePoints(edison.test, cps = 1:6, numNodes = 3))

# Marginal posteriors
(marginal_posteriors <- calculateEdgeProbabilities(edison.test)$probs.segs[[1]])

      [,1]  [,2]  [,3]
[1,] 0.000 1.000 0.014
[2,] 0.022 0.000 0.000
[3,] 1.000 0.097 1.000

# The marginal posteriors match the homogenous network above

#### Plotting ####
library(tidyverse)
vnames <- c("bmi", "sleep", "mental")

# Time series plot
EdgeList <- data.frame()

for (i in seq(length(posteriors))) {
  post_i <- posteriors[[i]]
  colnames(post_i) <- paste0(vnames, i+1)
  rownames(post_i) <- paste0(vnames, i)
  post_i %>%
    as.data.frame() %>%
    tibble::rownames_to_column("from") %>%
    tidyr::gather(to, weight, -from) %>%
    dplyr::bind_rows(EdgeList) -> EdgeList
}


nodes <- c(paste0("bmi", seq(1,6)), paste0("sleep", seq(1,6)), paste0("mental", seq(1,6)))
x <- rep(seq(0,5), 3)
y <- rep(seq(0,2), each = 6)
NodeList <- data.frame(nodes, x, y)

library(igraph)
var1_net <- graph_from_data_frame(vertices = NodeList, d= EdgeList, directed = FALSE)
plot(var1_net, vertex.size = 30, edge.width = E(var1_net)$weight*2)

# Marginal plot
# (to be done)
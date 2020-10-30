#### Set up ####
# setwd('~/Dropbox/HILDA/causal-models')
# datadir <- ("~/Dropbox/HILDA/data/")

library(tidyverse)
# library(haven)
library(EDISON)

#### and some helper functions ####
source("~/Dropbox/HILDA/src/JoinFrames.R")
source("~/Dropbox/HILDA/src/GetVarsLong.R")
source("~/Dropbox/HILDA/src/ZapLabel.R")

#### Load the HILDA dataset ####
filepaths <- list.files(
  path = datadir,
  pattern = '^Combined.*.dta$',
  full.names = TRUE
)

cat('Importing dataframes')
hilda <- list()
for (pathtofile in filepaths) {
  df <- read_dta(pathtofile)
  hilda <- append(hilda, list(df))
  cat('.')
}
cat('Done')

#### Create the variable matrices ####
predictors <- c('ghmh',       # mental health
                'losat',      # life satisfaction
                'losatyh',    # health satisfaction
                'wscei',      # weekly gross wages/salary
                'lssupac',    # I don't have anyone I can confide in
                'lssupvl',    # I often feel very lonely
                'lsrush',     # feeling rushed
                'lshrhw',     # hours per week housework
                'lsshare'     # share of housework
)

df.long <- GetVarsLong(hilda, predictors)

matrixlist <- list()
for (code in predictors) {
  df.long %>%
    filter(hildacode == code) %>%
    spread(waveid, value) %>%
    select(-hildacode, -xwaveid) %>%
    as.matrix() -> mat
  matrixlist[[code]] <- mat
}

hilda.array <- simplify2array(matrixlist)
write_rds(hilda.array, 'data/hildadata.rds')
hilda.array <- read_rds('data/hildadata.rds')
summary(hilda.array) # no NAs
hist(hilda.array)

# Specify inference with no changepoints
edison.options = defaultOptions()
edison.options$cp.fixed = TRUE

# Run EDISON
# randomidx <- sample.int(n = 2917, size = 1000) # random selection of row
# hilda.random <- hilda.array[randomidx, , ]

start.time <- Sys.time()
# hilda.test = EDISON.run(hilda.random, num.iter=10000, options=edison.options)
hilda.test = EDISON.run(hilda.array, num.iter=50000, options=edison.options)
end.time <- Sys.time()
print(end.time - start.time)
# Time difference of 1.248522 days (original EDISON)

# Calculate posterior probabilities of the edges in the network
ans <- calculateEdgeProbabilities(hilda.test)$probs.segs[[1]]
colnames(ans) = rownames(ans) = predictors
ans

# Full dataset, num.iter = 10000, developer version, time difference of 4.242534 mins
        ghmh losat losatyh wscei lssupac lssupvl lsrush lshrhw lsshare
ghmh       0     0       0     0       1       0      1      0       0
losat      1     0       1     0       1       0      1      1       0
losatyh    0     1       1     1       1       0      0      0       1
wscei      1     1       1     1       0       1      1      1       1
lssupac    0     1       0     0       1       1      0      0       0
lssupvl    1     1       1     0       0       1      0      0       1
lsrush     1     0       0     1       0       1      1      0       1
lshrhw     0     0       0     1       0       0      0      1       0
lsshare    0     0       0     0       0       0      0      1       0


network <- igraph::graph_from_adjacency_matrix(ans)
plot(network)

# Calculate posterior probabilities of the edges in the network
ans2 <- calculateEdgeProbabilities(hilda.test)$probs.segs[[1]]
colnames(ans2) = rownames(ans2) = predictors
ans2

# Full dataset, num.iter = 50000, developer version, time difference of 20.81051 mins
        ghmh losat losatyh wscei lssupac lssupvl lsrush lshrhw lsshare
ghmh       0     0       1     0       1       1      1      0       0
losat      1     0       1     0       1       1      0      1       0
losatyh    0     1       0     1       1       1      0      0       0
wscei      0     0       0     1       0       0      1      1       0
lssupac    1     1       0     0       0       0      0      0       0
lssupvl    0     1       0     0       1       0      0      0       1
lsrush     1     1       0     1       0       0      1      0       1
lshrhw     1     0       1     1       0       0      1      1       1
lsshare    0     0       1     0       0       1      0      1       1

network2 <- igraph::graph_from_adjacency_matrix(ans2)
plot(network2)

# Calculate posterior probabilities of the edges in the network
ans3 <- calculateEdgeProbabilities(hilda.test)$probs.segs[[1]]
colnames(ans3) = rownames(ans3) = predictors
ans3

# Full dataset, num.iter = 50000, developer version, time difference of 20.26309 mins
        ghmh losat losatyh wscei lssupac lssupvl lsrush lshrhw lsshare
ghmh       1     0       1     0       0       1      0   0.00       0
losat      1     0       1     0       1       0      0   1.00       0
losatyh    1     1       1     1       1       0      1   0.63       0
wscei      1     0       1     1       0       0      0   0.37       1
lssupac    0     0       0     0       1       1      0   0.00       0
lssupvl    0     1       0     0       1       1      1   0.00       0
lsrush     1     0       0     1       0       0      0   0.00       1
lshrhw     0     1       0     0       0       1      1   1.00       1
lsshare    0     1       0     1       0       0      1   1.00       1

# The results with the full dataset are not stable. I wonder if stationarity is
# causing this? We could try mean differencing for each individual. Before then,
# let's test a smaller set of variables to increase stability...

hilda.array.subset <- hilda.array[, , c('ghmh', 'losat', 'losatyh', 'lssupvl')]
start.time <- Sys.time()
hilda.test = EDISON.run(hilda.array.subset, num.iter=50000, options=edison.options)
end.time <- Sys.time()
print(end.time - start.time)

# Calculate posterior probabilities of the edges in the network
ans <- calculateEdgeProbabilities(hilda.test)$probs.segs[[1]]
colnames(ans) = rownames(ans) = c('ghmh', 'losat', 'losatyh', 'lssupvl')
ans

# Subset dataset, num.iter = 50000, developer version, time difference of 15.61926 mins
        ghmh losat losatyh lssupvl
ghmh       1     1       1       1
losat      1     1       1       1
losatyh    1     1       1       1
lssupvl    1     1       1       1

library(tidyverse)

# Import the different parametric response-groups identified in for-better-or-
# worse
#
# func to load the subgroups classified by the largest response component
find_subgroups <- function(parametric_fit) {
  
  coefs <- coef(parametric_fit)[[1]] # a 3D array 
  individualCoefs <- data.frame(
    "xwaveid" = rownames(coefs[, , 1]),
    "comp1" = unname(coefs[, 1, 2]),
    "comp2" = unname(coefs[, 1, 3]),
    "comp3" = unname(coefs[, 1, 4])
  )
  
  # Note: ranef() and fixef() can be used to extract the random and fixed effects
  # 
  # Find the maximum component and its sign
  individualCoefs$maxc <- max.col(abs(individualCoefs[, -1]))
  idx = individualCoefs$maxc+1
  individualCoefs$maxc <- individualCoefs$maxc*
    sign(as.numeric(individualCoefs[cbind(seq(length(idx)), idx)]))
  
  return(individualCoefs)
}

# event names
event_key <- c(
  `ledsc` = "widowed",
  `lefnw` = "bankruptcy",
  `lesep` = "divorce"
)

# import subgroups
subgroups <- map_dfr(
  .x = read_rds("data/parametric_fits.rds"), 
  .f = ~find_subgroups(.x), 
  .id = "eventcode"
) %>%
  mutate(eventcode = recode(eventcode, !!!event_key))

subgroups <- filter(subgroups, maxc %in% c(1, -1, 2)) %>%
  mutate(response = recode_factor(maxc, 
                                  `2` = "resilient",
                                  `-1` = "vulnerable", 
                                  `1` = "relieved", 
                                  .default = "undefined"))

# Import outcomes from for-better-or-worse
outcomes <- read_rds("data/outcomes.rds") %>%
  rowwise() %>%
  mutate(year = which(letters == wave) + 2000) %>%
  ungroup()

# Import covariates from for-better-or-worse and bind with outcomes
covariates <- read_rds("data/covariates.rds") %>%
  select(-female) %>%
  bind_rows(outcomes)

# > head(covariates)
# # A tibble: 6 x 5
#   xwaveid code     val wave   year
#   <chr>   <chr>  <dbl> <chr> <dbl>
# 1 0100001 hhda10     3 a      2001
# 2 0100002 hhda10     3 a      2001
# 3 0100003 hhda10     1 a      2001
# 4 0100004 hhda10     1 a      2001
# 5 0100005 hhda10     1 a      2001
# 6 0100006 hhda10     1 a      2001

> unique(covariates$code)
[1] "hhda10"  "hgage"   "esbrd"   "mrcurr"  "edfts"   "edhigh1" "losatlc"
[8] "losatyh" "lspact"  "lsrush"  "lsstime" "lsclub"  "lssocal" "lssuppv"
[15] "lssupnh" "lssuplf" "lssupac" "lssuplt" "lssupcd" "lssupvl" "lssuppi"
[22] "lssuptp" "lssupsh" "lsvol"   "lssefh"  "lssepa"  "ctbds"   "ctwps"  
[29] "ctsds"   "losat"   "gh9_sum"

# Import event details
read_rds("data/event_years.rds") %>%
  head()

# target variables
vars_of_interest <- c(
  lonely = "lssupvl", # Agreement with I often feel very lonely
  confide = "lssupac", # Agreement with I don't have anyone to confide in
  lean = "lssuplt", # Agreement with I have no one to lean on in times of trouble
  helpless = "lssupnh", # Agremeent with I often need help from people but can't get it
  helped = "lssupsh", # Agreement with When I need someone to help me out I can usually find someone
  health = "losatyh", # Satisfaction with health
  happy = "gh9_sum"
  )


# Bankruptcy
bankrupt.df <- subgroups %>%
  filter(eventcode == "bankruptcy") %>% 
  dplyr::select(xwaveid, response) %>%
  left_join(filter(covariates, code %in% vars_of_interest)) %>%
  mutate(val = ifelse(val < 0, NA_real_, val)) %>% 
  arrange(xwaveid, code, year) %>% 
  group_by(xwaveid, code) %>%
  fill(val) %>% 
  ungroup() %>%
  dplyr::select(response, xwaveid, wave, code, val) %>% 
  na.omit()


create_input <- function(df, subgroup, variables) {
  
  subgroup_df <- filter(df, response == subgroup) %>% dplyr::select(-response)
  
  matrixlist <- list()
  for (variable in variables) {
    subgroup_df %>%
      filter(code == variable) %>% 
      spread(wave, val) %>%
      gather(wave, val, a:p) %>% 
      group_by(xwaveid) %>%
      fill(val, .direction = "downup") %>%
      mutate(sval = as.numeric(scale(val))) %>%
      ungroup() %>%
      mutate(sval = if_else(is.na(sval), 0, sval)) %>%
      dplyr::select(-val, -code) %>%
      spread(wave, sval) %>%
      column_to_rownames("xwaveid") %>% 
      as.matrix() -> mat
    matrixlist[[variable]] <- mat
  }
  
  input.array <- simplify2array(matrixlist)
  return(input.array)
}

responsetype = "relieved"
input.array <- create_input(bankrupt.df, subgroup = responsetype, vars_of_interest)
summary(input.array)
hist(input.array, main = paste("Histogram of", responsetype))

library(EDISON)

# Specify inference with no changepoints
edison.options = defaultOptions()
edison.options$cp.fixed = TRUE

start.time <- Sys.time()
edison.result = EDISON.run(input.array, num.iter=50000, options=edison.options)
end.time <- Sys.time()
print(end.time - start.time) # Time difference of 2.443896 mins

# Calculate posterior probabilities of the edges in the network
ans <- calculateEdgeProbabilities(edison.result)$probs.segs[[1]]
colnames(ans) = rownames(ans) = names(vars_of_interest)
ans

network <- igraph::graph_from_adjacency_matrix(ans > 0.95)
sort(degree(network), decreasing = T)
plot(network, 
     main = paste("Network", responsetype, "after bankruptcy"),
     vertex.size = 40, 
     edge.arrow.size=.75, 
     # edge.curved = .35,
     layout = igraph::layout_in_circle(network)
     # layout = igraph::layout_on_grid(network)
     # layout = igraph::layout_as_star(network)
     )


sort(betweenness(network), decreasing = T)
sort(degree(network, mode = "in"), decreasing = T)
sort(closeness(network), decreasing = T)

sort(
  (degree(network, mode = "in") / 
    degree(network, mode = "out")) * degree(network, mode = "in"), 
  decreasing = T)

ans.relieved <- ans


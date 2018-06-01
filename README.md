Causal network analysis of happiness and wellbeing
================
Richard Morris
01 June, 2018

Introduction
------------

Answering causal questions in the social sciences (e.g., economics, public health) is difficult due to the inability to intervene in the causal variable and measure the impact of the intervention. However with the growing availability of large observational datasets and advances in statistical methods, other approaches to exploring causal questions are increasingly possible. Here we explore the application of a Bayesian network to answer a causal question in the Household Income and Labour Dynamics in Australia (HILDA) survey.

Bayesian networks are a Graphical model which represent the joint probability distribution between the variables of interest. Graphical models represent variables of interest as nodes and the associations between them as arcs or edges. In a Bayesian network, the edges represent the conditional probability of a node given it's parent. That is, a node may be a child of an upstream variable and a parent of a downstream variable. In a Bayesian network, a node may be conditionally dependent on a number of parent nodes, and also be associated with a number of child nodes. From an intuitive point of view, it can be argued that a “good” Bayesian network should represent the causal structure of the data it is describing. Such networks are usually fairly sparse, and their interpretation is at the same time clear and meaningful (Pearl, 2009). However, learning causal models, especially from observational data, presents significant challenges. Three additional assumptions are needed compared to non-causal Bayesian network learning:

-   There must exist a network structure that is faithful to the dependence structure of the causal system
-   There must be no latent variables (unobserved variables in the network) acting as confounding factors
-   Each variable is conditionally independent of its non-effects given its direct causes (this is the Causal Markov Assumption)

It is worth noting that even when dealing with interventional data (a true experiment), these assumptions may not be met. For instance, there are usually multiple equivalent network structures that could explain the data resulting from a true experiment (e.g., a latent variable could be mediating the effect of the intervention observed or a confound may exist). However one way to ensure the third assumption is met (Causal Markov Assumption) in observational data is to incorporate a unidirectional effect of time in the model. This obviously requires a dataset that is collected over multiple timepoints, and the unidirectional effect of time ensures that associations between variables at t2 cannot influence variables at t1. One such dataset is provided by the HILDA survey, which has collected responses from more than 17,000 Australians over 16 years. See <https://melbourneinstitute.unimelb.edu.au/hilda>

The HILDA survey provides a relatively comprehensive assessment of labour force participation, income, family characteristics, as well as personal life events, health and satisfaction. Included in this is the SF-36 health survey, which in turn provides a mental health score derived from agreement with five statements:

> How much of the time during the past 4 weeks have you :
> 1. Been a nervous person
> 2. Felt so down in the dumps nothing could cheer you up
> 3. Felt calm and peaceful
> 4. Felt down
> 5. Been a happy person

This dataset allows us to ask what variables may influence mental health, and vice versa, in the Australian population. Here we apply a **dynamic Bayesian network** analysis to test different causal networks around mental health.

Experiment 1
------------

------------------------------------------------------------------------

There are a variety of variables in the HILDA survey that we could include in our causal model. Initially we will focus on nine variables that represent indicators of **happiness, health, finance, social support, personal control, and relationship satisfaction**. These variables have been selected to represent each category while minimising selection bias as much as possible (e.g., none of the variables depend upon employment or relationship status)

> 1.  'ghmh', MH-I5
> 2.  'losat', life satisfaction
> 3.  'losatyh', health satisfaction
> 4.  'wscei', weekly gross wages/salary
> 5.  'lssupac', I don't have anyone I can confide in
> 6.  'lssupvl', I often feel very lonely
> 7.  'lsrush', feeling rushed
> 8.  'lshrhw', hours per week housework
> 9.  'lsshare', share of housework

In this experiment we will use the ARTIVA package (in R) to estimate the causal structure that may exist between the variables.

``` r
library(tidyverse)
```

    ## ── Attaching packages ────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 2.2.1     ✔ purrr   0.2.4
    ## ✔ tibble  1.4.2     ✔ dplyr   0.7.4
    ## ✔ tidyr   0.8.0     ✔ stringr 1.2.0
    ## ✔ readr   1.1.1     ✔ forcats 0.3.0

    ## ── Conflicts ───────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(ARTIVA)
```

    ## Loading required package: MASS

    ## 
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## Loading required package: igraph

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     as_data_frame, groups, union

    ## The following objects are masked from 'package:purrr':
    ## 
    ##     compose, simplify

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     crossing

    ## The following object is masked from 'package:tibble':
    ## 
    ##     as_data_frame

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

    ## Loading required package: gplots

    ## 
    ## Attaching package: 'gplots'

    ## The following object is masked from 'package:stats':
    ## 
    ##     lowess

``` r
source("src/traceNetwork.R")
source("src/PlotBubbleGraph.R")
```

### Methods

------------------------------------------------------------------------

In Experiment 1 we will determine the causal network between nine variables, among the 3206 people with a complete set of variable scores in each wave. However due to computational limitations, we will randomly sample ***n* = 1000** from the total sample. In this way, a unique set of n = 1000 can be obtained up to three times to repeat the experiment and thus test the **replicability** of the obtained causal network.

``` r
#### load the HILDA dataset ####
library(tidyverse)
library(haven)
source("~/Dropbox/HILDA/src/JoinFrames.R")
source("~/Dropbox/HILDA/src/GetVarsLong.R")
source("~/Dropbox/HILDA/src/ZapLabel.R")

filepaths <- list.files(
  path = '~/Dropbox/HILDA/data',
  pattern = '^Combined.*.dta$',
  full.names = TRUE
)

hildadata <- list()
for (pathtofile in filepaths) {
  df <- read_dta(pathtofile)
  hildadata <- append(hildadata, list(df))
}

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

df.long <- GetVarsLong(hildadata, predictors)
respondants <- unique(df.long$xwaveid)
df.sample <- filter(df.long, xwaveid %in% sample(respondants, 1000))

df.sample %>%
  unite(newvar, waveid, xwaveid) %>%
  spread(newvar, value) %>%
  remove_rownames() %>%
  column_to_rownames('hildacode') %>%
  as.matrix() -> datamatrix

#### Test sample in ARTIVA ####
targets <- predictors
parents <- predictors
N <- length(unique(df.sample$xwaveid))

start_time <- Sys.time()
results <- ARTIVAnet(
  targetData = datamatrix[targets, ],
  parentData = datamatrix[parents, ],
  targetNames = targets,
  parentNames = parents,
  dataDescription = rep(seq(14), each = N),
  dyn = 1, # minimum lag (0 or 1)
  edgesThreshold = 0.5,
  niter = 10000,
  nbCPinit = 0,
  maxCP = 0,
  PSRFactor = FALSE,
  saveEstimations = TRUE,
  saveIterations = TRUE,
  savePictures = TRUE)

end_time <- Sys.time()
print(end_time - start_time)
```

Results
-------

------------------------------------------------------------------------

The estimation of this network with nine variables takes around 6 days to run. The result table is an 81 row by 7 column table listing the posterior probability of each edge `PostProb` (thresholded to edges with a Pr &gt; 0.5), and the causal strength of the edge `CoeffMean`. Negative values indicate an inhibitory or preventative causal influence:

    ##     Parent  Target CPstart CPend PostProb CoeffMean edgesThreshold
    ## 1     ghmh    ghmh       2    14   1.0000   0.53313            0.5
    ## 2    losat    ghmh       2    14   1.0000   0.04855            0.5
    ## 3  losatyh    ghmh       2    14   1.0000   0.08389            0.5
    ## 4    wscei    ghmh       2    14   0.9623   0.02303            0.5
    ## 5  lssupac    ghmh       2    14   1.0000  -0.04273            0.5
    ## 6  lssupvl    ghmh       2    14   1.0000  -0.07836            0.5
    ## 7   lsrush    ghmh       2    14   1.0000   0.05881            0.5
    ## 8   lshrhw    ghmh       2    14   0.5217  -0.01163            0.5
    ## 9  lsshare    ghmh       2    14   0.2458   0.00000            0.5
    ## 10    ghmh   losat       2    14   1.0000   0.08240            0.5
    ## 11   losat   losat       2    14   1.0000   0.47482            0.5
    ## 12 losatyh   losat       2    14   1.0000   0.09849            0.5
    ## 13   wscei   losat       2    14   0.7623  -0.02217            0.5
    ## 14 lssupac   losat       2    14   1.0000  -0.05451            0.5
    ## 15 lssupvl   losat       2    14   1.0000  -0.04611            0.5
    ## 16  lsrush   losat       2    14   1.0000   0.03924            0.5
    ## 17  lshrhw   losat       2    14   0.4721   0.00000            0.5
    ## 18 lsshare   losat       2    14   0.1981   0.00000            0.5
    ## 19    ghmh losatyh       2    14   1.0000   0.08800            0.5
    ## 20   losat losatyh       2    14   1.0000   0.04447            0.5
    ## 21 losatyh losatyh       2    14   1.0000   0.61126            0.5
    ## 22   wscei losatyh       2    14   1.0000   0.03526            0.5
    ## 23 lssupac losatyh       2    14   0.6031  -0.02181            0.5
    ## 24 lssupvl losatyh       2    14   0.8443  -0.02766            0.5
    ## 25  lsrush losatyh       2    14   0.0783   0.00000            0.5
    ## 26  lshrhw losatyh       2    14   0.0376   0.00000            0.5
    ## 27 lsshare losatyh       2    14   0.4413   0.00000            0.5
    ## 28    ghmh   wscei       2    14   0.2710   0.00000            0.5
    ## 29   losat   wscei       2    14   0.4175   0.00000            0.5
    ## 30 losatyh   wscei       2    14   1.0000   0.02922            0.5
    ## 31   wscei   wscei       2    14   1.0000   0.84731            0.5
    ## 32 lssupac   wscei       2    14   0.0204   0.00000            0.5
    ## 33 lssupvl   wscei       2    14   0.0365   0.00000            0.5
    ## 34  lsrush   wscei       2    14   1.0000  -0.03818            0.5
    ## 35  lshrhw   wscei       2    14   1.0000  -0.03452            0.5
    ## 36 lsshare   wscei       2    14   0.2157   0.00000            0.5
    ## 37    ghmh lssupac       2    14   1.0000  -0.04135            0.5
    ## 38   losat lssupac       2    14   1.0000  -0.05670            0.5
    ## 39 losatyh lssupac       2    14   0.9900  -0.02459            0.5
    ## 40   wscei lssupac       2    14   0.5727   0.00325            0.5
    ## 41 lssupac lssupac       2    14   1.0000   0.33432            0.5
    ## 42 lssupvl lssupac       2    14   1.0000   0.08268            0.5
    ## 43  lsrush lssupac       2    14   0.6542  -0.00582            0.5
    ## 44  lshrhw lssupac       2    14   0.6300  -0.00545            0.5
    ## 45 lsshare lssupac       2    14   0.8984   0.01564            0.5
    ## 46    ghmh lssupvl       2    14   1.0000  -0.16776            0.5
    ## 47   losat lssupvl       2    14   1.0000  -0.06720            0.5
    ## 48 losatyh lssupvl       2    14   0.1816   0.00000            0.5
    ## 49   wscei lssupvl       2    14   0.0275   0.00000            0.5
    ## 50 lssupac lssupvl       2    14   1.0000   0.11466            0.5
    ## 51 lssupvl lssupvl       2    14   1.0000   0.34780            0.5
    ## 52  lsrush lssupvl       2    14   0.0221   0.00000            0.5
    ## 53  lshrhw lssupvl       2    14   0.0689   0.00000            0.5
    ## 54 lsshare lssupvl       2    14   1.0000  -0.06211            0.5
    ## 55    ghmh  lsrush       2    14   1.0000   0.05740            0.5
    ## 56   losat  lsrush       2    14   0.3953   0.00000            0.5
    ## 57 losatyh  lsrush       2    14   0.0483   0.00000            0.5
    ## 58   wscei  lsrush       2    14   1.0000  -0.07695            0.5
    ## 59 lssupac  lsrush       2    14   0.0192   0.00000            0.5
    ## 60 lssupvl  lsrush       2    14   0.0201   0.00000            0.5
    ## 61  lsrush  lsrush       2    14   1.0000   0.61607            0.5
    ## 62  lshrhw  lsrush       2    14   1.0000  -0.03994            0.5
    ## 63 lsshare  lsrush       2    14   1.0000   0.04108            0.5
    ## 64    ghmh  lshrhw       2    14   0.0056   0.00000            0.5
    ## 65   losat  lshrhw       2    14   0.2585   0.00000            0.5
    ## 66 losatyh  lshrhw       2    14   0.0344   0.00000            0.5
    ## 67   wscei  lshrhw       2    14   1.0000  -0.08327            0.5
    ## 68 lssupac  lshrhw       2    14   0.0227   0.00000            0.5
    ## 69 lssupvl  lshrhw       2    14   0.0417   0.00000            0.5
    ## 70  lsrush  lshrhw       2    14   0.8903  -0.02724            0.5
    ## 71  lshrhw  lshrhw       2    14   1.0000   0.60045            0.5
    ## 72 lsshare  lshrhw       2    14   1.0000  -0.10479            0.5
    ## 73    ghmh lsshare       2    14   0.0147   0.00000            0.5
    ## 74   losat lsshare       2    14   0.0247   0.00000            0.5
    ## 75 losatyh lsshare       2    14   0.1185   0.00000            0.5
    ## 76   wscei lsshare       2    14   0.9871   0.02794            0.5
    ## 77 lssupac lsshare       2    14   0.0543   0.00000            0.5
    ## 78 lssupvl lsshare       2    14   1.0000  -0.03183            0.5
    ## 79  lsrush lsshare       2    14   1.0000   0.04074            0.5
    ## 80  lshrhw lsshare       2    14   1.0000  -0.09049            0.5
    ## 81 lsshare lsshare       2    14   1.0000   0.67610            0.5

The edges and nodes represent the model structure, while the coefficients represent the model parameters. An easier way to visualize the results is a bubble graph.

Bubble graphs take various forms. Shown below are two variations: **circle** and **Reingold-Tilford**. The circle graph below shows all the nodes and the edges between them with a posterior probability of Pr &gt; 0.95. I.e., we can be very confident that these edges exist in the current dataset.

``` r
# layout = geneLines, random, circle, kamada.kawai, 
#          reingold.tilford, lgl, graphopt,
#          mds, fruchterman.reingold 
traceNetwork(results1, 0.95, layout = 'circle') 
```

![](figures/Exp%201%20circle-1.png)

We can see by the number of arrows pointing to the mental health score `ghmh` that it occupies a central position (i.e., a hub) in this causal network. Network metrics such as centrality and assortivity can be calculated to confirm this impression (to be done). However the circle layout doesn't reveal a great deal of structure in the network. For that, we must turn to other variants such as the Reingold-Tilford.

``` r
traceNetwork(results1, 0.95, layout = 'reingold.tilford') 
```

    ## Warning in layout_as_tree(structure(list(9, TRUE, c(0, 1, 2, 3, 4, 5, 6, :
    ## At structural_properties.c:3338 :graph contains a cycle, partial result is
    ## returned

![](figures/Exp%201%20Reingold-Tilford-1.png)

The Reingold-Tilford graph reveals the mental health score `ghmh` occupies the sole central position in this causal network, under which every other node can be laid in a hierarchy of layers. At the bottom of the network, the relationship variables (based on housework `lsshare`, `lshrhw`) have no direct association with mental health, but every other variable does. However both housework scores indirectly affect mental health via personal control `lsrush` (as does weekly income `wscei`).

#### Experiment 1 replication

------------------------------------------------------------------------

Here we replicate Experiment 1 with a new (independent) sample. Time taken was 6.2 days.

    ##     Parent  Target CPstart CPend PostProb CoeffMean edgesThreshold
    ## 1     ghmh    ghmh       2    14   1.0000   0.55333            0.5
    ## 2    losat    ghmh       2    14   1.0000   0.04077            0.5
    ## 3  losatyh    ghmh       2    14   1.0000   0.08101            0.5
    ## 4    wscei    ghmh       2    14   0.9620   0.02293            0.5
    ## 5  lssupac    ghmh       2    14   1.0000  -0.04094            0.5
    ## 6  lssupvl    ghmh       2    14   1.0000  -0.08220            0.5
    ## 7   lsrush    ghmh       2    14   1.0000   0.05421            0.5
    ## 8   lshrhw    ghmh       2    14   0.1037   0.00000            0.5
    ## 9  lsshare    ghmh       2    14   0.0929   0.00000            0.5
    ## 10    ghmh   losat       2    14   1.0000   0.06279            0.5
    ## 11   losat   losat       2    14   1.0000   0.47140            0.5
    ## 12 losatyh   losat       2    14   1.0000   0.08857            0.5
    ## 13   wscei   losat       2    14   0.9992  -0.02861            0.5
    ## 14 lssupac   losat       2    14   1.0000  -0.03390            0.5
    ## 15 lssupvl   losat       2    14   1.0000  -0.03834            0.5
    ## 16  lsrush   losat       2    14   1.0000   0.04095            0.5
    ## 17  lshrhw   losat       2    14   0.9973   0.02686            0.5
    ## 18 lsshare   losat       2    14   0.9871   0.02148            0.5
    ## 19    ghmh losatyh       2    14   1.0000   0.07433            0.5
    ## 20   losat losatyh       2    14   0.9921   0.03429            0.5
    ## 21 losatyh losatyh       2    14   1.0000   0.59687            0.5
    ## 22   wscei losatyh       2    14   1.0000   0.03521            0.5
    ## 23 lssupac losatyh       2    14   0.0333   0.00000            0.5
    ## 24 lssupvl losatyh       2    14   1.0000  -0.03662            0.5
    ## 25  lsrush losatyh       2    14   0.0227   0.00000            0.5
    ## 26  lshrhw losatyh       2    14   0.0228   0.00000            0.5
    ## 27 lsshare losatyh       2    14   0.7315  -0.02000            0.5
    ## 28    ghmh   wscei       2    14   0.0425   0.00000            0.5
    ## 29   losat   wscei       2    14   0.9777  -0.02540            0.5
    ## 30 losatyh   wscei       2    14   1.0000   0.03360            0.5
    ## 31   wscei   wscei       2    14   1.0000   0.87094            0.5
    ## 32 lssupac   wscei       2    14   0.0112   0.00000            0.5
    ## 33 lssupvl   wscei       2    14   0.0071   0.00000            0.5
    ## 34  lsrush   wscei       2    14   1.0000  -0.02614            0.5
    ## 35  lshrhw   wscei       2    14   1.0000  -0.02839            0.5
    ## 36 lsshare   wscei       2    14   0.0272   0.00000            0.5
    ## 37    ghmh lssupac       2    14   1.0000  -0.06476            0.5
    ## 38   losat lssupac       2    14   1.0000  -0.05883            0.5
    ## 39 losatyh lssupac       2    14   0.0995   0.00000            0.5
    ## 40   wscei lssupac       2    14   0.0119   0.00000            0.5
    ## 41 lssupac lssupac       2    14   1.0000   0.40788            0.5
    ## 42 lssupvl lssupac       2    14   1.0000   0.10407            0.5
    ## 43  lsrush lssupac       2    14   0.0145   0.00000            0.5
    ## 44  lshrhw lssupac       2    14   0.0156   0.00000            0.5
    ## 45 lsshare lssupac       2    14   0.0176   0.00000            0.5
    ## 46    ghmh lssupvl       2    14   1.0000  -0.16988            0.5
    ## 47   losat lssupvl       2    14   1.0000  -0.06119            0.5
    ## 48 losatyh lssupvl       2    14   0.9884  -0.03273            0.5
    ## 49   wscei lssupvl       2    14   0.9983  -0.02952            0.5
    ## 50 lssupac lssupvl       2    14   1.0000   0.10222            0.5
    ## 51 lssupvl lssupvl       2    14   1.0000   0.33876            0.5
    ## 52  lsrush lssupvl       2    14   0.1334   0.00000            0.5
    ## 53  lshrhw lssupvl       2    14   0.1494   0.00000            0.5
    ## 54 lsshare lssupvl       2    14   0.9989  -0.03427            0.5
    ## 55    ghmh  lsrush       2    14   1.0000   0.04731            0.5
    ## 56   losat  lsrush       2    14   0.7912   0.02487            0.5
    ## 57 losatyh  lsrush       2    14   0.0709   0.00000            0.5
    ## 58   wscei  lsrush       2    14   1.0000  -0.07205            0.5
    ## 59 lssupac  lsrush       2    14   0.0369   0.00000            0.5
    ## 60 lssupvl  lsrush       2    14   0.0239   0.00000            0.5
    ## 61  lsrush  lsrush       2    14   1.0000   0.60876            0.5
    ## 62  lshrhw  lsrush       2    14   1.0000  -0.04099            0.5
    ## 63 lsshare  lsrush       2    14   0.0915   0.00000            0.5
    ## 64    ghmh  lshrhw       2    14   0.0819   0.00000            0.5
    ## 65   losat  lshrhw       2    14   1.0000   0.03722            0.5
    ## 66 losatyh  lshrhw       2    14   0.0363   0.00000            0.5
    ## 67   wscei  lshrhw       2    14   1.0000  -0.07973            0.5
    ## 68 lssupac  lshrhw       2    14   0.0303   0.00000            0.5
    ## 69 lssupvl  lshrhw       2    14   0.0467   0.00000            0.5
    ## 70  lsrush  lshrhw       2    14   1.0000  -0.04073            0.5
    ## 71  lshrhw  lshrhw       2    14   1.0000   0.62489            0.5
    ## 72 lsshare  lshrhw       2    14   1.0000  -0.09090            0.5
    ## 73    ghmh lsshare       2    14   0.6468  -0.00774            0.5
    ## 74   losat lsshare       2    14   0.7664   0.01116            0.5
    ## 75 losatyh lsshare       2    14   0.8871  -0.01512            0.5
    ## 76   wscei lsshare       2    14   0.8255   0.01291            0.5
    ## 77 lssupac lsshare       2    14   0.5914   0.00256            0.5
    ## 78 lssupvl lsshare       2    14   0.8584  -0.01300            0.5
    ## 79  lsrush lsshare       2    14   0.6954   0.00874            0.5
    ## 80  lshrhw lsshare       2    14   1.0000  -0.06189            0.5
    ## 81 lsshare lsshare       2    14   1.0000   0.46443            0.5

![](figures/Experiment%201b%20results-1.png)

    ## Warning in layout_as_tree(structure(list(9, TRUE, c(0, 1, 2, 3, 4, 5, 6, :
    ## At structural_properties.c:3338 :graph contains a cycle, partial result is
    ## returned

![](figures/Experiment%201b%20results-2.png)

In comparison to the original sample, there are some differences apparent in the centrality/betweenness of `losat` and `lssupvl`. In the original sample, `ghmh` was clearly the most central variable, receiving the majority of inputs from the other variables. However here `losat` and `lssupvl` have the same number or more inputs than `ghmh`.

However the Reingold-Tilford layout appears to show that `ghmh` occupies the same root position as before, shown here at the top of the graph. In this respect at least, we may have replicated the graph structure with respect to mental health.

Experiment 2
------------

------------------------------------------------------------------------

One assumption of Experiment 1 is the network around mental health is stable across the 14 years of HILDA among our sample. This obviously may not be true for an individual, especially if some personal event has occurred which challenges their life or mental health. A unique feature of ARTIVA is that it can estimate networks which change over time (time-varying networks) by identifying *changepoints* and therefore the different subnetworks on either side of the changepoint. The HILDA dataset offers a unique opportunity to exploit this since it contains self-reported details of ' major life-events' for each respondant in each year. In the present experiment, we identify respondents in HILDA who report a single major life event in the 14 years of the survey, and test whether ARTIVA can detect a changepoint around the time of that life event.

### Methods

------------------------------------------------------------------------

Experiment 2 will gather the individuals who report a single negative life event in the 14 years of the survey and test for changepoints in each individual.

The variables included in the analysis were:

> 1.  'losat', How satisfied are you with your life
> 2.  'losatyh', Your health
> 3.  'losateo', Your employment opportunities
> 4.  'losatfs', Your financial situation
> 5.  'losatft', The amount of free time you have
> 6.  'losathl', The home in which you live
> 7.  'losatlc', Feeling part of your local community
> 8.  'losatnl', The neighbourhood in which you live
> 9.  'losatsf', How safe you feel
> 10. 'lssupnh', I often need help from others
> 11. 'lssupac', I don't have anyone I can confide in
> 12. 'lssuplt', I have no one to lean on
> 13. 'lssupvl', I often feel very lonely
> 14. 'ghmh', mental health

This includes all the life satisfaction variables, as well as social support variables, but not the financial, personal control or relationship variables in Experiment 1 (to be done)

``` r
#### Get the sample of respondents with life events #####
lefnw.df <- JoinFrames(hildadata,'lefnw') # finances worse
lefnw.df[which(lefnw.df < 2, arr.ind=TRUE)] <- NA
lefnw.df[is.na(lefnw.df)] <- 0
lefnw.df[which(lefnw.df == 2, arr.ind=TRUE)] <- 1
colSums(lefnw.df)

ledsc.df <- JoinFrames(hildadata,'ledsc') # death of spouse/child
ledsc.df[which(ledsc.df < 2, arr.ind=TRUE)] <- NA
ledsc.df[is.na(ledsc.df)] <- 0
ledsc.df[which(ledsc.df == 2, arr.ind=TRUE)] <- 1
colSums(ledsc.df)

leins.df <- JoinFrames(hildadata,'leins') # injury to self
leins.df[which(leins.df < 2, arr.ind=TRUE)] <- NA
leins.df[is.na(leins.df)] <- 0
leins.df[which(leins.df == 2, arr.ind=TRUE)] <- 1
colSums(leins.df)

# Collate the events
collated.df <- lefnw.df + leins.df + ledsc.df  # Add all the bad events together
collated.df[collated.df > 0] <- 1 # recode the events to single
colSums(collated.df)
collated.df$sums <- rowSums(collated.df) # Count the events for each person
collated.df <- rownames_to_column(collated.df, 'xwaveid')
sample.df <- filter(collated.df, sums == 1) # Find people with just one bad year

#### Create the variable dataframe ####
predictors <- c('losat',    # How satisfied are you with your life
                'losatyh',  # Your health
                'losateo',  # Your employment opportunities
                'losatfs',  # Your financial situation
                'losatft',  # The amount of free time you have
                'losathl',  # The home in which you live
                'losatlc',  # Feeling part of your local community
                'losatnl',  # The neighbourhood in which you live
                'losatsf',  # How safe you feel
                'lssupnh',  # I often need help from others
                'lssupac',  # I don't have anyone I can confide in
                'lssuplt',  # I have no one to lean on
                'lssupvl',  # I often feel very lonely
                'ghmh'       # mental health
)

df.long <- GetVarsLong(hildadata, predictors, scale = TRUE)
respondants <- unique(df.long$xwaveid)

start_time <- Sys.time()
for (w in letters[3:12]) {
  
  sample.id <- sample.df$xwaveid[sample.df[w] > 0]
  df.sample <- filter(df.long, xwaveid %in% sample.id)
  selected <- unique(df.sample$xwaveid)
  print(paste0('****New wave ', w, ' with ', length(selected), ' people****'))
  
  #### Test sample in ARTIVA ####
  targets <- predictors
  parents <- predictors
  CPresults <- as.tibble()
  
  for (i in seq(length(selected))) {
    
    df.subject <- filter(df.sample, xwaveid == selected[i])
    
    datamatrix <- df.subject %>%
      spread(waveid, value) %>%
      remove_rownames() %>%
      column_to_rownames('hildacode') %>%
      dplyr::select(-xwaveid) %>%
      as.matrix()
    
    results <- ARTIVAnet(
      targetData = datamatrix[targets, ],
      parentData = datamatrix[parents, ],
      targetNames = targets,
      parentNames = parents,
      dyn = 1, # minimum lag (0 or 1)
      segMinLength = 2,
      edgesThreshold = 0.5,
      niter = 10000,
      PSRFactor = FALSE,
      saveEstimations = TRUE,
      savePictures = TRUE
    )
    
    #### Store CP position ####
    CPpaths <- list.files(
      path = 'ARTIVAnet/Estimations',
      pattern = '^CPpositionPostDist.*.txt$',
      full.names = TRUE
    )
    
    for (path in CPpaths) {
      CPtable <- read.csv(path, sep = " ")
      idx_start <- unlist(gregexpr('_', path)) + 1
      idx_end <- unlist(gregexpr('.txt', path)) - 1
      CPtable$var <- substring(path, idx_start, idx_end)
      CPtable$xwaveid <- selected[i]
      CPtable <- rownames_to_column(CPtable, 'wave')
      CPtable$wave <- as.numeric(CPtable$wave) - 7
      CPresults <- rbind(CPresults, CPtable)
    }
    
  }
}
end_time <- Sys.time()
print(end_time - start_time)
write_csv(CPresults, 'results/CPresults.csv')
```

### Results

------------------------------------------------------------------------

The analysis took X days, Y hrs, Z min to complete (to be done). The posterior probability of a change point in each wave for each individual and for each variable was calculated.

    ##    wave Probability    sem   n
    ## 1    -3      0.0874 0.0042 544
    ## 2    -2      0.0604 0.0026 540
    ## 3    -1      0.0690 0.0035 540
    ## 4     0      0.0682 0.0034 538
    ## 5     1      0.0678 0.0035 537
    ## 6     2      0.0702 0.0036 538
    ## 7     3      0.0621 0.0031 543
    ## 8     4      0.0659 0.0032 538
    ## 9     5      0.0604 0.0027 539
    ## 10    6      0.0788 0.0037 543

There is clearly no peak in the average posterior probability of a changepoint near the time of negative life event (wave 0). The *n* in the table above represents the number of datapoints contributing to each mean score. It is a product of the number of variables (e.g., 14) and number of respondents in that wave (e.g., 39).

We can group the results by each variable to see if a peak appears near the life event (wave 0) for any specific variables.

``` r
# Figure 
CPresults %>%
  ggplot(aes(x = wave, y = posterior)) +
    geom_point() +
    geom_smooth(method = 'loess', span = 0.25) +
    facet_wrap(~var) +
    labs(title = 'Change point posterior distribution',
         subtitle = 'Grouped by variable',
         caption = 'ARTIVAnet_bysubject.r')
```

![](figures/Exp%202%20by%20variable-1.png)

There is no specific variable with peaks at or near the time of the life event (wave 0).

We can also filter on the maximum posterior probability of a change point in each variable (and in each individual)

``` r
CPresults %>%
  filter(posterior > 0) %>%
  group_by(xwaveid, var) %>%
  filter(posterior == max(posterior)) %>%
  ggplot(aes(x = wave, y = posterior)) + 
    geom_point() + 
    geom_smooth(method = 'loess', span = 0.45)
```

![](figures/Exp%202%20by%20maximum-1.png)

Here there looks like there is a slight (very slight!) peak in the second wave after the negative life event. To be continued...

#### Set up some working directories ####
setwd("~/Work/HILDA/Causal modelling")
#resultsdir <- ("../results/")
datadir <- ("~/Work/HILDA/super-well/data/")

#### load some libraries ####
library(tidyverse)
library(haven)
library(ARTIVA)
source("~/Work/HILDA/super-well/src/HILDA_fun.r")

#### Load the HILDA dataset ####
filepaths <- list.files(
  path = datadir,
  pattern = '^Combined.*.dta$',
  full.names = TRUE
)

cat('Importing dataframes')
hildadata <- list()
for (pathtofile in filepaths) {
  df <- read_dta(pathtofile)
  hildadata <- append(hildadata, list(df))
  cat('.')
}
cat('Done')

##### Get a list of complete variable dataframes ####
get.vars <- function(df.list, varnames) {
  
  # Collate the variables into a list
  var.list <- list()
  for (variable in varnames) {
    
    df <- joinframes(df.list, variable)
    df[which(df <= 0, arr.ind = TRUE)] <- NA
    df <- df[complete.cases(df),]
    df[,2:length(df)] <- scale(df[,-1])
    colnames(df) <- c('xwaveid',letters[1:length(df)-1])
    var.list <- append(var.list, list(df))
  }
  
  # Find the common xwaveids
  subj.list <- list()
  for (i in seq(length(var.list))) {
    subj.list <- append(subj.list, list(var.list[[i]]['xwaveid']))
  }
  subj.complete <- Reduce(intersect, subj.list)
  
  # Filter the common xwaveids
  for (i in seq(length(var.list))) {
    df <- var.list[[i]]
    subj.ind <- df[['xwaveid']] %in% subj.complete[['xwaveid']]
    df <- df[subj.ind,]
    var.list[i] <- list(df)
  }
  
  # update the list element names
  names(var.list) <- varnames
  
  return(var.list)
}

varname <- c('losatyh', 'ghmh', 'lsrush', 'losat', 'lssupvl', 'lspact')
df.list <- get.vars(hildadata, varname)
complete.subj <- df.list[[1]]['xwaveid'] # table of xwaveids with complete records

##### Calculate change scores for each variable ####

#### Test each subject in ARTIVA ####
targets <- c("ghmh")
parents <- names(df.list)
results <- list()
for (s in seq(3798)) {
  print(s)
  #### recast the variable list into a subject.df ####
  subj.mat <- matrix(nrow = 6, ncol = 14)
  for (i in seq(length(df.list))) {
    subj.mat[i,] <- as.matrix(df.list[[i]][s,2:15])
  }
  row.names(subj.mat) <- names(df.list)
  
  #### ARTIVA #####
  DBN <- ARTIVAnet(
    targetData = subj.mat[targets, ],
    parentData = subj.mat[parents, ],
    targetNames = targets,
    parentNames = parents,
    niter = 50000,
    nbCPinit = 0,  # 0 or 1
    maxCP = 2, # 1 or 2
    savePictures = TRUE)
  
  print(DBN)
  results <- append(results,list(DBN))
  
}

names(results) <- unlist(complete.subj[1:length(results),]) # keep track of xwaveids

#### Store results ####
complete.subj$CPstart <- as.numeric(NA)
complete.subj$CPstart1 <- as.numeric(NA)
complete.subj$CPstart2 <- as.numeric(NA)
complete.subj$CPstart3 <- as.numeric(NA)

for (subj in complete.subj$xwaveid) {
  newrow <- unlist(unique(results[[subj]]['CPstart']))
  complete.subj[complete.subj$xwaveid==subj,names(newrow)] <- newrow
}

#### Get the life events #####
ledsc.df <- joinframes(hildadata,'ledsc') # death of spouse/child
ledsc.df <- column_to_rownames(as.data.frame(ledsc.df), 'xwaveid')
ledsc.df[which(ledsc.df < 2, arr.ind=TRUE)] <- NA
ledsc.df[which(ledsc.df == 2, arr.ind=TRUE)] <- 1
colSums(ledsc.df, na.rm = TRUE)

lefnw.df <- joinframes(hildadata,'lefnw') # finances worse
lefnw.df <- column_to_rownames(as.data.frame(lefnw.df), 'xwaveid')
lefnw.df[which(lefnw.df < 2, arr.ind=TRUE)] <- NA
lefnw.df[which(lefnw.df == 2, arr.ind=TRUE)] <- 1
colSums(lefnw.df, na.rm = TRUE)

leins.df <- joinframes(hildadata,'leins') # injury to self
leins.df <- column_to_rownames(as.data.frame(leins.df), 'xwaveid')
leins.df[which(leins.df < 2, arr.ind=TRUE)] <- NA
leins.df[which(leins.df == 2, arr.ind=TRUE)] <- 1
colSums(leins.df, na.rm = TRUE)

lepcm.df <- joinframes(hildadata,'lepcm') # property crime
lepcm.df <- column_to_rownames(as.data.frame(lepcm.df), 'xwaveid')
lepcm.df[which(lepcm.df < 2, arr.ind=TRUE)] <- NA
lepcm.df[which(lepcm.df == 2, arr.ind=TRUE)] <- 1
colSums(lepcm.df, na.rm = TRUE)

ledhm.df <- joinframes(hildadata,'ledhm') # destroyed home
ledhm.df <- column_to_rownames(as.data.frame(ledhm.df), 'xwaveid')
ledhm.df[which(ledhm.df < 2, arr.ind=TRUE)] <- NA
ledhm.df[which(ledhm.df == 2, arr.ind=TRUE)] <- 1
colSums(ledhm.df, na.rm = TRUE)

#### Collate the negative life events ####
lefnw.df[is.na(lefnw.df)] <- 0
ledsc.df[is.na(ledsc.df)] <- 0
leins.df[is.na(leins.df)] <- 0
lepcm.df[is.na(lepcm.df)] <- 0
ledhm.df[is.na(ledhm.df)] <- 0
collated.df <- lefnw.df + ledsc.df + leins.df# + lepcm.df
collated.df[,8:13] <- collated.df[,8:13] + ledhm.df
colnames(collated.df) <- letters[2:14]
collated.df <- zap_label(collated.df)

#### Compare the occurrence of change points with life events ####
cp.subj <- complete.subj[!is.na(complete.subj$CPstart2),] # get records with change points
cp.subj$waveid <- letters[cp.subj[['CPstart2']]] # get the change point occurence (waveid)
ans <- collated.df[cbind(cp.subj$xwaveid, cp.subj$Col)] # generate answers

cp.subj$colstart <- letters[cp.subj[['CPstart2']]-1]
cp.subj$colend <- letters[cp.subj[['CPstart2']]+1]

ans <- data.frame()
for (i in seq(nrow(cp.subj))) {
  newrow <- collated.df[cp.subj$xwaveid[i], (cp.subj$CPstart2[i]-2):cp.subj$CPstart2[i]]
  names(newrow) <- c('before','during','after')
  ans <- rbind(ans, newrow)
}



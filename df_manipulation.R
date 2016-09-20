setwd('/Users/luxinyi/Desktop/STA440/case1')

# load  & clean data
Brain_network <- load('Brain_network.RData')
na_sub <- as.vector(unique(which(is.na(covariate), arr.ind = TRUE)[,1]))
covariate <- covariate[-na_sub, ] # remove subjects with NA in covariate
W <- W[ , , -na_sub] # remove connectivity matrices associated with those subjects
n_reg <- 68 # total number of brain regions
n_sub <- nrow(covariate) # total number of subjects

# calculate pairwise Euclidean distances
source('Coord_Brain.R')
get_d <- function(u, v, Coord_Brain){
  return (sqrt(sum((Coord_Brain[u, ] - Coord_Brain[v, ])^2)))
}
d <- mapply(get_d, rep(1:n_reg, each = n_reg), rep(1:n_reg, n_reg), MoreArgs = 
              list(Coord_Brain = Coord_Brain))

# calculate indicator for same hem
same_hem <- function(u, v){
  if((u<=34) & (v<=34)) ind = 1
  else if ((u>34) & (v>34)) ind = 1
  else ind = 0
  return (ind)
}
ind <- mapply(same_hem, rep(1:n_reg, each = n_reg), rep(1:n_reg, n_reg))

# reshape dataframe
con <- matrix(0, n_reg^2*n_sub, 9)
colnames(con) <- c('i', 'u', 'v', 'y', 'sex', 'age', 'open', 'd', 'ind')
con[, 'i'] <- rep(1:n_sub, each = n_reg^2)
con[, 'u'] <- rep(1:n_reg, each = n_reg)
con[, 'v'] <- rep(1:n_reg, n_reg)
con[, 'y'] <- as.vector(W)
con[con[, 'y'] > 0, 'y'] = 1 # binary connectivity
con[, 'sex'] <- rep(covariate$Sex, each = n_reg^2)
con[, 'age'] <- rep(covariate$Age, each = n_reg^2)
con[, 'open'] <- rep(covariate$Openness, each = n_reg^2)
con[, 'd'] <- rep(d, n_sub)
con[, 'ind'] <- rep(ind, n_sub)
con <- con[-which(con[,'u']<=con[,'v']), ] # remove upper triangles of matrices of W
setwd('/Users/luxinyi/Desktop/STA440/case1')

df <- load('Brain_network.RData')

get_con <- function(slice) {
  slice[slice > 0] = 1
  n_left <- sum(slice[1:34, 1:34])/2
  n_right <- sum(slice[35:68, 35:68])/2
  n_cross <- sum(slice[1:34, 35:68])
  n_total <- sum(slice)/2
  
  return(c(n_left, n_right, n_cross, n_total))
}


con <- as.data.frame(t(apply(W, 3, get_con)))
colnames(con) <- c('n_left', 'n_right', 'n_cross', 'n_total')

# Optional to include these statistics
# summary(con)

hist(con$n_left, main = 'Number of Connections Within Left-Hemisphere', xlab = 'Count')
hist(con$n_right, main = 'Number of Connections Within Right-Hemisphere', xlab = 'Count')
hist(con$n_cross, main = 'Number of Connections Across Hemispheres', xlab = 'Count')

# Each subject's brain connectivity is represented by a 68 by 68 matrix, of which every 
# element corresponds to the number of fiber connections between two brain regions. 
# A first approach to study the connectivities is to transform the matrices to be binary.
# Since only about 20% of all entries are non-zero, 0 is set as the threshold above which
# an entry is counted as 1. After the transformation, 4 counts are collected for each 
# subject: the number of connections within left hemisphere, within right hemisphere, 
# between hemispheres and overall. Examining the distributions of these values, it can be
# seen that within-hemisphere connections appear normally distributed, whereas the distribution
# of across-hemisphere connections is right-skewed and has a greater variance. 
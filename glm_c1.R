setwd('/Users/luxinyi/Desktop/STA440/case1')

# load data
Brain_network <- load('Brain_network.RData')
na_sub <- as.vector(unique(which(is.na(covariate), arr.ind = TRUE)[,1]))
covariate <- covariate[-na_sub, ] # remove subjects with NA in covariate
W <- W[ , , -na_sub] # remove connectivity matrices associated with those subjects
n_reg <- 68
n_sub <- nrow(covariate)

# calculate pairwise Euclidean distances
source('Coord_Brain.R')
get_d <- function(u, v, Coord_Brain){
   return (sqrt(sum((Coord_Brain[u, ] - Coord_Brain[v, ])^2)))
}
d <- mapply(get_d, rep(1:n_reg, each = n_reg), rep(1:n_reg, n_reg), MoreArgs = 
              list(Coord_Brain = Coord_Brain))

# reshape dataframe
con <- matrix(NA, n_reg^2*n_sub, 8)
colnames(con) <- c('i', 'u', 'v', 'y', 'sex', 'age', 'open', 'd')
con[, 'i'] <- rep(1:n_sub, each = n_reg^2)
con[, 'u'] <- rep(1:n_reg, each = n_reg)
con[, 'v'] <- rep(1:n_reg, n_reg)
con[, 'y'] <- as.vector(W)
con[con[, 'y'] > 0, 'y'] = 1 # binary connectivity
con[, 'sex'] <- rep(covariate$Sex, each = n_reg^2)
con[, 'age'] <- rep(covariate$Age, each = n_reg^2)
con[, 'open'] <- rep(covariate$Openness, each = n_reg^2)
con[, 'd'] <- rep(d, n_sub)
con <- con[-which(con[,'u']<=con[,'v']), ] # remove duplicate pair- and self-connections

# fit model on training set
n_train <- floor(n_sub*9/10) # number of training subjects
set.seed(1234)
i_train <- sample(1:n_sub, n_train, replace = FALSE)
train <- con[which(con[, 'i'] %in% i_train), ]

y <- train[, 'y']
d <- train[, 'd']
sex <- train[, 'sex']
age <- train[, 'age']
open <- train[, 'open']
glm1 <- glm(y ~ d, 
            family = binomial(link = 'logit'))
glm2 <- glm(y ~ sex + age + d, 
            family = binomial(link = 'logit'))
glm3 <- glm(y ~ sex + age + open + d, 
            family = binomial(link = 'logit'))

# AUC (glm2 & 3 are really slow though...)
auc(train[,'y'], glm1$fitted.values) # 0.8281
auc(train[,'y'], glm2$fitted.values) # 0.8281
auc(train[,'y'], glm3$fitted.values) # 0.8282


# fit model on test set
test <- con[-which(con[, 'i'] %in% i_train), ]
pred1 <- predict(glm1, data.frame(d = test[, 'd']), type = 'response')
pred2 <- predict(glm2, data.frame(d = test[, 'd'], sex = test[, 'sex'], age = test[, 'age']), type = 'response')
pred3 <- predict(glm3, data.frame(d = test[, 'd'], sex = test[, 'sex'], age = test[, 'age'], open = test[, 'open']), type = 'response')

auc(test[, 'y'], pred1) # 0.831
auc(test[, 'y'], pred2) # 0.8309
auc(test[, 'y'], pred3) # 0.8308

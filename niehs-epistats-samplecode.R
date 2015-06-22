################
# Setup
# Read in dataset
dat <- read.csv("DataSet2.csv")

# Load libraries
library(psych)
library(splines)
library(rpart)

# Specify covariates
covar <- which(substr(colnames(dat), 1, 1) == "x")




################
# PCA approach
pr1 <- principal(dat[, covar], nfactors = 6)

# Get scores
scores <- pr1$scores
colnames(scores) <- paste0("rotPC", seq(1, 6))
# Scale results
scores <- sweep(scores, 2, apply(scores, 2, sd), "/")
# Create new data
datPC <- data.frame(dat, scores)

# Adjusted linear model
lm1 <- lm(y ~ rotPC1 + rotPC2 + rotPC3 + rotPC4 + rotPC5 + rotPC6 + z1 + z3 + ns(z2, 3), data = datPC)
summary(lm1)$coef





################
# C&RT approach

# First regress out confounders
residxy <- dat
residxy <- apply(residxy, 2, function(x) lm(x ~ z1 + z3 + ns(z2, 3), data = dat)$resid)

# Reformat data
residxy <- residxy[, -which(colnames(residxy) %in% c("z1", "z2", "z3"))]
residxy <- data.frame(residxy)
colnames(residxy) <- paste0("r", colnames(residxy))

# Specify control for C&RT model
fit.control <- rpart.control(xval = 100, cp = 0, minbucket = 5, maxcompete = 4)

# Run C&RT model
whcov <- which(substr(colnames(residxy), 1, 2) == "rx")
eqn1 <- paste("ry ~", paste(colnames(residxy)[whcov], collapse = "+"))
rpart1 <- rpart(eval(eqn1), data = residxy, control = fit.control)


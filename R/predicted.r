# Compute predicted y values

# Notes:
# 1. Takes in results of fitCoeffs()

# Input:
# design = design matrix
# betas = matrix of betas

# Output:
# matrix of predicted values

predicted = function (design, betas) {
  design %*% t(betas)
}

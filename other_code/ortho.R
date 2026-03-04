####
#### R Script to Show that Linear Regression and Projection Approaches
#### Provide the Same Orthogonalized Covariates
####

###
### Simulate (Simple) Deposition Data 
###

## define short-term deposition
n = 30 # number of plots
S = runif(n, 1, 10)

## define long-term deposition
L = S*runif(n, 0.8, 1) + runif(n, 5, 15)*runif(n, S/max(S), 1)
cor(L, S) # simulating with correlation between S and L

## define historic deposition ; just to show lags can still be correlated
H = L - S
cor(H, S) # still correlated

###
### Obtain O_L using regression
###

fit = lm(L ~ S)
O_lm = unname(resid(fit))

###
### Obtain O_L using projection
###

X = cbind(1, S) # design matrix needs to include intercept [1, S]
Phat = X %*% solve(t(X) %*% X, t(X)) # projection onto span{1,S}
M = diag(n) - Phat
O_proj = as.numeric(M %*% L)

###
### Equality Check Between Two Methods
###

all.equal(O_lm, O_proj) # TRUE (up to ~1e-12)
max(abs(O_lm - O_proj)) # ~ machine epsilon

###
### Back-Transformation (if desired)
###

## using regression approach
a_hat = coef(fit)[1]
b_hat = coef(fit)[2]
L_lm = a_hat + b_hat * S + O_lm
all.equal(L_lm, L) # TRUE (up to ~1e-12)

## using projection approach
theta = solve(t(X) %*% X, t(X) %*% L) # same (a_hat, b_hat) via normal equations
L_proj = as.numeric(X %*% theta + O_proj) # fitted + residual
all.equal(L_proj, L) # TRUE (up to ~1e-12)
all.equal(L_proj, L_lm) # TRUE (up to ~1e-12)

## might need to think about fitting this differently to different regions
## if different regions have different degrees of correlation to make sure the
## residuals are 0 for each region - exploratory analysis to see what the strata might be

## would need to do different orthogonalized relationships for every interval b/c have
## different 'historic' values for every interval

## make sure to double-check that correlations are 0 post-orthogonalization
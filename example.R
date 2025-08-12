########################################################################################
#### Example code to train the multivariate Poisson-LogNormal Mixture Representation Model
########################################################################################
library(MASS)
library(lsa)

##### 1. Create simulated dataset ####

# Parameters for simulation
n = 100
N = 10000
p = 800
q = 20
b = 0
set.seed(b)

ar_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}

cor.mat = ar_cor(p, rho=0.5)  # pxp
ev = eigen(cor.mat)
V = ev$vectors
V.tild = V[,1:q]
var.mat = diag(4,q)
L = sqrt(var.mat) %*% ar_cor(q, rho=0.1) %*% sqrt(var.mat)
W = mvrnorm(n+N, rep(0,q), Sigma=L)
B = cbind(rep(0, q), 0.2, rep(0.8, q))
U = log(rpois(n+N, 2)*10 + 1)
group = rbinom(n+N, 1, 0.4)
X = cbind(1, U, group); colnames(X) = c("X1", "X2", "X3")
(X %*% t(B) %*% t(V.tild))[1:3, 1:3]
(W %*% t(V.tild))[1:3, 1:3]
Z = (X %*% t(B) + W) %*% t(V.tild)
Y = apply(Z, c(1,2), FUN=function(z){rpois(1, exp(z))})
X.l = X[1:n,]; Y.l = Y[1:n,]; dat.l = prepare_data(Y.l, X.l)
X.u = X[(n+1):(n+N),]; Y.u = Y[(n+1):(n+N),]; dat.u = prepare_data(Y.u, X.u)
Z.l = Z[1:n,]; Z.u = Z[(n+1):(n+N),]


#### 2. Run PLNfixedVunsup function to train the model ####
lb = list(B = matrix(-Inf, nrow=3, ncol=q),
          L = matrix(0, nrow=q, ncol=q), 
          M = matrix(-Inf, nrow=N, ncol=q),
          S = matrix(0, nrow=N, ncol=q))

pln.fixedv.unsup = PLNfixedVunsup(Abundance ~ X1 + X2, data  = dat.u, control=PLN_param(V=V.tild, rank=q, 
                                                                                        config_optim = list("algorithm" = "MMA",
                                                                                                            "lower_bounds" = lb,
                                                                                                            "ftol_out" = 1e-6,
                                                                                                            "maxit_out" = 10,
                                                                                                            "L0" = L,
                                                                                                            "tau0" = matrix(0.5, ncol=2, nrow=N))))
#### 3. Making predictions of cluster membership ####
fixedv.unsup.pred = predict(pln.fixedv.unsup, dat.u, type="posterior")[,2]


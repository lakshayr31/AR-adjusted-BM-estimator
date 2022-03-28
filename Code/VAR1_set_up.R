library("mvtnorm")
library("vars")
library("fastmatrix")
library("matrixcalc")
#### setting some things upfront
p <- 3
val <- 0.7
T <- 1e4
alpha <- 0.5

phi <- diag(val,p)
W <- diag(alpha,p)

#### First I need to generate a series of samples from a VAR (1)


X0 <- numeric(p)

X <- matrix(0, nrow = p, ncol = T)

X[,1] <- X0


for(i in 2:T){
    X[,i] <- phi %*% X[,i-1] + t(rmvnorm(1,numeric(p),W)) 
}

#### Calculating the actual sigma

V <- solve(diag(p^2) - kronecker.prod(phi,phi)) %*% c(diag(alpha,p))
V_mat <- matrix(V,p,p)

one_mat <- matrix(1,p,p)
act_sigma <- solve(one_mat-phi)%*%V_mat + V_mat%*%solve(one_mat-t(phi)) - V_mat

#### Implementing batch means

b <- floor(T^(1/2))
a <- floor(T/b)
Y <- matrix(0,nrow = p, ncol = a)

for(i in 0:(a-1)){
    Y[,i+1] <- rowMeans(X[,((b*i)+1):(b*(i+1))])
}

theta_hat <- rowMeans(X)

temp_sum <- matrix(0,p,p)

for(i in 1:a){
    temp_sum <- (Y[,i] - theta_hat) %*% t(Y[,i] - theta_hat)
}

bm_cov_matrix <- (b/(a-1))*(temp_sum)

#### Implementing VAR(1) 
data <- t(Y)

model <- VAR(y = data, p = 1,lag.max = NULL)

W_est <- cov(residuals(model))


coeff <- Acoef(model)
coeff_with_error <- Bcoef(model)

phi_est <- coeff[[1]]

V_est <- solve(diag(p^2) - kronecker.prod(phi_est,phi_est)) %*% c(W_est)

V_est_mat <- matrix(V_est,p,p)

one_mat <- matrix(1,p,p)
est_sigma <- solve(one_mat-phi_est)%*%V_est_mat + V_est_mat %*% solve(one_mat-t(phi_est)) - V_est_mat

error <- frobenius.norm(est_sigma-act_sigma)/frobenius.norm(act_sigma)


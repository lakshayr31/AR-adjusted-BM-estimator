library("mvtnorm")
library("vars")
library("fastmatrix")
library("matrixcalc")
library("pracma")

#### setting some things upfront

p <- 5
T <- 1e4
alpha <- 0.4

W <- diag(alpha,p)

#### First I need to generate a series of samples from a VAR (1)

#### Generating a phi

u <- pracma(p, "orthonormal")
sigma <- diag(p)

for(i in 1:p){
    if(i < p/2){
        sigma[i,i] <- runif(1,0,0.5)
    }
    else {
        sigma[i,i] <- runif(1,0.5,1)
    }
}

u <- randortho(p, "orthonormal")
phi <- u %*% sigma %*% t(u)
eigen(phi)$values
# for(x in 1:p)
# {
#     sum = 1
#     for(y in 1:p){
#         curr = runif(1,0,sum)
#         phi[x,y] = curr
#         sum = sum - curr
#     }
# }

iter <- 100

det_sigma_bm <- numeric(iter)
det_sigma_est <- numeric(iter)
det_sigma_act <- numeric(iter)
error_bm <- numeric(iter)
error_est <- numeric(iter)

for(x in 1:iter){
    print(paste0("Iter : ", x))

    X0 <- rnorm(p)

    X <- matrix(0, nrow = p, ncol = T)

    X[,1] <- X0


    for(i in 2:T){
        X[,i] <- phi %*% X[,i-1] + t(rmvnorm(1,numeric(p),W)) 
    }

    V <- solve(diag(p^2) - kronecker.prod(phi,phi)) %*% c(diag(alpha,p))
    V_mat <- matrix(V,p,p)

    one_mat <- diag(1,p)
    act_sigma <- solve(one_mat-phi)%*%V_mat + V_mat%*%solve(one_mat-t(phi)) - V_mat    

    ## Batch Means

    b <- floor(T^(1/2))
    a <- floor(T/b)
    Y <- matrix(0,nrow = p, ncol = a)

    for(i in 0:(a-1)){
        Y[,i+1] <- rowMeans(X[,((b*i)+1):(b*(i+1))])
    }

    theta_hat <- rowMeans(X)

    temp_sum <- matrix(0,p,p)

    for(i in 1:a){
        temp_sum <- temp_sum + (Y[,i] - theta_hat) %*% t(Y[,i] - theta_hat)
    }

    bm_cov_matrix <- (b/(a-1))*(temp_sum)

    bm_error <- frobenius.norm(bm_cov_matrix-act_sigma)/frobenius.norm(act_sigma)

    error_bm[x] <- bm_error

    print(paste0("Batch Means Error ",bm_error))

    data <- t(Y)

    model <- VAR(y = data, p = 1,lag.max = NULL)

    W_est <- summary(model)$covres

    coeff <- Acoef(model)[[1]]

    phi_est <- coeff

    V_est <- solve(diag(p^2) - kronecker.prod(phi_est,phi_est)) %*% c(W_est)

    V_est_mat <- matrix(V_est,p,p)

    one_mat <- diag(1,p)

    est_sigma <- b*solve(one_mat-phi_est)%*%V_est_mat + V_est_mat %*% solve(one_mat-t(phi_est)) - V_est_mat
    act_sigma
    error_est[x] <- frobenius.norm(est_sigma-act_sigma)/frobenius.norm(act_sigma)

    print(paste0("Estimated Error ",error_est[x]))

    det_sigma_act[x] <- determinant(act_sigma, logarithm=TRUE)$modulus
    det_sigma_bm[x] <- determinant(bm_cov_matrix, logarithm=TRUE)$modulus
    det_sigma_est[x] <- determinant(est_sigma, logarithm=TRUE)$modulus
}

plot(density(det_sigma_est), lwd = 1, col ="green", main = "Determinants")
# abline(v = sigma_true[x/2 + 1])
lines(density(det_sigma_bm), lwd = 1, col = "blue")
lines(density(det_sigma_act), lwd = 1, col = "red")

legend("topright", legend=c("ACT ","BM ", "EST "), col=c("red","blue", "green"), lty=1:1:1, cex=0.8)

plot(density(error_bm), main = "Error", lwd = 1, col = "red")
lines(density(error_est), lwd = 1, col = "blue")

legend("topright", legend=c("BM ","EST"), col=c("red","blue"), lty=1:1, cex=0.8)

# A lot of questions in this code 
# 1 Does the Phi matrix need to be similar to a covaraince matrix?
# 2 I am not getting similar values, and I can't figure out what am I doing wrong?
# 3 Is frobenius norm the onlu way to measure out performance, no other comparision or evalutaion metric available?


library("mvtnorm")
library("vars")
library("fastmatrix")
library("matrixcalc")

alpha <- 0.3
p <- 3
W <- diag(alpha,p)

phi <- matrix(0, nrow = p, ncol = p)
for(i in 1:3)
{
    for(x in 1:p)
    {
        sum = 1
        for(y in 1:p){
            curr = runif(1,0,sum)
            phi[x,y] = curr
            sum = sum - curr
        }
    }
    X0 <- numeric(p)

    X <- matrix(0, nrow = p, ncol = T)

    X[,1] <- X0


    for(i in 2:T){
        X[,i] <- phi %*% X[,i-1] + t(rmvnorm(1,numeric(p),W)) 
    }

    V <- solve(diag(p^2) - kronecker.prod(phi,phi)) %*% c(diag(alpha,p))
    V_mat <- matrix(V,p,p)

    one_mat <- diag(1,p)
    act_sigma <- solve(one_mat-phi)%*%V_mat + V_mat%*%solve(one_mat-t(phi)) - V_mat    

    print(paste0("Iteration No ", i))
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

    error <- frobenius.norm(est_sigma-act_sigma)/frobenius.norm(act_sigma)

    print(paste0("Estimated Error ",error))
    
}

# A lot of questions in this code 
# 1 Does the Phi matrix need to be similar to a covaraince matrix?
# 2 I am not getting similar values, and I can't figure out what am I doing wrong?
# 3 Is frobenius norm the onlu way to measure out performance, no other comparision or evalutaion metric available?


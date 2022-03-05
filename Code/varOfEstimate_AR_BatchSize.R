
# When I use this function then it doesnt give the right answer
calculate_est <- function(b, alpha_sq, rho) {
    return ((b*alpha_sq)/((1-rho)^2))
}

# b_dist <- c(100,200,500,1000);
# b_dist <- c(10,20,50,100);
b_dist <- numeric(2);
iter <- 100

sample_mean_dist <- numeric(length = iter)

bm_sigma_dist <- matrix(0, nrow = length(b_dist), ncol = iter)

est_sigma_dist <- matrix(0, nrow = length(b_dist), ncol = iter)
rho_dist <- matrix(0,nrow = length(b_dist), ncol = iter)
alpha_dist <- matrix(0,nrow = length(b_dist), ncol = iter)

est_sigma_dist_aic <- matrix(0, nrow = length(b_dist), ncol = iter)
rho_dist_aic <- matrix(0,nrow = length(b_dist), ncol = iter)
alpha_dist_aic <- matrix(0,nrow = length(b_dist), ncol = iter)

T <- 1e5
b_dist[2] <- floor(T^(1/2))
b_dist[1] <- floor(T^(1/3))

rho <- 0.95
alpha <- 1

sigma_true = (alpha^2)/(1-rho^2)

# This is to the code for which model has to be made for an AR1 process
for(it in 1:iter){
    print(paste("Iter : ", it))

    x <- numeric(T)
    x[1] <- rnorm(1,0,1)

    for(i in 2:T){
        x[i] <- x[i-1]*rho + rnorm(1,0,alpha)
    }

    n <- T

    for(b in b_dist){
        a <- floor(n/b)
        Y <- numeric(a)

        for(i in 0:(a-1)){
            Y[i+1] <- mean(x[((b*i)+1):(b*(i+1))])
        }

        mu <- mean(x)
        bm_sigma <- (b/(a-1))*(sum((Y-mu)^2))

        # AIC is FALSE
        model <- ar(Y, aic = FALSE, order.max = 1)

        rho_c <- model$ar
        alpha_sq_c <- model$var.pred
        
        est_sigma <- calculate_est(b, alpha_sq_c, rho_c)

        rho_dist[which(b == b_dist), it] <- rho_c
        alpha_dist[which(b == b_dist), it] <- alpha_sq_c

        bm_sigma_dist[which(b == b_dist), it] <- bm_sigma
        est_sigma_dist[which(b == b_dist), it] <- est_sigma

        # AIC is TRUE
        model_aic <- ar(Y, aic = TRUE, order.max = 1)

        if(length(model_aic$ar) == 0){
            rho_c_aic = 0
        } else {
            rho_c_aic <- model_aic$ar
        }
    
        alpha_sq_c_aic <- model_aic$var.pred

        rho_dist_aic[which(b == b_dist), it] <- rho_c_aic
        alpha_dist_aic[which(b == b_dist), it] <- alpha_sq_c_aic

        est_sigma_aic <- calculate_est(b, alpha_sq_c_aic, rho_c_aic)

        est_sigma_dist_aic[which(b == b_dist), it] <- est_sigma_aic
    }
}

cbind(apply(bm_sigma_dist,1,mean),apply(est_sigma_dist,1,mean),apply(est_sigma_dist_aic,1,mean))
cbind(apply(bm_sigma_dist,1,var),apply(est_sigma_dist,1,var),apply(est_sigma_dist_aic,1,var)) 

for(i in 1:length(b_dist)){
    mse_bm = (1/n)*sum((sigma_true-bm_sigma_dist[i,])^2);
    mse_est = (1/n)*sum((sigma_true-est_sigma_dist[i,])^2);
    mse_est_aic = (1/n)*sum((sigma_true-est_sigma_dist_aic[i,])^2);
    print(paste("Batch Size ",b_dist[i]," MSE BM: ",mse_bm," MSE EST: ",mse_est,"MSE AIC EST: ", mse_est_aic));
}

par(mfrow=c(1,2))    

for(i in 1:length(b_dist)) {
    plot(density(bm_sigma_dist[i,]), xlim = c(min(bm_sigma_dist[i,],est_sigma_dist[i,],est_sigma_dist_aic[i,]) - 10, max(bm_sigma_dist[i,],est_sigma_dist[i,],est_sigma_dist_aic[i,]) + 10), main = paste("Batch Size : ",b_dist[i]), lwd = 1, col = "red")
    abline(v = 400)
    lines(density(est_sigma_dist[i,]), lwd = 1, col = "blue")
    lines(density(est_sigma_dist_aic[i,]), lwd = 1, col ="green")
    legend("topright", legend=c("BM ","Estimated ", "Estimated AIC "), col=c("red","blue", "green"), lty=1:1:1, cex=0.8)
}

for(i in 1:length(b_dist)) {
    plot(density(rho_dist[i,]), col = "red", main = paste("Rho ", mean(rho_dist[i,]),"BS : ", b_dist[i]))
}

par(mfrow=c(1,1))



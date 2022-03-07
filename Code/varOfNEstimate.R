calculate_est <- function(b, alpha_sq, rho) {
    return ((b*alpha_sq)/((1-rho)^2))
}
# running the procedure of using MCMC to generate samples
# then finding the variance and plotting there density

# defining different targets

# exponential
rate = 2
exponential <- function(x,rate){
    return(dexp(x,rate))
}
# gamma
shape = 9
rate = 0.5
gamma <- function(x,shape,rate){
    return(dgamma(x,shape = shape, rate = rate))
}

# laplace
theta = 2
s = 4    
laplace <- function(x,theta,s){
    return(0.5*(1/s)*exp(-abs(x-theta)/s))
}
# chi-squared

k = 5
chisq <- function(x,k){
    return(dchisq(x,k))
}

# Using normal distribution as proposal 
normal <- function(x,mu,h){
	return (rnorm(x,mean = mu, sd = h))
}

#setting up parameters for distributions
target <- laplace
proposal <- normal

# setting up for running multiple mcmc
iter <- 100

b_dist <- c(100,200,500,1000);
# b_dist <- c(10,20,50,100);

#for the sigma of each run
est_sigma_dist <- matrix(0, nrow = length(b_dist), ncol = iter)
bm_sigma_dist <- matrix(0, nrow = length(b_dist), ncol = iter)
#for the sample mean of each run
sample_mean_dist = numeric(length = iter)

rho_dist <- matrix(0,nrow = length(b_dist), ncol = iter)
alpha_dist <- matrix(0,nrow = length(b_dist), ncol = iter)

est_sigma_dist_aic <- matrix(0, nrow = length(b_dist), ncol = iter)
rho_dist_aic <- matrix(0,nrow = length(b_dist), ncol = iter)
alpha_dist_aic <- matrix(0,nrow = length(b_dist), ncol = iter)


acc_dist = numeric(length = iter)
# running mulitple mcmc
# setting up changing parameters

T <- 1e4
# h <- 1
h <- 5

ini <- 2

for(it in 1:iter){
    print(paste("Iter : ", it))

    x <- numeric(T)
    acc <- 0
    x[1] <- ini + rnorm(1) 

    for(t in 2:T){
        y = proposal(1,x[t-1],h)
        alpha = min((target(y,theta,s)/target(x[t-1],theta,s)),1)
        u = runif(1,0,1)
        if(u <= alpha){
            x[t] = y
            acc <- acc + 1
        } else {
            x[t] = x[t-1]
        }
    }
    acc_dist[it] <- ((acc/T)*100)
    sample_mean_dist[it] = mean(x)

    n <- T

    for(b in b_dist){
        a <- n/b
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

mean(acc_dist)
var(sample_mean_dist)*T

cbind(apply(bm_sigma_dist,1,mean),apply(est_sigma_dist,1,mean),apply(est_sigma_dist_aic,1,mean))
cbind(apply(bm_sigma_dist,1,var),apply(est_sigma_dist,1,var),apply(est_sigma_dist_aic,1,var)) 
 
par(mfrow=c(2,2))    

for(i in 1:length(b_dist)) {
    plot(density(bm_sigma_dist[i,]), xlim = c(min(bm_sigma_dist[i,],est_sigma_dist[i,],est_sigma_dist_aic[i,]) - 10, max(bm_sigma_dist[i,],est_sigma_dist[i,],est_sigma_dist_aic[i,]) + 10), main = paste("Batch Size : ",b_dist[i]), lwd = 1, col = "red")
    abline(v = 400)
    lines(density(est_sigma_dist[i,]), lwd = 1, col = "blue")
    lines(density(est_sigma_dist_aic[i,]), lwd = 1, col ="green")
    legend("topright", legend=c("BM ","Estimated ", "Estimated AIC "), col=c("red","blue", "green"), lty=1:1:1, cex=0.8)
}

par(mfrow=c(2,2))
acf(Y_1,lag.max = 100,main=paste("Batch Size ",b_dist[1]))
acf(Y_2,lag.max = 100,main=paste("Batch Size ",b_dist[2]))
acf(Y_3,lag.max = 100,main=paste("Batch Size ",b_dist[3]))
acf(Y_4,lag.max = 100,main=paste("Batch Size ",b_dist[4]))
par(mfrow=c(1,1))


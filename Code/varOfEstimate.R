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
iter <- 400

b_dist <- c(50,100,200,500);

#for the sigma of each run
est_sigma_dist <- matrix(0, nrow = length(b_dist), ncol = iter/length(b_dist))
bm_sigma_dist <- matrix(0, nrow = length(b_dist), ncol = iter/length(b_dist))
#for the sample mean of each run
sample_mean_dist = numeric(length = iter)

acc_dist = numeric(length = iter)
# running mulitple mcmc
# setting up changing parameters

T <- 1e4

Y_1 = numeric(length = T/b_dist[1])
Y_2 = numeric(length = T/b_dist[2])
Y_3 = numeric(length = T/b_dist[3])
Y_4 = numeric(length = T/b_dist[4])

h <- 10

ini <- 2
it <- 1

while(it <= iter){
    x <- numeric(T)
    acc <- 0
    x[1] <- ini + rnorm(1) #rnorm for a little variance in starting value
    while(x[1] < 0){
        x[1] <- ini + rnorm(1)
    }
    for(t in 2:T){
        y = proposal(1,x[t-1],h)
        # sample space restrictions
        # while(y < 0){ 
        #     y = proposal(1,x[t-1],h)
        # }
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
    b <- b_dist[if (it%%length(b_dist) == 0) length(b_dist) else it%%length(b_dist)]

    a <- n/b
    mu = mean(x)
    Y = numeric(a)
    for(i in 0:(a-1)){
        Y[i+1] = mean(x[((b*i)+1):(b*(i+1))])
    }
    if(b == b_dist[1]){
        Y_1 = Y
    }
    else if (b == b_dist[2]){
        Y_2 = Y
    }
    else if (b == b_dist[3]){
        Y_3 = Y
    }
    else {
        Y_4 = Y
    }
    bm_sigma= (b/(a-1))*(sum((Y-mu)^2))

    model = ar(Y, aic = FALSE,order.max = 1)

    rho <- model$ar
    alpha_sq <- model$var.pred

    est_sigma <- (b*alpha_sq)/((1-rho^2))
    bm_sigma_dist[if (it%%length(b_dist) == 0) length(b_dist) else it%%length(b_dist), if(floor(it/length(b_dist)) == iter/length(b_dist)) iter/length(b_dist) else floor((it-1)/length(b_dist)) + 1] = bm_sigma
    est_sigma_dist[if (it%%length(b_dist) == 0) length(b_dist) else it%%length(b_dist), if(floor(it/length(b_dist)) == iter/length(b_dist)) iter/length(b_dist) else floor((it-1)/length(b_dist)) + 1] = est_sigma
    # print(paste("Row, Column Value", if (it%%length(b_dist) == 0) length(b_dist) else it%%length(b_dist), " ", if(floor(it/length(b_dist)) == iter/length(b_dist)) iter/length(b_dist) else floor((it-1)/length(b_dist)) + 1 , " ", b_dist[if (it%%length(b_dist) == 0) length(b_dist) else it%%length(b_dist)]))
    it <- it + 1
}
plot(density(x))
mean(x)
mean(acc_dist)

var(sample_mean_dist)
var(sample_mean_dist)*T

cbind(apply(bm_sigma_dist,1,mean),apply(est_sigma_dist,1,mean))
cbind(apply(bm_sigma_dist,1,var),apply(est_sigma_dist,1,var))
par(mfrow=c(1,1))
acf(x,lag.max = 100)

par(mfrow=c(2,2))    

for(i in 1:4) {
    plot(density(est_sigma_dist[i,]),xlim=c(0,600), main = paste("Batch Size : ",b_dist[i]), lwd = 2, col = "red")
    lines(density(bm_sigma_dist[i,]), lwd = 2, col = "blue")
    legend("topright", legend=c("Estimated ","BM "), col=c("red","blue"), lty=1:1, cex=0.8)
}
par(mfrow=c(2,2))
acf(Y_1,lag.max = 100,main=paste("Batch Size ",b_dist[1]))
acf(Y_2,lag.max = 100,main=paste("Batch Size ",b_dist[2]))
acf(Y_3,lag.max = 100,main=paste("Batch Size ",b_dist[3]))
acf(Y_4,lag.max = 100,main=paste("Batch Size ",b_dist[4]))



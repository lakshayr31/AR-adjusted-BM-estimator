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

#for the sigma of each run
sigma_dist = numeric(length = iter)
bm_sigma_dist = numeric(length = iter)
#for the sample mean of each run
sample_mean_dist = numeric(length = iter)

acc_dist = numeric(length = iter)
# running mulitple mcmc
# setting up changing parameters

b_dist <- c(50,100,200,400);

b_sigma <- matrix(0, nrow = length(b_dist), ncol = iter/length(b_dist))
T <- 1e4

h <- 8

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
        while(y < 0){ 
            y = proposal(1,x[t-1],h)
        }
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
    bm_sigma_dist[it] = (b/(a-1))*(sum((Y-mu)^2))

    model = ar(Y, aic = FALSE,order.max = 1)
    if(model$order == 0){
        print("Got a model with order 0")
        break
    }
    rho <- model$ar
    alpha_sq <- model$var.pred
    sigma_dist[it] <- (b*alpha_sq)/((1-rho^2))
    b_sigma[if (it%%length(b_dist) == 0) length(b_dist) else it%%length(b_dist), if(floor(it/length(b_dist)) == 25) 25 else floor(it/length(b_dist)) + 1] = (b*alpha_sq)/((1-rho^2))
    it <- it + 1
}

mean(acc_dist)

var(sample_mean_dist)
var(sample_mean_dist)*T

mean(sigma_dist)
mean(bm_sigma_dist)
var(sigma_dist)
var(bm_sigma_dist)


plot(density(sigma_dist), main = "Density Plot")
lines(density(bm_sigma_dist),col = "red")
abline(v=mean(sigma_dist), col="black")
abline(v=mean(bm_sigma_dist), col="red")
abline(v=var(sample_mean_dist)*T, col = "blue")
legend("topright", legend=c("Variance Estimated", "Batch Means Variance"),
       col=c("black", "red"), lty=1:1, cex=0.8)
legend("right", legend=c("mean(sigma_dist)","mean(bm_sigma_dist)","var(sample_mean_dist)*T"),col=c("black","red","blue"), lty=1:1:1, cex=0.8)

col = c("red","blue","green","black")
plot(density(b_sigma[1,]), lwd = 1, col = col[1])
for(x in 2:4) {
    lines(density(b_sigma[x,]), lwd = x, col = col[x])
}

legend("topright", legend=c("B = 50", "B = 100", "B = 200", "B = 400"),
       col=col, lty=1:1:1:1, cex=0.8)
for(x in 1:4) {
    abline(v=mean(b_sigma[x,]), col=col[x])
}

for(x in 1:4) {
    cat(b_dist[x], " ", mean(b_sigma[x,])," \n")
}
for(x in 1:4) {
    cat(b_dist[x], " ", var(b_sigma[x,])," \n")
}
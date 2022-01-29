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

theta = 0
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

{ iter <- 100

#for the sigma of each run
sigma_dist = numeric(length = iter)
bm_sigma_dist = numeric(length = iter)
#for the sample mean of each run
sample_mean_dist = numeric(length = iter)

acc_dist = numeric(length = iter)
# running mulitple mcmc
# setting up changing parameters
b_dist <- c(50,100,200,500);

T <- 1e4

h <- 1.2

ini <- 0
count <- 0

for(it in 1:iter){
    x <- numeric(T)
    acc <- 0
    x[1] <- ini + rnorm(1) #rnorm for a little variance in starting value
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
    b <- b_dist[it%%length(b_dist) + 1]
    a <- n/b
    mu = mean(x)
    Y = numeric(a)

    for(i in 0:(a-1)){
        Y[i+1] = mean(x[((b*i)+1):(b*(i+1))])
    }
    bm_sigma_dist[it] = (b/(a-1))*(sum((Y-mu)^2))

    model = ar(Y, order.max = 1)
    if(model$order == 0){
        count = count + 1
        sigma_dist[it] = bm_sigma_dist[it]
    } else {
        rho <- model$ar
        alpha_sq <- model$var.pred

        sigma_dist[it] <- (b*alpha_sq)/((1-rho^2))
    }
}
}

var(sample_mean_dist)

plot(density(sigma_dist))
lines(density(bm_sigma_dist),col = "red")

var(sigma_dist)
mean(sigma_dist)
var(bm_sigma_dist)
mean(bm_sigma_dist)

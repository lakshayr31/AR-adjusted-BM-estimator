# Setting up general code for MCMC error analysis

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

#normal
mu = 10
s_2 = 2
normal_d <- function(x,mu,s_2){
    return(dnorm(x,mean = mu, sd = s_2^(0.5)))
}
# Using normal distribution as proposal 
normal <- function(x,mu,h){
	return (rnorm(x,mean = mu, sd = h))
}

#setting up parameters for distributions
target <- exponential
proposal <- normal

# Samples 
T <- 1e5
# Step size
h <- 1
# Starting value
acc <- 0 
x <- numeric(T)
x[1] <- 2

for(t in 2:T){
    y = rnorm(1,x[t-1],h)
    while(y <= 0){
        y = rnorm(1,x[t-1],h)
    }
    alpha = min((target(y,rate)/target(x[t-1],rate)),1)
    u = runif(1,0,1)
    if(u <= alpha){
        x[t] = y
        acc <- acc + 1
    } else {
        x[t] = x[t-1]
    }
}
accP = (acc/T)*100
print(paste0("Accuracy ",accP))

# Analysis of Samples
# sample_mean <- mean(x)
# sample_mean
# sample_var <- var(x)
# sample_var
# plot(density(x), lwd = 5)
# lines((density(rexp(T,rate = rate))), lwd= 2, col = 'red')
# legend("topright", legend=c("MCMC Sample Density", "Actual Sample Density"),
#        col=c("black", "red"), lty=1:1, cex=0.8)
# Batch Means 
n <- T
b <- 100
a <- n/b
mu = mean(x)
Y = numeric(length = a)

for(i in 0:(a-1)){
    Y[i+1] = mean(x[((b*i)+1):(b*(i+1))])
}
# plot(density(Y))
# mean(Y)
# var(Y)
#Analysis of Batch Means
bm_sample_var <- (b/(a-1))*(sum((Y-mu)^2))
bm_sample_var
print(paste0("Batch Means Sample Var: ",bm_sample_var))

# AR(1) process
model = ar(Y, aic = FALSE, order.max = 1)
model
rho <- model$ar
rho
alpha_sq <- model$var.pred
alpha_sq
acf(Y,main = "Batch Means Correlation")
# estimator estimate of variance
estimated_sample_var = (b*alpha_sq)/((1-rho^2))

h
accP
bm_sample_var
estimated_sample_var
rho
# Printing Values Properly
cat(" Sample Variance : ", sample_var, "\n", "Batch Means Sample Variance ", bm_sample_var, "\n", "Estimator Variance ", estimated_sample_var, "\n")


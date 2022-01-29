# Setting up general code for MCMC error analysis

# defining different targets

# exponential
rate = 2
exp <- function(x,rate){
    return(dexp(x,rate))
}
# gamma
shape = 9
rate = 0.5
gamma <- function(x,shape,rate){
    return(dgamma(x,shape = shape, rate = rate))
}
# laplace
theta = 5
laplace <- function(x,theta = theta){
    return(0.5*exp(-abs(x-theta)))
}
# chi-squared

k = 5
chisq <- function(x,k){
    return(dchisq(x,k))
}

# Using normal distribution as proposal 
normal <- function(x,mu,h){
	return (rnorm(x,mean = mu, sd =     h))
}

#setting up parameters for distributions
target <- chisq
proposal <- normal

# Samples 
T <- 1e5


# Starting value
x[1] <- k

# Step size
h <- 1

acc <- 0 

for(t in 2:T){
    y = proposal(1,x[t-1],h)
    alpha = min((target(y,k)/target(x[t-1],k)),1)
    u = runif(1,0,1)
    if(u <= alpha){
        x[t] = y
        acc <- acc + 1
    } else {
        x[t] = x[t-1]
    }
}

print(paste0("Accuracy ",(acc/T)*100))

# Analysis of Samples
sample_mean <- mean(x)
sample_mean
sample_var <- var(x)
sample_var
plot(density(x))
lines((density(rchisq(T,df = k))), col = 'red')

# Batch Means 
n <- T
b <- 100
a <- n/b
mu = mean(x)
Y = numeric(length = a)

for(i in 0:(a-1)){
    Y[i+1] = mean(x[((b*i)+1):(b*(i+1))])
}
plot(density(Y))

#Analysis of Batch Means
bm_sample_var <- (b/(a-1))*(sum((Y-mu)^2))
bm_sample_var
print(paste0("Batch Means Sample Var: ",bm_sample_var))

# AR(1) process
model = ar(Y, order.max = 1)
model
rho <- model$ar
rho
alpha_sq <- model$var.pred
alpha_sq
acf(Y,main = "Batch Means Correlation")
# estimator estimate of variance
estimated_sample_var = (b*alpha_sq)/((1-rho^2))
bm_sample_var
estimated_sample_var

# Printing Values Properly
cat(" Sample Variance : ", sample_var, "\n", "Batch Means Sample Variance ", bm_sample_var, "\n", "Estimator Variance ", estimated_sample_var, "\n")


library("languageserver")
library("httpgd")

#simulating the metropolis hasting's algorithm for 
# target distribution cauchy - C(s,t) where s - scale paramter and t - location parameter 
# proposal distribution normal - N(x,h) where x is the sample taken and h is the step size, or the amount of variability needed in the proposal 

cauchy <- function(x,t,s) {
	return (s*(1/(pi*(s^2 + (x-t)^2))))
}

normal <- function(x,mu,h){
	term1 = 1/((2*pi*(h^2))^(0.5))
	term2 = exp((-1/(2*(h^2)))*(x-mu)^2)
	return (term1*term2)
}

# cauchy(1,0,1)
# dcauchy(1,0,1)

# normal(1,0,2^(0.5))
# dnorm(1,0,2^(0.5))

target <- cauchy
t = 0
s = 2

proposal <- normal
h = 1/(2^(0.5))
#metropolis hastings implementation

T = 500
x = numeric(length = T)
x[1] = 0

for(t in 2:T){
	y = rnorm(1, x[t-1], h)
	alpha = min((dcauchy(y,t,s)/dcauchy(x[t-1],t,s)),1)
	u = runif(1, 0, 1)
	if(u <= alpha){
		x[t] = y
	} else {
		x[t] = x[t-1]
	}
}
plot(density(x))
x_test = rcauchy(T,t,s)
x_test


# try with an exponential distribution with rate = 2, and uniform distribution (x-h,x+h)

h = 1
rate = 1

T = 1e4
x = numeric(length = T)
x[1] = 3
acc = 0

for(t in 2:T){
	# y = runif(1,x[t-1]-h,x[t-1]+h)
	y = x[t-1] + rnorm(1,mean = 0, sd = h)
	alpha = min((dexp(y,rate)/dexp(x[t-1],rate)),1)
	u = runif(1,0,1)
	if(u < alpha){
		x[t] = y
		acc = acc + 1
	} else{
		x[t] = x[t-1]
	}
}

(acc/T)*100
mean(x)
x_test = rexp(T,rate)
plot(density(x_test))
plot(density(x))


acf(x, main = "ACF Plot")
plot.ts(x, main = "Trace Plot")



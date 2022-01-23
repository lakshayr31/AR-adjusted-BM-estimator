# MH Algorithm
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
		x[t-1] = y
	}
}

(acc/T)*100
mean(x)
1/rate
x_test = rexp(T,rate)
plot(density(x))
plot(density(x_test))

acf(x, main = "ACF Plot")
plot.ts(x, main = "Trace Plot")

# The mean seems to be really off, whereas the ACF and Trace plots are nice.

# BatchMeans

a = 100
b = 100
N = a*b

mu = mean(x)
mu
Y = numeric(length = b) 

for(i in 0:(a-1)){
    Y[i+1] = mean(x[(b*i+1):(b*(i+1))])
}

var = (b/(a-1))*(sum((Y-mu)**2))
var

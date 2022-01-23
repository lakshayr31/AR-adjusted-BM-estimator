a = 5
b = 10
N = a*b

X = rnorm(N, mean = 5, sd = 2)
X
mu = mean(X)
mu
Y = numeric(length = 5) 
Y
for(i in 0:(a-1)){
    print(mean(X[(b*i+1):(b*(i+1))]))
    Y[i+1] = mean(X[(b*i+1):(b*(i+1))])
}
Y
var = (b/(a-1))*(sum((Y-mu)**2))
var

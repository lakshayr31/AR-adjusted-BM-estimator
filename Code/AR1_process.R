# AR1 process generation

# The parameters for an AR1 process are : rho and alpha


# The formula is X_n+1 = rho*X_n + Y_n ~ N(0,alpha^2)

# How to choose the first term, through sampling from a distribution

# Using a normal distribution ~ N(0,1)

# Setting up 

T <- 1e3
rho_dist <-  c(0.5,0.6,0.7,0.8,0.9)
rho_dist_high <- c(0.95,0.96,0.97,0.98,0.99)


par(mfrow=c(3,2))
# for(rho in rho_dist){
for(rho in rho_dist_high){   
    x <- numeric(T)
    x[1] <- rnorm(1,0,1)

    for(i in 2:T){
        x[i] <- x[i-1]*rho + rnorm(1,0,alpha)
    }

    # acf(x, main = paste("Rho : ", rho))

    # plot.ts(x, main = paste("Rho : ", rho))

    # model = ar(x, aic = FALSE, order.max = 1)
    # print(paste("Rho : ",model$ar,"Alpha^2 : ",model$var.pred))
}

rho <- 0.9
alpha <- 1

x <- numeric(T)
x[1] <- rnorm(1,0,1)

for(i in 2:T){
    x[i] <- x[i-1]*rho + rnorm(1,0,alpha)
}

# acf(x, main = paste("Rho : ", rho))

# plot.ts(x, main = paste("Rho : ", rho))

model = ar(x, aic = FALSE, order.max = 1)
print(paste("Rho : ",model$ar,"Alpha^2 : ",model$var.pred))

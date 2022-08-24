log_f <- function(beta)
{
  one.minus.yx <- (1 - y)*X	
  return(-sum(beta^2)/2 - sum(log(1 + exp(-X%*%beta))) - sum(one.minus.yx%*%beta))
}

bayes_logit_mh <- function(y, X, N = 1e4, prop.sd = .35)
{
  p <- dim(X)[2]
  one.minus.yx <- (1 - y)*X	
  
  # starting value is the MLE
  foo <- glm(y ~X - 1, family = binomial("logit"))$coef
  beta <- as.matrix(foo, ncol = 1)
  beta.mat <- matrix(0, nrow = N, ncol = p)
  beta.mat[1, ] <- as.numeric(beta)
  accept <- 0

  for(i in 2:N)
  {
    #symmetric density
    prop <- rnorm(p, mean = beta, sd = prop.sd)
    
    # log of the MH ratio
    log.rat <- log_f(prop) - log_f(beta)
    if(log(runif(1)) < log.rat)
    {
      beta <- prop
      accept <- accept + 1
    }
    beta.mat[i, ] <- beta
  }
  print(paste("Acceptance Prob = ", accept/N))
  return(beta.mat)
}

titanic <- read.csv("https://dvats.github.io/assets/titanic.csv")

y <- titanic[,1]
X <- as.matrix(titanic[, -1])

## acceptance is too low! we want 23%
## so decrease proposal variance
chain <- bayes_logit_mh(y = y, X = X, N = 1e4, prop.sd = .005) 

# par(mfrow=c(3,2))
# plot.ts(chain[,1])
# plot.ts(chain[,2])
# plot.ts(chain[,3])
# plot.ts(chain[,4])
# plot.ts(chain[,5])
# plot.ts(chain[,6])

b_dist <- numeric(2);

# n_sd = 3
n_T = 3

iter <- 100

#for plotting everything on the same graph
bm_sigma_dist_phi <- matrix(0, nrow = n_T*length(b_dist), ncol = iter)
est_sigma_dist_phi <- matrix(0, nrow = n_T*length(b_dist), ncol = iter)
est_sigma_dist_aic_phi <- matrix(0, nrow = n_T*length(b_dist), ncol = iter)

rho_dist_phi <- matrix(0, nrow = n_T*length(b_dist), ncol = iter)
rho_dist_aic_phi <- matrix(0, nrow = n_T*length(b_dist), ncol = iter)


sample_mean_dist <- numeric(length = iter)

bm_sigma_dist <- matrix(0, nrow = length(b_dist), ncol = iter)
est_sigma_dist <- matrix(0, nrow = length(b_dist), ncol = iter)
rho_dist <- matrix(0,nrow = length(b_dist), ncol = iter)
alpha_dist <- matrix(0,nrow = length(b_dist), ncol = iter)

est_sigma_dist_aic <- matrix(0, nrow = length(b_dist), ncol = iter)
rho_dist_aic <- matrix(0,nrow = length(b_dist), ncol = iter)
alpha_dist_aic <- matrix(0,nrow = length(b_dist), ncol = iter)


T_dist = c(1e3,1e4,1e5)
# b_dist[2] <- floor(T^(1/2))
# b_dist[1] <- floor(T)
sd = 0.005
# sd_dist = c(0.0060, 0.0045, 0.0030)
# This is to the code for which model has to be made for an AR1 process
y <- titanic[,1]
X <- as.matrix(titanic[, -1])

for(T in T_dist){
    print(paste("For Sample Size : ",T))
    b_dist[2] <- floor(T^(1/2))
    b_dist[1] <- floor(T^(1/3))
    
    for(it in 1:iter){
        print(paste("Iter : ", it))

        x <- bayes_logit_mh(y = y, X = X, N = T, prop.sd = sd) 
        
        x <- x[,6]

        n <- T

        for(b in b_dist){
            a <- floor(n/b)
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
            
            est_sigma <- (b*alpha_sq_c)/((1-rho_c)^2)

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

            est_sigma_aic <- (b*alpha_sq_c_aic)/((1-rho_c_aic)^2)

            est_sigma_dist_aic[which(b == b_dist), it] <- est_sigma_aic
        }
    }

    for(i in 1:length(b_dist)){
        bm_sigma_dist_phi[(which(T_dist == T)- 1)*length(b_dist) + i , ] <- bm_sigma_dist[i,]
        est_sigma_dist_phi[(which(T_dist == T)- 1)*length(b_dist) + i, ] <- est_sigma_dist[i,]
        est_sigma_dist_aic_phi[(which(T_dist == T)- 1)*length(b_dist) + i,]  <- est_sigma_dist_aic[i,]

        rho_dist_phi[(which(T_dist == T)- 1)*length(b_dist) + i, ] <- rho_dist[i,]
        rho_dist_aic_phi[(which(T_dist == T)- 1)*length(b_dist) + i, ] <- rho_dist_aic[i,]
    }
}
 
iter <- 100

for (T in 1:length(T_dist)){
    for(i in 1:length(b_dist)){

        bm_sigma_dist[i,] <- bm_sigma_dist_phi[(T- 1)*length(b_dist) + i , ] 
        est_sigma_dist[i,] <- est_sigma_dist_phi[(T- 1)*length(b_dist) + i, ] 
        est_sigma_dist_aic[i,] <- est_sigma_dist_aic_phi[(T- 1)*length(b_dist) + i,]  

        rho_dist[i,] <- rho_dist_phi[(T- 1)*length(b_dist) + i, ] 
        rho_dist_aic[i,] <- rho_dist_aic_phi[(T- 1)*length(b_dist) + i, ] 
        
    }

    mean_var = cbind(apply(bm_sigma_dist,1,mean),apply(est_sigma_dist,1,mean),apply(est_sigma_dist_aic,1,mean))
    var_var = cbind(apply(bm_sigma_dist,1,var),apply(est_sigma_dist,1,var),apply(est_sigma_dist_aic,1,var))

    data <- data.frame("mean"=mean_var, "var"=var_var)
    
    rownames(data) = c(b_dist[1],b_dist[2])
    
    write.csv(data, file = paste0("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/StatisticCSVFIle/T_",as.character(T_dist[T]),".csv"))

}

pdf(file = paste("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/Report_Images/T_SD_BLR.pdf",sep=""), height = 8, width = 8, pointsize = 10)


par(mfrow=c(3,2))   

#choosing the scale of the graph 
min_lim_glo = 0
max_lim_glo = 0

for(x in seq.int(0,5,2)){
    T <- T_dist[x/2 + 1]
    max_lim_glo = -1e4
    for(i in 1:length(b_dist)){
        max_lim_glo = max(max_lim_glo, max(bm_sigma_dist_phi[x+i,],est_sigma_dist_phi[x+i,],est_sigma_dist_aic_phi[x+i,]))
    }
}

# max_lim_glo = 0.003

print(paste0("Max Lim Glo : ", max_lim_glo))

for(x in seq.int(0,5,2)){

    T <- T_dist[x/2 + 1]
    b_dist[2] <- floor(T^(1/2))
    b_dist[1] <- floor(T^(1/3))

    for(i in 1:length(b_dist)) {
        bm_sigma_dist[i,] = bm_sigma_dist_phi[x+i,];
        est_sigma_dist[i,] = est_sigma_dist_phi[x+i,];
        est_sigma_dist_aic[i,] = est_sigma_dist_aic_phi[x+i,];

        # plot(density(bm_sigma_dist[i,]), xlim = c(min_lim - 0.5*(range), max_lim + 0.5*(range)), main = paste("T : ", T_dist[x/2 + 1],"Batch Size : ",b_dist[i]), lwd = 1, col = "red")
        plot(density(bm_sigma_dist[i,]), xlim = c(0, max_lim_glo), main = paste("T : ", T_dist[x/2 + 1], "Batch Size : ",b_dist[i] ), lwd = 1, col = "red")
        lines(density(est_sigma_dist[i,]), lwd = 1, col = "blue")
        lines(density(est_sigma_dist_aic[i,]), lwd = 1, col ="green")
        legend("topright", legend=c("BM ","Estimated ", "Estimated AIC "), col=c("red","blue", "green"), lty=1:1:1, cex=0.8)
    }
}

dev.off()

# pdf(file = paste("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/Report_Images/RHO_T_BLR.pdf",sep=""), height = 10, width = 8, pointsize = 10)

# par(mfrow=c(3,2))    

# for( x in seq.int(0,5,2)){
#     for(i in 1:length(b_dist)){
#         rho_dist[i, ] = rho_dist_phi[x+i, ]
#         rho_dist_aic[i, ] = rho_dist_aic_phi[x+i, ]
#         min_x = min(min(density(rho_dist[i, ])$x), max(density(rho_dist_aic[i, ])$x))
#         max_x = max(max(density(rho_dist[i, ])$x), max(density(rho_dist_aic[i, ])$x))
#         # plot(density(rho_dist_aic[i,]), col = "blue")
#         plot(density(rho_dist_aic[i,]), ylim = c(0,max(max(density(rho_dist[i, ])$y), max(density(rho_dist_aic[i, ])$y))), xlim = c(min_x, max_x),  col = "blue", main = paste("BS : ", b_dist[i])) # "Rho ", mean(rho_dist[i,]) ,"AIC Rho ", mean(rho_dist_aic[i,])
#         lines(density(rho_dist[i,]), lwd = 1, col = "red")
#         legend("topright", legend=c("AIC = FALSE ","AIC = TRUE"), col=c("red","blue"), lty=1:1, cex=0.8)
#     }
# }

# dev.off()


# pdf(file = paste("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Latec Images/Week 8/RHO_AIC.pdf",sep=""), height = 8, width = 10, pointsize = 10)
# par(mfrow=c(2,1))    

# for(i in 1:length(b_dist)) {
#     plot(density(rho_dist_aic[i,]), col = "red", main = paste("Rho ", mean(rho_dist_aic[i,]),"BS : ", b_dist[i]))
#     abline(v = mean(rho_dist_aic[i,]))
# }

# dev.off()

vals <- matrix(0, nrow = 6, ncol = 6)
# vals[2,] <- as.numeric(data[2,])
T_val = c(1000,10000,1e+05)


for(x in 0:(length(T_val)-1)){
    data <- read.csv(paste0("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/StatisticCSVFIle/T_",as.character(T_val[x+1]),".csv"))
    for(y in 1:2){
        vals[(2*x) + y,] <- as.numeric(data[y,])
    }
}

pdf(file = paste("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/Report_Images/BLR_Stat_T_Mean.pdf",sep=""), height = 8, width = 8, pointsize = 10)

name <- c("Mean","Variance");

par(oma = c(4,1,1,1), mfrow = c(1, 3), mar = c(2, 2, 1, 1))

mean_max = max(vals[,1:3])
var_max = max(vals[,4:6])
stat_max = c(mean_max,var_max)

for(j in 0:0){
    for(i in seq.int(0,5,2)){
        curr_data <- t(vals[(i+1):(i+2),])
        b_dist[2] <- floor(T_val[(i/2)+1]^(1/2))
        b_dist[1] <- floor(T_val[(i/2)+1]^(1/3))
    
        plot_data <- curr_data[((j*3)+1):((j*3)+3),]
        barplot(plot_data, main = paste0(name[j+1], " T :",T_val[(i/2)+1]), ylim = c(0,stat_max[j+1]), col = c("red","blue","green"), names.arg = c(b_dist[1],b_dist[2]),beside = TRUE)
    }
}


for(j in 0:0){
    for(i in seq.int(0,5,2)){
        curr_data <- t(vals[(i+1):(i+2),])
        b_dist[2] <- floor(T_val[(i/2)+1]^(1/2))
        b_dist[1] <- floor(T_val[(i/2)+1]^(1/3))
    
        plot_data <- curr_data[((j*3)+1):((j*3)+3),]
        barplot(plot_data, main = paste0(name[j+1], " T :",T_val[(i/2)+1]), ylim = c(0,stat_max[j+1]), col = c("red","blue","green"), names.arg = c(b_dist[1],b_dist[2]),beside = TRUE)
    }
}

# par(oma = c(4,1,1,1), mfrow = c(1, 3), mar = c(2, 2, 1, 1))
# for(i in seq.int(0,5,2)){
#     curr_data <- t(vals[(i+1):(i+2),])
#     b_dist[2] <- floor(T_val[(i/2)+1]^(1/2))
#     b_dist[1] <- floor(T_val[(i/2)+1]^(1/3))

#     for(j in 0:0){
#         plot_data <- curr_data[((j*3)+1):((j*3)+3),]
#         barplot(plot_data, main = paste0(name[j+1], " T :",T_val[(i/2)+1]), ylim = range(pretty(c(0,plot_data))), col = c("red","blue","green"), names.arg = c(b_dist[1],b_dist[2]),beside = TRUE)
#     }
# }

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(1, 1, 1, 1), new = TRUE)
plot(10, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c("red","blue", "green")
legend(x = "bottom",inset = 0,
        legend = c("BM ","Estimated ", "Estimated AIC "), 
        col=plot_colors, lwd=5, cex=.8, horiz = TRUE)

dev.off()

curr_data <- vals[1:2,]
t_curr_data <- t(curr_data)
# reshape_data <- reshape(curr_data, )
broken_data <- t_curr_data[1:3,]

barplot(broken_data, main = "Mean", ylim = range(pretty(c(0,broken_data))), col = c("red","blue","green"), names.arg = c(b_dist[1],b_dist[2]),beside = TRUE)
legend("topright", legend=c("BM ","Estimated ", "Estimated AIC "), col=c("red","blue", "green"), lty=1:1:1, cex=0.75)



par(mfrow=c(1,1))

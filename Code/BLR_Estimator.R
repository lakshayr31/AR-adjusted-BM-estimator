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
chain <- bayes_logit_mh(y = y, X = X, N = 1e3, prop.sd = .35) 


b_dist <- numeric(2);
iter <- 100
n_sd = 3
#for plotting everything on the same graph
bm_sigma_dist_phi <- matrix(0, nrow = n_sd*length(b_dist), ncol = iter)
est_sigma_dist_phi <- matrix(0, nrow = n_sd*length(b_dist), ncol = iter)
est_sigma_dist_aic_phi <- matrix(0, nrow = n_sd*length(b_dist), ncol = iter)

rho_dist_phi <- matrix(0, nrow = n_sd*length(b_dist), ncol = iter)
rho_dist_aic_phi <- matrix(0, nrow = n_sd*length(b_dist), ncol = iter)


sample_mean_dist <- numeric(length = iter)

bm_sigma_dist <- matrix(0, nrow = length(b_dist), ncol = iter)
est_sigma_dist <- matrix(0, nrow = length(b_dist), ncol = iter)
rho_dist <- matrix(0,nrow = length(b_dist), ncol = iter)
alpha_dist <- matrix(0,nrow = length(b_dist), ncol = iter)

est_sigma_dist_aic <- matrix(0, nrow = length(b_dist), ncol = iter)
rho_dist_aic <- matrix(0,nrow = length(b_dist), ncol = iter)
alpha_dist_aic <- matrix(0,nrow = length(b_dist), ncol = iter)


T <- 1e4
b_dist[2] <- floor(T^(1/2))
b_dist[1] <- b_dist[2]/2

sd_dist = c(0.0060, 0.0045, 0.0030)
# This is to the code for which model has to be made for an AR1 process
for(sd in 1:length(sd_dist)){
    for(it in 1:iter){
        print(paste("Iter : ", it))

        y <- titanic[,1]
        X <- as.matrix(titanic[, -1])

        ## acceptance is too low! we want 23%
        ## so decrease proposal variance
        x <- bayes_logit_mh(y = y, X = X, N = T, prop.sd = sd_dist[sd]) 
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
        bm_sigma_dist_phi[(sd- 1)*length(b_dist) + i , ] <- bm_sigma_dist[i,]
        est_sigma_dist_phi[(sd- 1)*length(b_dist) + i, ] <- est_sigma_dist[i,]
        est_sigma_dist_aic_phi[(sd- 1)*length(b_dist) + i,]  <- est_sigma_dist_aic[i,]

        rho_dist_phi[(sd- 1)*length(b_dist) + i, ] <- rho_dist[i,]
        rho_dist_aic_phi[(sd- 1)*length(b_dist) + i, ] <- rho_dist_aic[i,]
    }
}
 
iter <- 100

for (sd in 1:length(sd_dist)){
    for(i in 1:length(b_dist)){
        
        y <- titanic[,1]
        X <- as.matrix(titanic[, -1])

        ## acceptance is too low! we want 23%
        ## so decrease proposal variance
        var_x <- numeric(iter) 
        for(it in 1:iter){
            print(paste("Iter : ", it))

            x <- bayes_logit_mh(y = y, X = X, N = 1e4, prop.sd = sd_dist[sd]) 
            x <- x[,6]
            sigma_true <- var(x)
            var_x[it] <- sigma_true
        }
        sigma_true <- var(var_x)

        bm_sigma_dist[i,] <- bm_sigma_dist_phi[(sd- 1)*length(b_dist) + i , ] 
        est_sigma_dist[i,] <- est_sigma_dist_phi[(sd- 1)*length(b_dist) + i, ] 
        est_sigma_dist_aic[i,] <- est_sigma_dist_aic_phi[(sd- 1)*length(b_dist) + i,]  

        rho_dist[i,] <- rho_dist_phi[(sd- 1)*length(b_dist) + i, ] 
        rho_dist_aic[i,] <- rho_dist_aic_phi[(sd- 1)*length(b_dist) + i, ] 
        cbind(apply(bm_sigma_dist,1,mean),apply(est_sigma_dist,1,mean),apply(est_sigma_dist_aic,1,mean))
        cbind(apply(bm_sigma_dist,1,var),apply(est_sigma_dist,1,var),apply(est_sigma_dist_aic,1,var)) 

        mean_var = cbind(apply(bm_sigma_dist,1,mean),apply(est_sigma_dist,1,mean),apply(est_sigma_dist_aic,1,mean))
        var_var = cbind(apply(bm_sigma_dist,1,var),apply(est_sigma_dist,1,var),apply(est_sigma_dist_aic,1,var)) 

        mse_var <- matrix(0, nrow = length(b_dist), ncol = 3)
    }
    for(i in 1:length(b_dist)){
            mse_bm = (1/n)*sum((sigma_true-bm_sigma_dist[i,])^2);
            mse_est = (1/n)*sum((sigma_true-est_sigma_dist[i,])^2);
            mse_est_aic = (1/n)*sum((sigma_true-est_sigma_dist_aic[i,])^2);
            # mse_var[i,] = formatC(c(mse_bm, mse_est, mse_est_aic), format="e", digits = 4)
            mse_var[i,] = c(mse_bm, mse_est, mse_est_aic);
            print(paste("Batch Size ",b_dist[i]," MSE BM: ",mse_bm," MSE EST: ",mse_est,"MSE AIC EST: ", mse_est_aic));
        }

    data <- data.frame("mean"=mean_var, "var"=var_var, "mse"=mse_var)
    rownames(data) = c(b_dist[1],b_dist[2])
    write.csv(data, file = paste("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/StatisticCSVFIle/sd_",sd_dist[sd]))
}
 

pdf(file = paste("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/Report_Images/SD_AR_BLR.pdf",sep=""), height = 10, width = 8, pointsize = 10)


par(mfrow=c(4,2))   


for(x in seq.int(0,5,2)){
    for(i in 1:length(b_dist)) {
        bm_sigma_dist[i,] = bm_sigma_dist_phi[x+i,];
        est_sigma_dist[i,] = est_sigma_dist_phi[x+i,];
        est_sigma_dist_aic[i,] = est_sigma_dist_aic_phi[x+i,];
        min_lim = min(bm_sigma_dist[i,],est_sigma_dist[i,],est_sigma_dist_aic[i,])
        max_lim = max(bm_sigma_dist[i,],est_sigma_dist[i,],est_sigma_dist_aic[i,])
        range = max_lim - min_lim
        plot(density(bm_sigma_dist[i,]), xlim = c(min_lim - 0.5*(range), max_lim + 0.5*(range)), main = paste("SD : ", sd_dist[x/2 + 1],"Batch Size : ",b_dist[i]), lwd = 1, col = "red")
        abline(v = sigma_true[x/2 + 1])
        lines(density(est_sigma_dist[i,]), lwd = 1, col = "blue")
        lines(density(est_sigma_dist_aic[i,]), lwd = 1, col ="green")
        legend("topright", legend=c("BM ","Estimated ", "Estimated AIC "), col=c("red","blue", "green"), lty=1:1:1, cex=0.8)
    }
}

dev.off()

pdf(file = paste("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/Report_Images/Rho_BLR.pdf",sep=""), height = 10, width = 8, pointsize = 10)

par(mfrow=c(4,2))    
for( x in seq.int(0,5,2)){
    for(i in 1:length(b_dist)){
        rho_dist[i, ] = rho_dist_phi[x+i, ]
        rho_dist_aic[i, ] = rho_dist_aic_phi[x+i, ]
        min_x = min(min(density(rho_dist[i, ])$x), max(density(rho_dist_aic[i, ])$x))
        max_x = max(max(density(rho_dist[i, ])$x), max(density(rho_dist_aic[i, ])$x))
        plot(density(rho_dist_aic[i,]), col = "blue")
        # plot(density(rho_dist_aic[i,]), ylim = c(0,max(max(density(rho_dist[i, ])$y), max(density(rho_dist_aic[i, ])$y))), xlim = c(min_x, max_x),  col = "blue", main = paste("BS : ", b_dist[i])) # "Rho ", mean(rho_dist[i,]) ,"AIC Rho ", mean(rho_dist_aic[i,])
        # lines(density(rho_dist[i,]), lwd = 1, col = "red")
        legend("topright", legend=c("AIC = FALSE ","AIC = TRUE"), col=c("red","blue"), lty=1:1, cex=0.8)
    }
}

dev.off()


pdf(file = paste("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Latec Images/Week 8/RHO_AIC.pdf",sep=""), height = 8, width = 10, pointsize = 10)
par(mfrow=c(2,1))    

for(i in 1:length(b_dist)) {
    plot(density(rho_dist_aic[i,]), col = "red", main = paste("Rho ", mean(rho_dist_aic[i,]),"BS : ", b_dist[i]))
    abline(v = mean(rho_dist_aic[i,]))
}

dev.off()


par(mfrow=c(1,1))


# When I use this function then it doesnt give the right answer
calculate_est <- function(b, alpha_sq, rho) {
    return ((b*alpha_sq)/((1-rho)^2))
}

# b_dist <- c(100,200,500,1000);
# b_dist <- c(10,20,50,100);

b_dist <- numeric(2);
iter <- 100
n_phi = 4
#for plotting everything on the same graph
bm_sigma_dist_phi <- matrix(0, nrow = n_phi*length(b_dist), ncol = iter)
est_sigma_dist_phi <- matrix(0, nrow = n_phi*length(b_dist), ncol = iter)
est_sigma_dist_aic_phi <- matrix(0, nrow = n_phi*length(b_dist), ncol = iter)

rho_dist_phi <- matrix(0, nrow = n_phi*length(b_dist), ncol = iter)
rho_dist_aic_phi <- matrix(0, nrow = n_phi*length(b_dist), ncol = iter)


sample_mean_dist <- numeric(length = iter)

bm_sigma_dist <- matrix(0, nrow = length(b_dist), ncol = iter)

est_sigma_dist <- matrix(0, nrow = length(b_dist), ncol = iter)
rho_dist <- matrix(0,nrow = length(b_dist), ncol = iter)
alpha_dist <- matrix(0,nrow = length(b_dist), ncol = iter)

est_sigma_dist_aic <- matrix(0, nrow = length(b_dist), ncol = iter)
rho_dist_aic <- matrix(0,nrow = length(b_dist), ncol = iter)
alpha_dist_aic <- matrix(0,nrow = length(b_dist), ncol = iter)


T <- 1e5
b_dist[2] <- floor(T^(1/2))
b_dist[1] <- floor(T^(1/3))

phi <- 0.99
alpha <- 1


phi_dist = c(0.6,0.9,0.95,0.99)

sigma_true <- numeric(4)
for(i in 1:length(phi_dist)){
    sigma_true[i] <- (alpha^2)/((1-phi_dist[i])^2)
}


# This is to the code for which model has to be made for an AR1 process
for(it in 1:iter){
    print(paste("Iter : ", it))

    x <- numeric(T)
    x[1] <- rnorm(1,0,1)

    for(i in 2:T){
        x[i] <- x[i-1]*phi + rnorm(1,0,alpha)
    }

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
        
        est_sigma <- calculate_est(b, alpha_sq_c, rho_c)

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

        est_sigma_aic <- calculate_est(b, alpha_sq_c_aic, rho_c_aic)

        est_sigma_dist_aic[which(b == b_dist), it] <- est_sigma_aic
    }
}

cbind(apply(bm_sigma_dist,1,mean),apply(est_sigma_dist,1,mean),apply(est_sigma_dist_aic,1,mean))
cbind(apply(bm_sigma_dist,1,var),apply(est_sigma_dist,1,var),apply(est_sigma_dist_aic,1,var)) 

mean_var = cbind(apply(bm_sigma_dist,1,mean),apply(est_sigma_dist,1,mean),apply(est_sigma_dist_aic,1,mean))
var_var = cbind(apply(bm_sigma_dist,1,var),apply(est_sigma_dist,1,var),apply(est_sigma_dist_aic,1,var)) 

mse_var <- matrix(0, nrow = length(b_dist), ncol = 3)
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
data
write.csv(data, file = paste("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/StatisticCSVFIle/phi_",phi))
# 0.6, 0.9. 0.95, 0.99

x = 6
for(i in 1:length(b_dist)){
    print(x + i)
    bm_sigma_dist_phi[x + i, ] <- bm_sigma_dist[i,]
    est_sigma_dist_phi[x + i, ] <- est_sigma_dist[i,]
    est_sigma_dist_aic_phi[x + i,]  <- est_sigma_dist_aic[i,]

    rho_dist_phi[x + i, ] <- rho_dist[i,]
    rho_dist_aic_phi[x + i, ] <- rho_dist_aic[i,]
}


pdf(file = paste("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/Report_Images/Phi_AR_1.pdf",sep=""), height = 10, width = 8, pointsize = 10)


par(mfrow=c(4,2))   


for(x in seq.int(0,7,2)){
    for(i in 1:length(b_dist)) {
        bm_sigma_dist[i,] = bm_sigma_dist_phi[x+i,];
        est_sigma_dist[i,] = est_sigma_dist_phi[x+i,];
        est_sigma_dist_aic[i,] = est_sigma_dist_aic_phi[x+i,];
        min_lim = min(bm_sigma_dist[i,],est_sigma_dist[i,],est_sigma_dist_aic[i,])
        max_lim = max(bm_sigma_dist[i,],est_sigma_dist[i,],est_sigma_dist_aic[i,])
        range = max_lim - min_lim
        plot(density(bm_sigma_dist[i,]), xlim = c(min_lim - 0.5*(range), max_lim + 0.5*(range)), main = paste("Phi : ", phi_dist[x/2 + 1],"Batch Size : ",b_dist[i]), lwd = 1, col = "red")
        abline(v = sigma_true[x/2 + 1])
        lines(density(est_sigma_dist[i,]), lwd = 1, col = "blue")
        lines(density(est_sigma_dist_aic[i,]), lwd = 1, col ="green")
        legend("topright", legend=c("BM ","Estimated ", "Estimated AIC "), col=c("red","blue", "green"), lty=1:1:1, cex=0.8)
    }
}

dev.off()

pdf(file = paste("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/Report_Images/Rho_2.pdf",sep=""), height = 10, width = 8, pointsize = 10)

par(mfrow=c(4,2))    
for( x in seq.int(0,7,2)){
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


pdf(file = paste("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Latec Images/Week 8/",rho,"_RHO_AIC.pdf",sep=""), height = 8, width = 10, pointsize = 10)
par(mfrow=c(2,1))    

for(i in 1:length(b_dist)) {
    plot(density(rho_dist_aic[i,]), col = "red", main = paste("Rho ", mean(rho_dist_aic[i,]),"BS : ", b_dist[i]))
    abline(v = mean(rho_dist_aic[i,]))
}

dev.off()


par(mfrow=c(1,1))

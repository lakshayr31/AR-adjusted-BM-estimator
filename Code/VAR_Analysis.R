library("mvtnorm")
library("vars")
library("fastmatrix")
library("matrixcalc")
library("pracma")

#### setting some things upfront

p <- 5
T <- 1e5
alpha <- 0.4

W <- diag(alpha,p)

#### First I need to generate a series of samples from a VAR (1)

#### Generating a phi
sigma <- diag(p)

for(i in 1:p){
    if(i < p/2){
        sigma[i,i] <- runif(1,0,0.5)
    }
    else {
        sigma[i,i] <- runif(1,0.5,1)
    }
}

u <- randortho(p, "orthonormal")

phi <- u %*% sigma %*% t(u)
sigma
eigen(phi)$values

iter <- 100
b_dist <- numeric(2);

det_sigma_bm <- matrix(0, nrow = length(b_dist), ncol = iter)
det_sigma_est <- matrix(0, nrow = length(b_dist), ncol = iter)
error_bm <- matrix(0, nrow = length(b_dist), ncol = iter)
error_est <- matrix(0, nrow = length(b_dist), ncol = iter)


b_dist[1] <- floor(T^(1/3))
b_dist[2] <- floor(T^(1/2))

V <- solve(diag(p^2) - kronecker.prod(phi,phi)) %*% c(diag(alpha,p))
V_mat <- matrix(V,p,p)

one_mat <- diag(1,p)
act_sigma <- solve(one_mat-phi)%*%V_mat + V_mat%*%solve(one_mat-t(phi)) - V_mat    
det_act_sigma <-  determinant(act_sigma, logarithm=TRUE)$modulus[1]

for(x in 1:iter){
    print(paste0("Iter : ", x))
    X0 <- rnorm(p)

    X <- matrix(0, nrow = p, ncol = T)

    X[,1] <- X0

    for(i in 2:T){
        X[,i] <- phi %*% X[,i-1] + t(rmvnorm(1,numeric(p),W)) 
    }

    for(b in b_dist){

        a <- floor(T/b)
        Y <- matrix(0,nrow = p, ncol = a)

        for(i in 0:(a-1)){
            Y[,i+1] <- rowMeans(X[,((b*i)+1):(b*(i+1))])
        }

        theta_hat <- rowMeans(X)

        temp_sum <- matrix(0,p,p)

        for(i in 1:a){
            temp_sum <- temp_sum + (Y[,i] - theta_hat) %*% t(Y[,i] - theta_hat)
        }

        bm_cov_matrix <- (b/(a-1))*(temp_sum)

        bm_error <- frobenius.norm(bm_cov_matrix-act_sigma)/frobenius.norm(act_sigma)

        print(paste0("Batch Means Error ",bm_error))

        data <- t(Y)

        model <- VAR(y = data, p = 1,lag.max = NULL)

        W_est <- summary(model)$covres

        coeff <- Acoef(model)[[1]]

        phi_est <- coeff

        V_est <- solve(diag(p^2) - kronecker.prod(phi_est,phi_est)) %*% c(W_est)

        V_est_mat <- matrix(V_est,p,p)

        one_mat <- diag(1,p)

        est_sigma <- b*(solve(one_mat-phi_est)%*%V_est_mat + V_est_mat %*% solve(one_mat-t(phi_est)) - V_est_mat)

        est_error <- frobenius.norm(est_sigma-act_sigma)/frobenius.norm(act_sigma)

        print(paste0("Estimated Error ",est_error))

        error_bm[which(b_dist == b), x] <- bm_error
        error_est[which(b_dist == b),x] <- est_error

        det_sigma_bm[which(b_dist == b),x] <- determinant(bm_cov_matrix, logarithm=TRUE)$modulus[1]
        det_sigma_est[which(b_dist == b),x] <- determinant(est_sigma, logarithm=TRUE)$modulus[1]
        
    }
}


mean_var = cbind(apply(det_sigma_bm,1,mean),apply(det_sigma_est,1,mean))
var_var = cbind(apply(det_sigma_bm,1,var),apply(det_sigma_est,1,var))

data <- data.frame("mean"=mean_var, "var"=var_var)

rownames(data) = c(b_dist[1],b_dist[2])

write.csv(data, file = paste0("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/StatisticCSVFIle/VAR_GEN.csv"))

# Plotting the data 

vals <- matrix(0, nrow = 2, ncol = 4)
# vals[2,] <- as.numeric(data[2,])

data <- read.csv(paste0("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/StatisticCSVFIle/VAR_GEN.csv"))
vals <- data

pdf(file = paste("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/Report_Images/VAR_GEN_MeanAndVar.pdf",sep=""), height = 4, width = 8, pointsize = 10)

name <- c("Mean","Variance");

par(oma = c(4,1,1,1), mfrow = c(1, 2), mar = c(2, 2, 1, 1))

mean_max = min(vals[,1:2])
var_max = max(vals[,3:4])
stat_max = c(mean_max,var_max)

# for(j in 0:0){
    # for(i in seq.int(0,5,2)){
i <- 0
curr_data <- t(vals[(i+1):(i+2),])
b_dist[2] <- floor(T_val[(i/2)+1]^(1/2))
b_dist[1] <- floor(T_val[(i/2)+1]^(1/3))

j <- 0
plot_data <- curr_data[((j*2)+1):((j*2)+2),]
barplot(plot_data, main = paste0(name[j+1], " T :", 1e5), ylim = c(0,stat_max[j+1]), col = c("red","blue"), names.arg = c(b_dist[1],b_dist[2]),beside = TRUE)

j <- 1
plot_data <- curr_data[((j*2)+1):((j*2)+2),]
barplot(plot_data, main = paste0(name[j+1], " T :",1e5), ylim = c(0,stat_max[j+1]), col = c("red","blue"), names.arg = c(b_dist[1],b_dist[2]),beside = TRUE)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(1, 1, 1, 1), new = TRUE)
plot(10, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c("red","blue")
legend(x = "bottom",inset = 0,
        legend = c("BM ","Estimated "), 
        col=plot_colors, lwd=5, cex=.8, horiz = TRUE)

dev.off()

# pdf(file = paste("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/Report_Images/BLR_VAR_T_Variance.pdf",sep=""), height = 8, width = 8, pointsize = 10)

# for(j in 1:1){
#     for(i in seq.int(0,5,2)){
#         curr_data <- t(vals[(i+1):(i+2),])
#         b_dist[2] <- floor(T_val[(i/2)+1]^(1/2))
#         b_dist[1] <- floor(T_val[(i/2)+1]^(1/3))
    
#         plot_data <- curr_data[((j*2)+1):((j*2)+2),]
#         barplot(plot_data, main = paste0(name[j+1], " T :",T_val[(i/2)+1]), ylim = c(0,stat_max[j+1]), col = c("red","blue"), names.arg = c(b_dist[1],b_dist[2]),beside = TRUE)
#     }
# }

# dev.off()


pdf(file = paste("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/Report_Images/VAR_GEN_DET.pdf",sep=""), height = 5, width = 10, pointsize = 10)

par(mfrow = c(1, 2))

for(b in b_dist){
    max_y <- max(max_y,(c(density(det_sigma_bm[ which(b == b_dist), ])$y, density(det_sigma_est[which(b == b_dist), ])$y)))
    max_x <- max(max_x,(c(density(det_sigma_bm[ which(b == b_dist), ])$x, density(det_sigma_est[which(b == b_dist), ])$x)))
    # plot(density(det_sigma_est[ which(b == b_dist), ]), main = paste0(" Determinants T : ", T, " Batch Size : ", b) ,lwd = 1, col ="red")
    # plot(density(det_sigma_est[ which(b == b_dist), ]), xlim = c(0,max_x*1.20) , ylim = c(0,max_y*1.10) , main = paste0(" Determinants T : ", T, " Batch Size : ", b) ,lwd = 1, col ="red")
    # lines(density(det_sigma_bm[ which(b == b_dist), ]), lwd = 1, col = "blue")
    # abline(v = det_act_sigma)
    # legend("topright", legend=c("BM ", "EST "), col=c("blue","red"), lty=1:1:1, cex=0.8)

}

for(b in b_dist){
    # max_y <- max(c(density(det_sigma_bm[ which(b == b_dist), ])$y, density(det_sigma_est[which(b == b_dist), ])$y))
    # max_x <- max(c(density(det_sigma_bm[ which(b == b_dist), ])$x, density(det_sigma_est[which(b == b_dist), ])$x))
    # plot(density(det_sigma_est[ which(b == b_dist), ]), main = paste0(" Determinants T : ", T, " Batch Size : ", b) ,lwd = 1, col ="red")
    plot(density(det_sigma_est[ which(b == b_dist), ]), xlim = c(0,max_x*1.20) , ylim = c(0,max_y*1.10) , main = paste0(" Determinants T : ", T, " Batch Size : ", b) ,lwd = 1, col ="red")
    lines(density(det_sigma_bm[ which(b == b_dist), ]), lwd = 1, col = "blue")
    abline(v = det_act_sigma)
    # legend("topright", legend=c("BM ", "EST "), col=c("blue","red"), lty=1:1:1, cex=0.8)

}

dev.off()


pdf(file = paste("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/Report_Images/VAR_GEN_ERR.pdf",sep=""), height = 5, width = 10, pointsize = 10)


par(mfrow = c(1, 2))
max_x <- 0
max_y <- 0
for(b in b_dist){
    max_y <- max(max_y,c(density(error_bm[ which(b == b_dist), ])$y, density(error_est[which(b == b_dist), ])$y))
    max_x <- max(max_x,c(density(error_bm[ which(b == b_dist), ])$x, density(error_est[which(b == b_dist), ])$x))
    # plot(density(error_est[ which(b == b_dist), ]), xlim = c(0,max_x*1.50) , ylim = c(0,max_y*1.10) , main = paste0(" Error T : ", T, " Batch Size : ", b) ,lwd = 1, col ="red")
    # lines(density(error_bm[ which(b == b_dist), ]), lwd = 1, col = "blue")
    # legend("topright", legend=c("BM ", "EST "), col=c("blue","red"), lty=1:1:1, cex=0.8)
}

for(b in b_dist){
    # max_y <- max(max_y,c(density(error_bm[ which(b == b_dist), ])$y, density(error_est[which(b == b_dist), ])$y))
    # max_x <- max(max_x,c(density(error_bm[ which(b == b_dist), ])$x, density(error_est[which(b == b_dist), ])$x))
    plot(density(error_est[ which(b == b_dist), ]), xlim = c(0,max_x*1.50) , ylim = c(0,max_y*1.10) , main = paste0(" Error T : ", T, " Batch Size : ", b) ,lwd = 1, col ="red")
    lines(density(error_bm[ which(b == b_dist), ]), lwd = 1, col = "blue")
    legend("topright", legend=c("BM ", "EST "), col=c("blue","red"), lty=1:1:1, cex=0.8)
}

# A lot of questions in this code 
# 1 Does the Phi matrix need to be similar to a covaraince matrix?
# 2 I am not getting similar values, and I can't figure out what am I doing wrong?
# 3 Is frobenius norm the onlu way to measure out performance, no other comparision or evalutaion metric available?


dev.off()

library("mvtnorm")
library("vars")
library("fastmatrix")
library("matrixcalc")
library("pracma")

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

# titanic <- read.csv("https://dvats.github.io/assets/titanic.csv")

y <- titanic[,1]
X <- as.matrix(titanic[, -1])
p <- dim(X)[2]
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

bm_sigma_dist_T <- matrix(0, nrow = n_T*length(b_dist), ncol = iter)
est_sigma_dist_T <- matrix(0, nrow = n_T*length(b_dist), ncol = iter)

bm_sigma_dist <- matrix(0, nrow = length(b_dist), ncol = iter)
est_sigma_dist <- matrix(0, nrow = length(b_dist), ncol = iter)

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
        x <- t(x)
        n <- T

        for(b in b_dist){
            a <- floor(n/b)
            Y <- matrix(0,nrow = p, ncol = a)
        
            for(i in 0:(a-1)){
                Y[,i+1] <- rowMeans(x[,((b*i)+1):(b*(i+1))])
            }

            theta_hat <- rowMeans(x)
            temp_sum <- matrix(0,p,p)

            for(i in 1:a){
                temp_sum <- temp_sum + (Y[,i] - theta_hat) %*% t(Y[,i] - theta_hat)
            }

            bm_cov_matrix <- (b/(a-1))*(temp_sum)

            data <- t(Y)
            model <- VAR(y = data, p = 1,lag.max = NULL)

            W_est <- summary(model)$covres

            coeff <- Acoef(model)[[1]]

            phi_est <- coeff

            V_est <- solve(diag(p^2) - kronecker.prod(phi_est,phi_est)) %*% c(W_est)

            V_est_mat <- matrix(V_est,p,p)

            one_mat <- diag(1,p)

            est_sigma <- b*(solve(one_mat-phi_est)%*%V_est_mat + V_est_mat %*% solve(one_mat-t(phi_est)) - V_est_mat)

            bm_sigma_dist[which(b == b_dist), it] <- determinant(bm_cov_matrix, logarithm=TRUE)$modulus[1]
            est_sigma_dist[which(b == b_dist), it] <- determinant(est_sigma, logarithm=TRUE)$modulus[1]
        }
    }

    for(i in 1:length(b_dist)){
        bm_sigma_dist_T[(which(T_dist == T)- 1)*length(b_dist) + i , ] <- bm_sigma_dist[i,]
        est_sigma_dist_T[(which(T_dist == T)- 1)*length(b_dist) + i, ] <- est_sigma_dist[i,]
    }
}
 
iter <- 100

for (T in 1:length(T_dist)){
    for(i in 1:length(b_dist)){

        bm_sigma_dist[i,] <- bm_sigma_dist_T[(T- 1)*length(b_dist) + i , ] 
        est_sigma_dist[i,] <- est_sigma_dist_T[(T- 1)*length(b_dist) + i, ] 
    }

    mean_var = cbind(apply(bm_sigma_dist,1,mean),apply(est_sigma_dist,1,mean))
    var_var = cbind(apply(bm_sigma_dist,1,var),apply(est_sigma_dist,1,var))

    data <- data.frame("mean"=mean_var, "var"=var_var)
    
    # rownames(data) = c(b_dist[1],b_dist[2])
    
    write.csv(data, file = paste0("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/StatisticCSVFIle/BLR_VAR_T_",as.character(T_dist[T]),".csv"))

}

pdf(file = paste("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/Report_Images/T_SD_BLR_VAR.pdf",sep=""), height = 10, width = 8, pointsize = 10)

par(mfrow=c(3,2))   

#choosing the scale of the graph 
min_lim_glo = 1e4
max_lim_glo = -1e4

for(x in seq.int(0,5,2)){
    T <- T_dist[x/2 + 1]
    max_lim_glo = -1e4
    for(i in 1:length(b_dist)){
        max_lim_glo = max(max_lim_glo, max(bm_sigma_dist_T[x+i,],est_sigma_dist_T[x+i,]))
        min_lim_glo = min(min_lim_glo, min(bm_sigma_dist_T[x+i,],est_sigma_dist_T[x+i,]))
    }
}

print(paste0("Max Lim Glo : ", max_lim_glo))

print(paste0("Min Lim Glo : ", min_lim_glo))

for(x in seq.int(0,5,2)){

    T <- T_dist[x/2 + 1]
    b_dist[2] <- floor(T^(1/2))
    b_dist[1] <- floor(T^(1/3))
    
    for(i in 1:length(b_dist)) {
        bm_sigma_dist[i,] = bm_sigma_dist_T[x+i,];
        est_sigma_dist[i,] = est_sigma_dist_T[x+i,];
        # max_lim = max(c(bm_sigma_dist[i,],est_sigma_dist[i,]))
        # min_lim = min(c(bm_sigma_dist[i,],est_sigma_dist[i,]))
        # plot(density(bm_sigma_dist[i,]), xlim = c(min_lim - 0.5*(range), max_lim + 0.5*(range)), main = paste("T : ", T_dist[x/2 + 1],"Batch Size : ",b_dist[i]), lwd = 1, col = "red")
        plot(density(bm_sigma_dist[i,]), xlim = c(min_lim_glo, max_lim_glo), main = paste("T : ", T_dist[x/2 + 1], "Batch Size : ",b_dist[i] ), lwd = 1, col = "red")
        lines(density(est_sigma_dist[i,]), lwd = 1, col = "blue")
        legend("topright", legend=c("BM ","Estimated "), col=c("red","blue"), lty=1:1, cex=0.8)
    }
}

dev.off()

################################################################################

vals <- matrix(0, nrow = 6, ncol = 4)
# vals[2,] <- as.numeric(data[2,])
T_val = c(1e3,1e4,1e+05)


for(x in 0:(length(T_val)-1)){
    data <- read.csv(paste0("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/StatisticCSVFIle/BLR_VAR_T_",as.character(T_val[x+1]),".csv"))
    for(y in 1:2){
        vals[(2*x) + y,] <- as.numeric(data[y,])
    }
}

pdf(file = paste("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/Report_Images/BLR_VAR_T_Mean.pdf",sep=""), height = 6, width = 8, pointsize = 10)

name <- c("Mean","Variance");

par(oma = c(4,1,1,1), mfrow = c(1, 3), mar = c(2, 2, 1, 1))

mean_max = min(vals[,1:2])
var_max = max(vals[,3:4])
stat_max = c(mean_max,var_max)

for(j in 0:0){
    for(i in seq.int(0,5,2)){
        curr_data <- t(vals[(i+1):(i+2),])
        b_dist[2] <- floor(T_val[(i/2)+1]^(1/2))
        b_dist[1] <- floor(T_val[(i/2)+1]^(1/3))
    
        plot_data <- curr_data[((j*2)+1):((j*2)+2),]
        barplot(plot_data, main = paste0(name[j+1], " T :",T_val[(i/2)+1]), ylim = c(0,stat_max[j+1]), col = c("red","blue"), names.arg = c(b_dist[1],b_dist[2]),beside = TRUE)
        # legend("topright", legend=c("BM ", "EST "), col=c("red","blue"), lty=1:1:1, cex=0.8)
    }
}


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(1, 1, 1, 1), new = TRUE)
plot(10, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c("red","blue")
legend(x = "bottom",inset = 0,
        legend = c("BM ","Estimated "), 
        col=plot_colors, lwd=5, cex=.8, horiz = TRUE)


dev.off()

pdf(file = paste("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/Report_Images/BLR_VAR_T_Variance.pdf",sep=""), height = 6, width = 8, pointsize = 10)

par(oma = c(4,1,1,1), mfrow = c(1, 3), mar = c(2, 2, 1, 1))

for(j in 1:1){
    for(i in seq.int(0,5,2)){
        curr_data <- t(vals[(i+1):(i+2),])
        b_dist[2] <- floor(T_val[(i/2)+1]^(1/2))
        b_dist[1] <- floor(T_val[(i/2)+1]^(1/3))
    
        plot_data <- curr_data[((j*2)+1):((j*2)+2),]
        barplot(plot_data, main = paste0(name[j+1], " T :",T_val[(i/2)+1]), ylim = c(0,stat_max[j+1]), col = c("red","blue"), names.arg = c(b_dist[1],b_dist[2]),beside = TRUE)
        # legend("topright", legend=c("BM ", "EST "), col=c("red","blue"), lty=1:1:1, cex=0.8)

    }
}


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(1, 1, 1, 1), new = TRUE)
plot(10, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c("red","blue")
legend(x = "bottom",inset = 0,
        legend = c("BM ","Estimated "), 
        col=plot_colors, lwd=5, cex=.8, horiz = TRUE)


dev.off()

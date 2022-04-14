data <- read.csv("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/StatisticCSVFIle/sd_0_003.csv")
print(data)

vals <- matrix(0, nrow = 6, ncol = 9)
# vals[2,] <- as.numeric(data[2,])
sd_val = c(3,45,6)


for(x in 0:(length(sd_val)-1)){
    data <- read.csv(paste0("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/StatisticCSVFIle/sd_0_00",sd_val[x+1],".csv"))
    for(y in 1:2){
        vals[(2*x) + y,] <- as.numeric(data[y,])
    }
}

pdf(file = paste("/mnt/d/2022_Jan-Apr/Semester 8/Mth UGP/Code/Report_Images/BLR_Stats_2.pdf",sep=""), height = 10, width = 8, pointsize = 10)

name <- c("Mean","Variance","MSE");
par(oma = c(4,1,1,1), mfrow = c(3, 3), mar = c(2, 2, 1, 1))
for(i in seq.int(0,5,2)){
    curr_data <- t(vals[(i+1):(i+2),])
    for(j in 0:2){
        plot_data <- curr_data[((j*3)+1):((j*3)+3),]
        barplot(plot_data, main = paste0(name[j+1], " SD : 0.00",sd_val[(i/2)+1]), ylim = range(pretty(c(0,plot_data))), col = c("red","blue","green"), names.arg = c(b_dist[1],b_dist[2]),beside = TRUE)
    }
}

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


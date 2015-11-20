library(ggplot2)
library(gridExtra)

args<-commandArgs(TRUE)

dat <- read.delim(args[1], header=T)
theta <- as.numeric(args[2])
gamma <- as.numeric(args[3])

p_val <- dat$p.val
x <- seq(0, 1, by=0.001)
y <- numeric(length(x))
for (i in 1:length(x)) {y[i] <- length(p_val[p_val<=x[i]])/length(x)}

dat2 <- data.frame(x , y)
qq_plot <- ggplot(dat2, aes(x, y)) + geom_point(colour="red") + theme_bw() + geom_abline(intercept=0, slope=1, alpha=0.5) + theme(legend.position="None")
qq_plot <- qq_plot +  labs(x="Expected", y="Observed", title='QQplot of p-values')

# p_hist <- ggplot(dat, aes(p.value)) + geom_bar(binwidth=0.05, fill='white', colour="Black")
# p_hist <- p_hist + theme_bw() + labs(x="p-value") + ggtitle("Histogram of p-values")

theta_neutral <- ggplot(dat, aes(theta_neutral)) + geom_bar(binwidth=0.00005, fill='white', colour="Black")
theta_neutral <- theta_neutral + theme_bw() + labs(x="theta") + ggtitle("Histogram of theta neutral") + geom_vline(xintercept = theta, colour='Red')

theta_none <- ggplot(dat, aes(theta_none)) + geom_bar(binwidth=0.00005, fill='white', colour="Black")
theta_none <- theta_none + theme_bw() + labs(x="theta") + ggtitle("Histogram of theta none") + geom_vline(xintercept = theta, colour='Red')

gamma_plot <- ggplot(dat, aes(gamma_none)) + geom_bar(binwidth=0.05, fill='white', colour="Black")
gamma_plot <- gamma_plot + theme_bw() + labs(x="gamma") + ggtitle("Histogram of gamma") + geom_vline(xintercept = gamma, colour='Red')

pdf("neutral_test_1.pdf", width=8.2, height=6.5)
grid.arrange(qq_plot, theta_neutral, theta_none, gamma_plot, ncol=2)
dev.off()

sink('summ_stats.txt')
summary(dat[,c(1, 3, 4)])
sink()




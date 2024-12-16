####################################################### R Code for Real Data Application of non-spatial latent model###################################################
# Load necessary packages
library("LaplacesDemon")
library("geoR")
library("mvtnorm")
library("MCMCpack")
library("msm")
# Set seed and load data
set.seed(12345)
data <- read.csv('D:\\Myfiles\\Mythesisresearch\\birthweight.csv')
attach(data)
# Prepare covariates
X_tig <- as.matrix(subset(data[, 2:6], Regions == 1))
X_afa <- as.matrix(subset(data[, 2:6], Regions == 2))
X_amh <- as.matrix(subset(data[, 2:6], Regions == 3))
X_oro <- as.matrix(subset(data[, 2:6], Regions == 4))
X_som <- as.matrix(subset(data[, 2:6], Regions == 5))
X_ben <- as.matrix(subset(data[, 2:6], Regions == 6))
X_snn <- as.matrix(subset(data[, 2:6], Regions == 7))
X_gam <- as.matrix(subset(data[, 2:6], Regions == 8))
X_har <- as.matrix(subset(data[, 2:6], Regions == 9))
X_add <- as.matrix(subset(data[, 2:6], Regions == 10))
X_dir <- as.matrix(subset(data[, 2:6], Regions == 11))
# Rescale covariates
X_tig <- scale(X_tig, center = TRUE, scale = TRUE)
X_afa <- scale(X_afa, center = TRUE, scale = TRUE)
X_amh <- scale(X_amh, center = TRUE, scale = TRUE)
X_oro <- scale(X_oro, center = TRUE, scale = TRUE)
X_som <- scale(X_som, center = TRUE, scale = TRUE)
X_ben <- scale(X_ben, center = TRUE, scale = TRUE)
X_snn <- scale(X_snn, center = TRUE, scale = TRUE)
X_gam <- scale(X_gam, center = TRUE, scale = TRUE)
X_har <- scale(X_har, center = TRUE, scale = TRUE)
X_add <- scale(X_add, center = TRUE, scale = TRUE)
X_dir <- scale(X_dir, center = TRUE, scale = TRUE)
# Set row names
rownames(X_tig) <- 1:nrow(X_tig)
rownames(X_afa) <- 1:nrow(X_afa)
rownames(X_amh) <- 1:nrow(X_amh)
rownames(X_oro) <- 1:nrow(X_oro)
rownames(X_som) <- 1:nrow(X_som)
rownames(X_ben) <- 1:nrow(X_ben)
rownames(X_snn) <- 1:nrow(X_snn)
rownames(X_gam) <- 1:nrow(X_gam)
rownames(X_har) <- 1:nrow(X_har)
rownames(X_add) <- 1:nrow(X_add)
rownames(X_dir) <- 1:nrow(X_dir)
# Prepare X and Y lists
X <- list(X_tig, X_afa, X_amh, X_oro, X_som, X_ben, X_snn, X_gam, X_har, X_add, X_dir)
y_tig <- matrix(subset(data[, 1], Regions == 1), nrow(X_tig), 1)
y_afa <- matrix(subset(data[, 1], Regions == 2), nrow(X_afa), 1)
y_amh <- matrix(subset(data[, 1], Regions == 3), nrow(X_amh), 1)
y_oro <- matrix(subset(data[, 1], Regions == 4), nrow(X_oro), 1)
y_som <- matrix(subset(data[, 1], Regions == 5), nrow(X_som), 1)
y_ben <- matrix(subset(data[, 1], Regions == 6), nrow(X_ben), 1)
y_snn <- matrix(subset(data[, 1], Regions == 7), nrow(X_snn), 1)
y_gam <- matrix(subset(data[, 1], Regions == 8), nrow(X_gam), 1)
y_har <- matrix(subset(data[, 1], Regions == 9), nrow(X_har), 1)
y_add <- matrix(subset(data[, 1], Regions == 10), nrow(X_add), 1)
y_dir <- matrix(subset(data[, 1], Regions == 11), nrow(X_dir), 1)
Y <- list(y_tig, y_afa, y_amh, y_oro, y_som, y_ben, y_snn, y_gam, y_har, y_add, y_dir)
# Initial values for z
z_tig <- matrix(rep(0.3, nrow(X_tig)), nrow(X_tig), 1)
z_afa <- matrix(rep(0.3, nrow(X_afa)), nrow(X_afa), 1)
z_amh <- matrix(rep(0.3, nrow(X_amh)), nrow(X_amh), 1)
z_oro <- matrix(rep(0.3, nrow(X_oro)), nrow(X_oro), 1)
z_som <- matrix(rep(0.3, nrow(X_som)), nrow(X_som), 1)
z_ben <- matrix(rep(0.3, nrow(X_ben)), nrow(X_ben), 1)
z_snn <- matrix(rep(0.3, nrow(X_snn)), nrow(X_snn), 1)
z_gam <- matrix(rep(0.3, nrow(X_gam)), nrow(X_gam), 1)
z_har <- matrix(rep(0.3, nrow(X_har)), nrow(X_har), 1)
z_add <- matrix(rep(0.3, nrow(X_add)), nrow(X_add), 1)
z_dir <- matrix(rep(0.3, nrow(X_dir)), nrow(X_dir), 1)
z <- list(z_tig, z_afa, z_amh, z_oro, z_som, z_ben, z_snn, z_gam, z_har, z_add, z_dir)
# Other variables
m <- 11
n_obs <- matrix(c(length(z_tig), length(z_afa), length(z_amh), length(z_oro), length(z_som), length(z_ben), 
                  length(z_snn), length(z_gam), length(z_har), length(z_add), length(z_dir)), m, 1)
#Hyperparameters
lamda2_prior <- 2
omega_upsilon_prior <- 15
n_iter <- 30000
# Initialize variables
BETA_POST <- upsilon_POST <- sigmasq_upsilon_POST <- NULL
upsilon_posterior <- rep(0.1, m)
sigmasq_upsilon_post <- NULL
z_all <- NULL
# Gibbs Sampler
for (i in 1:n_iter) {
  for (j in 1:m) {
    omega_upsilon_posterior[j] <- (omega_upsilon_prior + 1) # updating degrees of freedom of variance of latent trait
    lamda2_posterior[j] <- (lamda2_prior +
                              ((upsilon_posterior[j])^2) / omega_upsilon_posterior[j]) # update scale parameter of variance of latent trait
    sigmasq_upsilon_post[j] <- rinvchisq(1, omega_upsilon_posterior[j],
                                         lamda2_posterior[j]) # variance of latent trait
  }
  
  # given variance, sample from upsilon (latent traits)
  mean_upsilon_posterior <- list()
  for (t in 1:m) {
    mean_upsilon_posterior[[t]] <- matrix(NA, length(z[[t]]), 1)
    for (t_inner in 1:length(z[[t]])) {
      mean_upsilon_posterior[[t]][t_inner, ] <- (z[[t]][t_inner, ] - X[[t]][t_inner, ] %*% beta)
    }
  }
  
  for (j in 1:m) {
    var_upsilon_posterior[j] <- 1 / (n_obs[j, ] + (1 / sigmasq_upsilon_post[j]))
  }
  # summing over time
  mean_upsilon_posterior_sum <- sapply(mean_upsilon_posterior, function(matrix) sum(matrix))
  
  # updating latent traits
  for (j in 1:m) {
    mean_upsilon_posterior_sum[j] <- mean_upsilon_posterior_sum[j] / (n_obs[j, ] + 1 / sigmasq_upsilon_post[j])
    upsilon_posterior[j] <- rnorm(1, mean_upsilon_posterior_sum[j],
                                  sqrt(var_upsilon_posterior[j]))
  }
  # As I have only one latent trait upsilon for each region
  upsilon_rep <- lapply(upsilon_posterior, function(x) rep(x, n_obs))
  # Convert the list into matrix
  upsilon_rep <- lapply(upsilon_rep, function(x) matrix(x, nrow = length(x), ncol = 1))
  a <- b <- list()
  for (j in 1:m) {
    a[[j]] <- b[[j]] <- matrix(NA, length(z[[j]]), 1)
  }
  # Sampling for each region
  for (j in 1:m) {
    for (t in 1:length(z[[j]])) {
      zt <- z[[j]]
      a[[j]][t, ] <- max(zt[Y[[j]][, 1] < Y[[j]][t, ]])
      b[[j]][t, ] <- min(zt[Y[[j]][t, ] < Y[[j]][, 1]])
      z[[j]][t, ] <- rtnorm(1, X[[j]][t, ] %*% beta + upsilon_rep[[j]][t, ], 1, a[[j]][t, ], b[[j]][t, ])
    }
  }
  muj1 <- array(NA, c(p, p, m))
  muj1_sum <- matrix(NA, p, p)
  muj2 <- array(NA, c(p, 1, m))
  muj2_sum <- matrix(NA, p, 1)
  muj_final <- matrix(NA, p, 1)
  varj_final <- matrix(NA, p, p)
  # Update BETA
  for (j in 1:m) {
    muj1[, , j] <- t(X[[j]]) %*% X[[j]]
    muj1_sum <- apply(muj1, c(1:2), sum)
    muj2[, , j] <- t(X[[j]]) %*% (z[[j]] - upsilon_rep[[j]])
    muj2_sum <- apply(muj2, c(1:2), sum)
    muj_final <- solve(muj1_sum) %*% muj2_sum
    varj_final <- solve(muj1_sum)
  }
  BETA <- mvrnorm(n = 1, mu = muj_final,
                  Sigma = varj_final) # sampling from a multivariate Normal with updated mean and variance 
  BETA_POST <- rbind(BETA_POST, BETA)
  upsilon_POST <- rbind(upsilon_POST, upsilon_posterior)
  sigmasq_upsilon_POST <- rbind(sigmasq_upsilon_POST, sigmasq_upsilon_post)
  print(paste("Iteration", i))
  z_all <- rbind(z_all, z)
}
# Checking the convergence with burn
library(coda)
library(bayesplot)
library(stableGR)
Beta_gibbis_burn <- BETA_POST[seq(1, nrow(BETA_POST[5001:30000, ]), 5), ]
upsilon_gibbis_burn <- upsilon_POST[seq(1, nrow(upsilon_POST[5001:30000, ]), 5), ]
par(mfrow = c(5, 1), mar = c(3, 3, 2, 2))
for (i in 1:5) {
  traceplot(as.mcmc(Beta_gibbis_burn[, i]), col = "blue", main = expression(paste(beta[i])))
}
par(mfrow = c(3, 2), mar = c(7, 3, 3, 1))
for (i in 1:5) {
  densplot(as.mcmc(Beta_gibbis_burn[, i]), col = "blue", type = "h", main = expression(paste(beta[i])))
}
par(mfrow = c(3, 2), mar = c(7, 4, 3, 1))
for (i in 1:5) {
  acf(as.mcmc(Beta_gibbis_burn[, i]), col = "blue", lwd = 2, ylab = "Autocorrelation", ci = FALSE, lag.max = 5000, 
      main = expression(paste(beta[i])))
}
# Summary Statistics
beta_means <- sapply(1:5, function(i) mean(BETA_POST[, i]))
beta_sds <- sapply(1:5, function(i) sd(BETA_POST[, i]))
print("Means:")
print(beta_means)
print("Standard Deviations:")
print(beta_sds)
# HPD Interval
Beta_mcmc_object <- as.mcmc(Beta_gibbis_burn)
Beta_credible_interval <- HPDinterval(Beta_mcmc_object, prob = 0.95)
# Potential Scale Reduction Factor (PSRF) Test
PSRF_beta_test <- stable.GR(Beta_gibbis_burn)
####For upsilon
# Set up the layout for multiple plots
par(mfrow=c(6,1), mar=c(3,3,2,2))
# Plot traceplots for each parameter
for (i in 1:6) {
  traceplot(as.mcmc(upsilon_gibbis_burn[,i]), col="blue", main=expression(paste(upsilon[i])))
}
par(mfrow=c(5,1), mar=c(3,3,2,2))
for (i in 7:11) {
  traceplot(as.mcmc(upsilon_gibbis_burn[,i]), col="blue", main=expression(paste(upsilon[i])))
}
# Plot density plots
par(mfrow=c(3,2), mar=c(7,3,3,1))
for (i in 1:6) {
  densplot(as.mcmc(upsilon_gibbis_burn[,i]), col="blue", type="h", main=expression(paste(upsilon[i])))
}
par(mfrow=c(3,2), mar=c(7,3,3,1))
for (i in 7:11) {
  densplot(as.mcmc(upsilon_gibbis_burn[,i]), col="blue", type="h", main=expression(paste(upsilon[i])))
}
# Plot autocorrelation functions
par(mfrow=c(3,2), mar=c(7,4,3,1))
for (i in 1:6) {
  acf(as.mcmc(upsilon_gibbis_burn[,i]), col="blue", lwd=3, ylab="Autocorrelation", ci=F, lag.max=5000, main=expression(paste(upsilon[i])))
}
par(mfrow=c(3,2), mar=c(7,4,3,1))
for (i in 7:11) {
  acf(as.mcmc(upsilon_gibbis_burn[,i]), col="blue", lwd=3, ylab="Autocorrelation", ci=F, lag.max=5000, main=expression(paste(upsilon[i])))
}
# Calculate posterior means of upsilon
for (i in 1:11) {
  cat("Posterior mean of upsilon[", i, "]: ", mean(upsilon_gibbis_burn[,i]), "\n")
}
# Calculate posterior standard deviations of upsilon
for (i in 1:11) {
  cat("Posterior standard deviation of upsilon[", i, "]: ", sd(upsilon_gibbis_burn[,i]), "\n")
}
# Calculate credible intervals
upsilon_mcmc_object <- as.mcmc(upsilon_gibbis_burn)
upsilon_credible_interval <- HPDinterval(upsilon_mcmc_object, prob = 0.95)
# Perform PSRF test
PSRF_upsilon_test <- stable.GR(upsilon_gibbis_burn)
#####Regional variation 
library(reshape2)
library(viridis)
library(ggplot2)
new_colnames <- c("1","2","3","4","5","6","7","8","9","10","11")
colnames(g) <- new_colnames
g=data.frame(sigma_gibbis_burn)
data <- var(g[sapply(g,is.numeric)])
data1 <- melt(data)
data2=read.csv(('D:\\Myfiles\\Mythesisresearch\\Variation.csv'))
data2=data.frame(data2)
colnames(data2)=c("Var1","Var2", "Variations")
x_labels <- c("Tigray", "Afar", "Amhara", "Oromia", "Somale", "Benishangul", "SNN", "Gambela", "Harari", "Addis Abeba", "Dire Dawa")
y_labels <- c("Tigray", "Afar", "Amhara", "Oromia", "Somale", "Benishangul", "SNN", "Gambela", "Harari", "Addis Abeba", "Dire Dawa")
ggplot(data2, aes(x = Var1, y = Var2, fill = Variations)) +
  geom_tile() +scale_fill_viridis(discrete = FALSE)+geom_tile()+
  labs(title = "",
       x = "Regions",
       y = "Regions")+scale_x_discrete(labels = x_labels)+
  scale_y_discrete(labels = y_labels)
#########Scaterplot
library(foreign)
library(ggplot2)
library(VGAM)
library(extraDistr)
library(ggpubr)
mydata=read.spss('D:\\Myfiles\\Mythesisresearch\\LBWregion.sav')
attach(mydata)
Regiondata<-data.frame(childweight,region,Percent)
r=ggplot(data=Regiondata,mapping=aes(x=region,Number,y=Percent,col=childweight))+geom_point(size=3)+ylim(0,95)+labs(x="Regions",y="Percentage of children",colour="Birth weight:")
agedata=read.spss('D:\\Myfiles\\Mythesisresearch\\LBWage.sav')
attach(agedata)
agedata<-data.frame(childweight,motherage,Percent)
a=ggplot(data=agedata,mapping=aes(x=motherage,Number,y=Percent,col=childweight))+geom_point(size=3)+labs(x="Mother Age",y="Percentage of children ",colour="Birth weight:")
ggarrange(r, a, ncol=1, nrow=2, common.legend = TRUE,legend="bottom")
ancdata=read.spss('D:\\Myfiles\\Mythesisresearch\\LBWanc.sav')
attach(ancdata)
ancdata<-data.frame(childweight,ancvisit,Percent)
an=ggplot(data=ancdata,mapping=aes(x=ancvisit,Number,y=Percent,col=childweight))+geom_point(size=3)+ylim(0,60)+labs(x="Number of ANC visits",y="Percentage of children",colour="Birth weight:")
birthorderdata=read.spss('D:\\Myfiles\\Mythesisresearch\\LBWbirthorder.sav')
attach(birthorderdata)
birthorderdata<-data.frame(childweight,birthorder,Percent)
br=ggplot(data=birthorderdata,mapping=aes(x=birthorder,Number,y=Percent,col=childweight))+geom_point(size=3)+ylim(0,50)+labs(x="Birth order",y="Percentage of children",colour="Birth weight:")
ggarrange(an, br, ncol=1, nrow=2, common.legend = TRUE,legend="bottom")
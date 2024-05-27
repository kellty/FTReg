source('func.R')
load("ADHD_prep.RData")
pT <- Xtsr_ADHD@modes[2]
rho_list_ADHD <- 10^seq(0.5,3.5,0.2)
GCV_ADHD <- AIC_ADHD <- BIC_ADHD <- rho_ADHD <- NULL
for (rT in 2:6) {
r_ADHD <- c(rT,2,2,2)
cv_ADHD <- ftreg.gcv(y_ADHD, Xtsr_ADHD, breaks_ADHD, r_ADHD, M_ADHD, rho_list_ADHD, halforder=2, iter_max=160)
idx <- which.min(cv_ADHD$GCV)
GCV_ADHD <- c(GCV_ADHD, cv_ADHD$GCV[idx])
AIC_ADHD <- c(AIC_ADHD, cv_ADHD$AIC[idx])
BIC_ADHD <- c(BIC_ADHD, cv_ADHD$BIC[idx])
rho_ADHD <- c(rho_ADHD, rho_list_ADHD[idx])
}
r_ADHD <- c((2:6)[which.min(GCV_ADHD)],2,2,2)
rho_cv_ADHD <- rho_ADHD[which.min(GCV_ADHD)]
cv10_ADHD <- ftreg.cv(y_ADHD, Xtsr_ADHD, breaks_ADHD, r_ADHD, M_ADHD, c(rho_cv_ADHD,1e-2), halforder=2, cv_fold=10, iter_max=160)
ftreg.results_ADHD <- ftreg(y_ADHD, Xtsr_ADHD, breaks_ADHD, r_ADHD, M_ADHD, rho_cv_ADHD, halforder=2, iter_max=160)
R2 <- 1 - sum(ftreg.results_ADHD$res^2) / sum((y_ADHD-mean(y_ADHD))^2)
Theta <- ftreg.results_ADHD$nsbasis(breaks_ADHD) %*% ftreg.results_ADHD$regcoef_mat
regcoef_tsr_sample <- k_fold(Theta, 1, c(pT,8,8,4))
rep_mat <- function(ntot, nsub){
  H <- matrix(0, nrow=ntot, ncol=ntot/nsub)
  for (j in 1:ncol(H)) {
    H[((j-1)*nsub+1):(j*nsub),j] <- 1/nsub
  }
  return(H)
}
rep_z <- cbind(rbind(rep_mat(24,8),matrix(0,9,3)), c(rep(0,24),rep(1/9,9))) # rep_mat(33,4)
rep_mat_list <- list(rep_mat(64,8),rep_mat(64,8),rep_z)
regcoef_tsr_sample <- ttl(regcoef_tsr_sample, rep_mat_list, 2:4)
regcoef_tsr_sample <- k_fold(k_unfold(regcoef_tsr_sample, 1), 4, c(64,64,33,pT))
library(RNifti)
nim <- readNifti("ADHD/Peking_1038415_1/rest_1/NIfTI/rest.nii.gz")[,,,101:(100+pT)]
library(oro.nifti)
x <- 37;  y <- 45;  z <- 21;  idx <- 3*(1:16)
tmp_x <- abs(regcoef_tsr_sample@data[x,,,])
tmp_x <- ifelse(tmp_x>quantile(tmp_x,0.8), regcoef_tsr_sample@data[x,,,], NA)
tmp_x[1,1,] <- ifelse(is.na(tmp_x[1,1,]), 0, tmp_x[1,1,])
tmp_y <- abs(regcoef_tsr_sample@data[,y,,])
tmp_y <- ifelse(tmp_y>quantile(tmp_y,0.8), regcoef_tsr_sample@data[,y,,], NA)
tmp_y[1,1,] <- ifelse(is.na(tmp_y[1,1,]), 0, tmp_y[1,1,])
tmp_z <- abs(regcoef_tsr_sample@data[,,z,])
tmp_z <- ifelse(tmp_z>quantile(tmp_z,0.8), regcoef_tsr_sample@data[,,z,], NA)
tmp_z[1,1,] <- ifelse(is.na(tmp_z[1,1,]), 0, tmp_z[1,1,])
setEPS()
postscript("ADHD_slice_x.eps")
overlay(nim[x,,,idx], tmp_x[,,idx])
dev.off()
setEPS()
postscript("ADHD_slice_y.eps")
overlay(nim[,y,,idx], tmp_y[,,idx])
dev.off()
setEPS()
postscript("ADHD_slice_z.eps")
overlay(nim[,,z,idx], tmp_z[,,idx])
dev.off()
# for (x in seq(5,64,8)) {
#   png(file=paste0("ADHD_slice_x",as.character(x),".png"))
#   tmp <- abs(regcoef_tsr_sample@data[x,,,])
#   tmp <- ifelse(tmp>quantile(tmp,0.8), regcoef_tsr_sample@data[x,,,], NA)
#   tmp[1,1,] <- ifelse(is.na(tmp[1,1,]), 0, tmp[1,1,])
#   overlay(nim[x,,,], tmp)
#   dev.off()
# }
# for (y in seq(5,64,8)) {
#   png(file=paste0("ADHD_slice_y",as.character(y),".png"))
#   tmp <- abs(regcoef_tsr_sample@data[,y,,])
#   tmp <- ifelse(tmp>quantile(tmp,0.8), regcoef_tsr_sample@data[,y,,], NA)
#   tmp[1,1,] <- ifelse(is.na(tmp[1,1,]), 0, tmp[1,1,])
#   overlay(nim[,y,,], tmp)
#   dev.off()
# }
# for (z in seq(5,33,8)) {
#   png(file=paste0("ADHD_slice_z",as.character(z),".png"))
#   tmp <- abs(regcoef_tsr_sample@data[,,z,])
#   tmp <- ifelse(tmp>quantile(tmp,0.8), regcoef_tsr_sample@data[,,z,], NA)
#   tmp[1,1,] <- ifelse(is.na(tmp[1,1,]), 0, tmp[1,1,])
#   overlay(nim[,,z,], tmp)
#   dev.off()
# }

setEPS()
postscript("ADHD_t.eps")
# png("ADHD_t.png")
x <- 37;  y <- 36;  z <- 21
coef_ADHD <- solve(t(ftreg.results_ADHD$nsbasis(breaks_ADHD))%*%ftreg.results_ADHD$nsbasis(breaks_ADHD),
                   t(ftreg.results_ADHD$nsbasis(breaks_ADHD))%*%regcoef_tsr_sample@data[x,y,z,])
plot((100:900)/1000, ftreg.results_ADHD$nsbasis((100:900)/1000)%*%coef_ADHD, xlab='Time', ylab="Effect on ADHD", type='l', lwd=2, col='red')
x <- 37;  y <- 36;  z <- 14
coef_ADHD <- solve(t(ftreg.results_ADHD$nsbasis(breaks_ADHD))%*%ftreg.results_ADHD$nsbasis(breaks_ADHD),
                   t(ftreg.results_ADHD$nsbasis(breaks_ADHD))%*%regcoef_tsr_sample@data[x,y,z,])
lines((100:900)/1000, ftreg.results_ADHD$nsbasis((100:900)/1000)%*%coef_ADHD, lwd=2, lty=2, col='orange')
x <- 37;  y <- 50;  z <- 14
coef_ADHD <- solve(t(ftreg.results_ADHD$nsbasis(breaks_ADHD))%*%ftreg.results_ADHD$nsbasis(breaks_ADHD),
                   t(ftreg.results_ADHD$nsbasis(breaks_ADHD))%*%regcoef_tsr_sample@data[x,y,z,])
lines((100:900)/1000, ftreg.results_ADHD$nsbasis((100:900)/1000)%*%coef_ADHD, lwd=2, lty=4, col='black')
legend('topright', c("cortical surfaces", "white matter", "cerebellum"), lwd=2, lty=c(1,2,4), col=c('red','orange','black'))
dev.off()

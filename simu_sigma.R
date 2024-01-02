set.seed(4321)
result_sigma <- NULL
n <- 500
r_fit <- r
sigma_list <- c(0.02,0.05,0.08,0.1)
rho_list <- c(0,10^seq(-8.6,-4,1)) # tuning parameter
for (sig in sigma_list) {
  RISE_s <- RISE_0 <- yabs <- NULL
for (montecarlo in 1:50) {
Xftsr <- list(basis=Xbasis, score=as.tensor(
              array(rnorm(n*K*prod(p[-1])), dim=c(n,K,p[-1]))) ) # Covariate
Xtsr <- ttm(Xftsr$score, Xftsr$basis(breaks), 2) + as.tensor(
        array(rnorm(n*prod(p), sd=0.05), dim=c(n,p)) )
M <- matrix(0, nrow=n, ncol=1)
gam <- rnorm(ncol(M), 0.666)
eps <- rnorm(n, sd=sig) # Error
y <- M %*% gam + eps # Response
for (idx in 1:ngrid) {
  Xtsr_at_grid <- k_unfold(ttm(Xftsr$score, Xftsr$basis(grid_point[idx]), 2), 1)@data
  regcoef_at_grid <- c(regcoef_vec$nsbasis(grid_point[idx]) %*% regcoef_vec$regcoef_mat)
  y <- Xtsr_at_grid %*% regcoef_at_grid /ngrid + y # Response
}
ISE <- NULL
for (rho in rho_list) {
  print("  Set tuning parameter ...")
  ftreg_tmp <- ftreg(y, Xtsr, breaks, r_fit, M, rho, halforder=2, regcoef_mat)
  err <- 0
  if (sum(regcoef_mat^2)) for (idx in 1:ngrid) {
    regcoef_at_grid <- c(regcoef_vec$nsbasis(grid_point[idx]) %*% regcoef_mat)
    regcoef_est_at_grid <- c(ftreg_tmp$nsbasis(grid_point[idx]) %*% ftreg_tmp$regcoef_mat)
    err <- err + sum((regcoef_at_grid - regcoef_est_at_grid)^2)/ngrid
  }
  ISE <- c(ISE, err)
}
yabs <- c(yabs, mean(abs(y)))
print(paste("mean of |y| =", mean(abs(y)), ",  σ =", sig))
print(paste0("lgρ=", log10(rho_list), ", RISE=", round(ISE/ground, 7)))
RISE_s <- c(RISE_s, min(ISE)/ground)
RISE_0 <- c(RISE_0, ISE[1]/ground)
}
result_sigma <- rbind(result_sigma, c(mean(RISE_s),sqrt(var(RISE_s)/50),mean(RISE_0),sqrt(var(RISE_0)/50),mean(yabs),sqrt(var(yabs)/50)))
}
SNR <- result_sigma[,5]/sigma_list
print(paste("σ =", sigma_list, ",  SNR =", SNR))
row.names(result_sigma) <- sigma_list
# save(result_sigma, file="simu_sigma.RData")
setEPS()
postscript("simu_RISE_sigma.eps")
# png("simu_RISE_sigma.png")
plot(sigma_list, result_sigma[,1], ylim=range(result_sigma[,c(1,3)]), xlab=expression(sigma), ylab='RISE', type='b', lwd=2)
lines(sigma_list, result_sigma[,3], type='b', lty=2, lwd=2)
# segments(sigma_list, result_sigma[,3]-result_sigma[,4], sigma_list, result_sigma[,3]+result_sigma[,4], lty=1)
# arrows(sigma_list, result_sigma[,3]-result_sigma[,4], sigma_list, result_sigma[,3]+result_sigma[,4], code=3, angle=90, length=0.1)
legend('topleft', c('functional','tabular'), lty=1:2, lwd=2)
dev.off()

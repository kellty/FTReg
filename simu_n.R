set.seed(4321)
result_n <- NULL
r_fit <- r
n_list <- c(150,200,300,500,800,1200,1700)
rho_list <- c(0,10^seq(-8.6,-4,1)) # tuning parameter
for (n in n_list) {
  RISE <- RISE_s <- RISE_0 <- NULL
for (montecarlo in 1:50) {
Xftsr <- list(basis=Xbasis, score=as.tensor(
              array(rnorm(n*K*prod(p[-1])), dim=c(n,K,p[-1]))) ) # Covariate
Xtsr <- ttm(Xftsr$score, Xftsr$basis(breaks), 2) + as.tensor(
        array(rnorm(n*prod(p), sd=0.05), dim=c(n,p)) )
M <- matrix(0, nrow=n, ncol=1)
gam <- rnorm(ncol(M), 0.666)
eps <- rnorm(n, sd=0.1) # Error
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
print(paste0("lgρ=", log10(rho_list), ", RISE=", round(ISE/ground, 7)))
RISE <- rbind(RISE, ISE/ground)
RISE_s <- c(RISE_s, min(ISE)/ground)
RISE_0 <- c(RISE_0, ISE[1]/ground)
}
result_n <- rbind(result_n, c(mean(RISE_s),sqrt(var(RISE_s)/50),mean(RISE_0),sqrt(var(RISE_0)/50)))
}
row.names(result_n) <- n_list
# save(result_n, file="simu_n.RData")
setEPS()
postscript("simu_RISE_n.eps")
# png("simu_RISE_n.png")
plot(n_list[-(1:2)], result_n[-(1:2),1], log='x', xlab=expression(n), ylab='RISE', type='b', lwd=2)
# lines(n_list, result_n[,3], type='b', lty=2, lwd=2)
segments(n_list[-(1:2)], result_n[-(1:2),1]-result_n[-(1:2),2], n_list[-(1:2)], result_n[-(1:2),1]+result_n[-(1:2),2], lty=1)
arrows(n_list[-(1:2)], result_n[-(1:2),1]-result_n[-(1:2),2], n_list[-(1:2)], result_n[-(1:2),1]+result_n[-(1:2),2], code=3, angle=90, length=0.1)
# fit_n <- lm(log(result_n[,1])~log(n_list))
lines(1:2000, 0.96*result_n[length(n_list),1]*(max(n_list)/1:2000)^0.5, lty=3, lwd=2)
legend('topright', expression(RISE %prop% n^{-1/2}), lty=3, lwd=1)
dev.off()

set.seed(4321)
result <- NULL
r_fit <- r
n_list <- c(300,500,800,1200,1700)
for (n in n_list) {
  GCV_s <- RISE_s <- GCV <- RISE <- NULL
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
cv <- ftreg.gcv(y, Xtsr, breaks, r_fit, M, rho_list, halforder=2, regcoef_mat, ns.results$nsbasis)
print(paste0("lgρ=", log10(rho_list), ", GCV=", round(cv$GCV, 7), ", RISE=", round(cv$ISE/ground, 7)))
GCV <- rbind(GCV, cv$GCV);  RISE <- rbind(RISE, cv$ISE/ground)
GCV_s <- c(GCV_s, min(cv$GCV));  RISE_s <- c(RISE_s, cv$ISE[which.min(cv$GCV)]/ground)
# rho_cv <- rho_list[which.min(cv$GCV)]
# ftreg.results <- ftreg(y, Xtsr, breaks, r_fit, M, rho_cv, halforder=2, regcoef_mat)
# rho_best <- rho_list[which.min(cv$ISE)]
# ftreg.results_best <- ftreg(y, Xtsr, breaks, r_fit, M, rho_best, halforder=2, regcoef_mat)
}
result <- rbind(result, c(mean(GCV_s),sqrt(var(GCV_s)/50),mean(RISE_s),sqrt(var(RISE_s)/50)))
}
# save(result, file="simu_n.RData")
setEPS()
postscript("simu_RISE_n.eps")
plot(n_list, result[,3], log='xy', xlab=expression(n), ylab='RISE', type='b', lwd=3)
segments(n_list, result[,3]-result[,4], n_list, result[,3]+result[,4], lty=1)
arrows(n_list, result[,3]-result[,4], n_list, result[,3]+result[,4], code=3, angle=90, length=0.1)
fit_n <- lm(log(result[,3])~log(n_list))
lines(n_list, exp(fit_n$fitted.values), lty=2)
legend('topright', paste('slope =',round(fit_n$coefficients[2],2)), lty=2)
dev.off()

set.seed(4321)
result <- result_IC <- NULL
n <- 500
for (r_fit in list(c(2,3,3),c(2,2,2),c(2,3,2),c(2,4,2),c(2,4,3),c(3,2,2),c(3,3,2),c(3,3,3),c(3,4,2))) {
  GCV_s <- RISE_s <- GCV <- RISE <- NULL
  AIC_s <- BIC_s <- AIC <- BIC <- RISE_sA <- RISE_sB <- NULL
for (montecarlo in 1:100) {
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
AIC <- rbind(AIC, cv$AIC);  BIC <- rbind(BIC, cv$BIC)
GCV_s <- c(GCV_s, min(cv$GCV)); RISE_s <- c(RISE_s, cv$ISE[which.min(cv$GCV)]/ground)
AIC_s <- c(AIC_s, min(cv$AIC)); RISE_sA <- c(RISE_sA, cv$ISE[which.min(cv$AIC)]/ground)
BIC_s <- c(BIC_s, min(cv$BIC)); RISE_sB <- c(RISE_sB, cv$ISE[which.min(cv$BIC)]/ground)
# rho_cv <- rho_list[which.min(cv$GCV)]
# ftreg.results <- ftreg(y, Xtsr, breaks, r_fit, M, rho_cv, halforder=2, regcoef_mat)
# rho_best <- rho_list[which.min(cv$ISE)]
# ftreg.results_best <- ftreg(y, Xtsr, breaks, r_fit, M, rho_best, halforder=2, regcoef_mat)
}
result <- rbind(result, c(mean(GCV_s),sqrt(var(GCV_s)/100),mean(RISE_s),sqrt(var(RISE_s)/100)))
result_IC <- rbind(result_IC, c(mean(AIC_s),sqrt(var(AIC_s)/100),mean(RISE_sA),sqrt(var(RISE_sA)/100),
                                mean(BIC_s),sqrt(var(BIC_s)/100),mean(RISE_sB),sqrt(var(RISE_sB)/100)))
}
save(result, file="simu_r.RData")

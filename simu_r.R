set.seed(4321)
result_r <- NULL
n <- 500
r_list <- list(c(2,3,3),c(2,2,2),c(2,3,2),c(2,4,2),c(2,4,3),c(3,2,2),c(3,3,2),c(3,3,3),c(3,4,2),c(3,4,3),c(3,4,4),c(4,4,3),c(4,4,4),c(5,5,5))
for (r_fit in r_list) {
  RISE_s <- GCV_s <- AIC_s <- BIC_s <- NULL
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
idx <- which.min(cv$GCV)
RISE_s <- c(RISE_s, cv$ISE[idx]/ground)
GCV_s <- c(GCV_s, cv$GCV[idx])
AIC_s <- c(AIC_s, cv$AIC[idx])
BIC_s <- c(BIC_s, cv$BIC[idx])
}
result_r <- rbind(result_r, c(mean(RISE_s),sqrt(var(RISE_s)/100),mean(GCV_s),sqrt(var(GCV_s)/100),
                              mean(AIC_s),sqrt(var(AIC_s)/100),mean(BIC_s),sqrt(var(BIC_s)/100)))
}
row.names(result_r) <- as.character(r_list)
# save(result_r, file="simu_r.RData")

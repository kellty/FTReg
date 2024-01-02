set.seed(1234)
result_p <- NULL
n <- 500
r_fit <- r
rho_list <- c(0,10^seq(-10.6,-3,0.5)) # tuning parameter
Stsr <- as.tensor(array(rnorm(prod(r)), dim=r)) # Core tensor
p_list <- c(3,6,9,13,17)
for (pT in p_list) {
p <- c(pT,8,8) # c(pT,p1,...,pD) with D>1
breaks <- (0.5:(p[1]-0.5))/p[1]
ns.results <- natspline(breaks)
Umat_list <- list()
for (d in 1:Dplus1) {
  Umat_tmp <- qr.Q(qr(matrix(rnorm(p[d]*r[d]), ncol=r[d])))
  Umat_list <- append(Umat_list, list(Umat_tmp))
}
THETAtsr <- ttl(Stsr, Umat_list, 1:Dplus1)
regcoef_mat <- solve(t(ns.results$ns_at_knots)%*%ns.results$ns_at_knots,
                     t(ns.results$ns_at_knots)%*%k_unfold(THETAtsr, 1)@data)
regcoef_vec <- list(nsbasis=ns.results$nsbasis, regcoef_mat=regcoef_mat)
ground <- 0
for (idx in 1:ngrid) {
  regcoef_at_grid <- c(regcoef_vec$nsbasis(grid_point[idx]) %*% regcoef_vec$regcoef_mat)
  ground <- ground + sum(regcoef_at_grid^2)/ngrid
}
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
result_p <- rbind(result_p, c(mean(RISE_s),sqrt(var(RISE_s)/50),mean(RISE_0),sqrt(var(RISE_0)/50)))
}
row.names(result_p) <- p_list
# save(result_p, file="simu_p.RData")
setEPS()
postscript("simu_RISE_p.eps")
# png("simu_RISE_p.png")
plot(p_list, result_p[,1], xlab=expression(p[0]), ylab='RISE', type='b', lwd=2)
# lines(p_list, result_p[,3], type='b', lty=3, lwd=2, pch=2)
segments(p_list, result_p[,1]-result_p[,2], p_list, result_p[,1]+result_p[,2], lty=1)
arrows(p_list, result_p[,1]-result_p[,2], p_list, result_p[,1]+result_p[,2], code=3, angle=90, length=0.1)
# legend('top', c('functional','tabular'), lty=c(1,3), lwd=2, pch=1:2)
dev.off()

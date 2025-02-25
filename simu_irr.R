source('func.R')
set.seed(1234)
K <- 30 # num of basis func
Xbasis <- function(u){
  basis_mat <- matrix(NA, nrow=length(u), ncol=K)
  for (k in 1:K) {
    basis_mat[,k] <- sin(k*pi*u)/k
  }
  return(basis_mat)
}
p <- c(12,8,8) # c(pT,p1,...,pD) with D>1
# r <- ceiling(p/10); r[r==1] <- 2 # prod(r)>=max(r)^2
r <- c(2,3,3) # true Tucker rank
rho_list <- 10^seq(-8.8,-4,1) # tuning parameter
Dplus1 <- length(r)
breaks <- c((1/3)*(0.5:(5-0.5))/5, (1/3)+(1/3)*(0.5:(2-0.5))/2, (2/3)+(1/3)*(0.5:(5-0.5))/5)
ns.results <- natspline(breaks)
Stsr <- as.tensor(array(rnorm(prod(r)), dim=r)) # Core tensor
Umat_list <- list()
for (d in 1:Dplus1) {
  Umat_tmp <- qr.Q(qr(matrix(rnorm(p[d]*r[d]), ncol=r[d])))
  Umat_list <- append(Umat_list, list(Umat_tmp))
}
THETAtsr <- ttl(Stsr, Umat_list, 1:Dplus1)
regcoef_mat <- solve(t(ns.results$ns_at_knots)%*%ns.results$ns_at_knots,
                     t(ns.results$ns_at_knots)%*%k_unfold(THETAtsr, 1)@data)
regcoef_vec <- list(nsbasis=ns.results$nsbasis, regcoef_mat=regcoef_mat)
ngrid <- 1e3
grid_point <- (0.5:(ngrid-0.5))/ngrid
ground <- 0
for (idx in 1:ngrid) {
  regcoef_at_grid <- c(regcoef_vec$nsbasis(grid_point[idx]) %*% regcoef_vec$regcoef_mat)
  ground <- ground + sum(regcoef_at_grid^2)/ngrid
}

set.seed(4321)
n <- 500
r_fit <- r
GCV <- RISE <- err <- NULL
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
cv <- ftreg.gcv(y, Xtsr, breaks, r_fit, M, c(0,rho_list), halforder=2, regcoef_mat, ns.results$nsbasis)
print(paste0("lgρ=", log10(c(0,rho_list)), ", GCV=", round(cv$GCV, 7), ", RISE=", round(cv$ISE/ground, 7)))
GCV <- rbind(GCV, cv$GCV);  RISE <- rbind(RISE, cv$ISE/ground)
rho_cv <- c(0,rho_list)[which.min(cv$GCV)]
ftreg.results <- ftreg(y, Xtsr, breaks, r_fit, M, rho_cv, halforder=2, regcoef_mat)
# rho_best <- rho_list[which.min(cv$ISE)]
# ftreg.results_best <- ftreg(y, Xtsr, breaks, r_fit, M, rho_best, halforder=2, regcoef_mat)
err[[montecarlo]] <- ftreg.results$error_matrix[,4]
}
pdf("simu_irr_GCV_rho.pdf", width=4, height=5)
plot(rho_list, colMeans(GCV[,-1]), log='x', xlab=expression(rho), ylab='GCV', type='b', lwd=3)
abline(h=mean(GCV[,1]), lty=2, lwd=2)
legend('topleft', expression(rho==0), lty=2, lwd=2)
# for (i in 1:nrow(GCV)) {
#   lines(rho_list, GCV[i,], type='l', lty=2)
# }
dev.off()
pdf("simu_irr_RISE_rho.pdf", width=4, height=5)
plot(rho_list, colMeans(RISE[,-1]), log='x', xlab=expression(rho), ylab='RISE', type='b', lwd=3)
abline(h=mean(RISE[,1]), lty=2, lwd=2)
legend('topleft', expression(rho==0), lty=2, lwd=2)
# for (i in 1:nrow(RISE)) {
#   lines(rho_list, RISE[i,], type='l', lty=2, col=8)
# }
dev.off()
pdf("simu_irr_conv.pdf", width=4, height=5)
plot(1:10, err[[1]][1:10], log='y', type='l', lwd=2,
     xlab='iteration k', ylab=expression(abs(abs(Theta^{k}-Theta))[F]/abs(abs(Theta))[F]))
for (montecarlo in 2:5) {
  lines(1:10, err[[montecarlo]][1:10], lwd=2, col=montecarlo)
}
dev.off()

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
rho_list <- 10^seq(-8.6,-4,1) # tuning parameter
Dplus1 <- length(r)
breaks <- (0.5:(p[1]-0.5))/p[1]
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

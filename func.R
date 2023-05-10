library(rTensor)

tpartunfold <- function(X, idx_preserved){# unfold the modes for X
  # Example: X: 5*3*4*10 tensor; tpartunfold(X, c(1,4)) yields 5*10*12 tensor
  p <- X@modes
  idx_unfolded <- setdiff(1:length(p), idx_preserved)
  idx_hat <- c(idx_preserved, idx_unfolded)
  p_unfolded <- prod(p[idx_unfolded])
  Y <- array(aperm(X@data, idx_hat), dim=c(p[idx_preserved], p_unfolded))
  return(as.tensor(Y))
}

tsrreg_RGN <- function(yc, Ztsr, r, rhoAmat=0, THETAtsr=as.tensor(0),
                       THETAtsr_ini=as.tensor(0), iter_max=60, tol=1e-8){
# yc: centered response n-vector;
# Ztsr: covariate n*pT*p1*...*pD tensor;
# r: Tucker rank c(rT,r1,...,rD)
# rhoAmat: roughness penalty pT*pT matrix;
# THETAtsr: true parameter pT*p1*...*pD tensor;
# THETAtsr_ini: initialized parameter pT*p1*...*pD tensor;
# iter_max: the maximum number of iterations;
# tol: tolerence for successful recovery /stoping criteria;
# Output estimated parameter, residual vector, running time and estimation error
  Dplus1 <- length(r) # Dplus1>2
  p <- Ztsr@modes # c(n,pT,p1,...,pD)
  n <- p[1]
  p <- p[-1] # c(pT,p1,...,pD)
  p.prod <- prod(p[-1])
  r.prod <- prod(r[-1])
  r.prod.T <- r.prod*r[1]
  pminusr <- p - r
  r.pminusr <- r * pminusr
  LstarAL <- diag(r.prod.T+sum(r.pminusr))
  idx_core <- 1:r.prod.T
  if (any(pminusr==0) || r.prod.T < max(r)^2) {
    print("Error: non-conformable rank")
    return(NA)
  }
  THETAfnorm <- fnorm(THETAtsr)
  if (THETAfnorm) {
    eps <- yc - k_unfold(Ztsr, 1)@data %*% c(THETAtsr@data)
    print(paste("infeasible R^2 =", 1-sum(eps^2)/sum(yc^2)))
  }
  if (fnorm(THETAtsr_ini)==0) {
    THETAtsr_ini <- as.tensor(array( t(t(yc)%*%k_unfold(Ztsr,1)@data), dim=p))
    # THETAtsr_ini <- as.tensor(array(rnorm(p[1]*p.prod), dim=p))
  }
  hosvd_result <- hosvd(THETAtsr_ini, r)
  THETAnew <- hosvd_result$est  # initialization
  Stsr <- hosvd_result$Z        # core tensor
  Umat_list <- hosvd_result$U   # loading matrices
  errmat <- NULL
  t0 <- proc.time()
  for (iter in 1:iter_max){
    Umatperp_list <- Umattran_list <- Vmat_list <- Wmat_list <- hatDlist <- list()
    for (d in 1:Dplus1) {
      Umatperp_tmp <- as.matrix(qr.Q(qr(Umat_list[[d]]),complete=TRUE)[,(r[d]+1):p[d]])
      Umatperp_list <- append(Umatperp_list, list(Umatperp_tmp))
      Umattran_list <- append(Umattran_list, list(t(Umat_list[[d]])))
      Vmat_tmp <- as.matrix(qr.Q(qr(t(k_unfold(Stsr, d)@data))))
      Vmat_list <- append(Vmat_list, list(Vmat_tmp))
      Wmat_tmp <- kronecker_list(Umat_list[setdiff(Dplus1:1, d)]) %*% Vmat_tmp
      Wmat_list <- append(Wmat_list, list( Wmat_tmp ))
    }
    ZL <- matrix(ttl(Ztsr, Umattran_list, 2:(Dplus1+1))@data, nrow=n)
    for (d in 1:Dplus1) {
      # ZL.d <- ttm(tpartunfold(ttm(Ztsr, t(Umatperp_list[[d]]), d+1),
      #                         c(1,d+1)), t(Wmat_list[[d]]), 3)
      tmp <- ttl(Ztsr, Umattran_list[-d], setdiff(2:(Dplus1+1), d+1))
      tmp <- ttm(tmp, t(Umatperp_list[[d]]), d+1)
      ZL.d <- ttm(tpartunfold(tmp, c(1,d+1)), t(Vmat_list[[d]]), 3)
      ZL <- cbind(ZL, k_unfold(ZL.d, 1)@data)
    }
    UTtAUT <- t(Umat_list[[1]]) %*% rhoAmat %*% Umat_list[[1]]
    LstarAL[idx_core,idx_core] <- kronecker(diag(r.prod), UTtAUT)
    idx_tmp <- r.prod.T + r.pminusr[1]
    idx <- (r.prod.T+1):idx_tmp
    tmp <- kronecker(Vmat_list[[1]], t(Umat_list[[1]])%*%rhoAmat%*%Umatperp_list[[1]])
    LstarAL[idx_core,idx] <- tmp
    LstarAL[idx,idx_core] <- t(tmp)
    tmp <- t(Umatperp_list[[1]]) %*% rhoAmat %*% Umatperp_list[[1]]
    LstarAL[idx,idx] <- kronecker(diag(r[1]), tmp)
    for (d in 2:Dplus1) {
      idx <- (idx_tmp+1):(idx_tmp+r.pminusr[d])
      tmp <- t(Vmat_list[[d]])%*%kronecker(diag(r.prod/r[d]), UTtAUT)%*%Vmat_list[[d]]
      LstarAL[idx,idx] <- kronecker(tmp, diag(pminusr[d]))
      idx_tmp <- idx_tmp + r.pminusr[d]
    }
    # penalized least squares estimates for the coefficients
    estCoef <- solve(t(ZL)%*%ZL+n*LstarAL, t(ZL)%*%yc)
    hatC <- as.tensor(array(estCoef[idx_core], r))
    idx_tmp <- r.prod.T
    for (d in 1:Dplus1) {
      idx <- (idx_tmp+1):(idx_tmp+r.pminusr[d])
      hatDd <- matrix(estCoef[idx], ncol=r[d])
      hatDlist <- append(hatDlist, list( hatDd ))
      idx_tmp <- idx_tmp + r.pminusr[d]
    }
    THETAold <- THETAnew
    THETAnew <- ttl(hatC, Umat_list, 1:Dplus1)
    for (d in 1:Dplus1) {
      tmp <- Umatperp_list[[d]] %*% hatDlist[[d]] %*% t(Wmat_list[[d]])
      THETAnew <- THETAnew + k_fold(tmp, d, modes=p)
      # m <- r; m[d] <- p[d]
      # tmp <- Umatperp_list[[d]] %*% hatDlist[[d]] %*% t(Vmat_list[[d]])
      # tmp <- k_fold(tmp, d, modes=m)
      # THETAnew <- THETAnew + ttl(tmp, Umat_list[-d], setdiff(1:Dplus1, d))
    }
    hosvd_result <- hosvd(THETAnew, r)
    Stsr <- hosvd_result$Z
    Umat_list <- hosvd_result$U
    THETAnew <- hosvd_result$est
    res <- yc - k_unfold(Ztsr, 1)@data %*% c(THETAnew@data)
    Rsquared <- 1 - sum(res^2)/sum(yc^2)
    THETA_rela_err <- fnorm(THETAnew-THETAold)/fnorm(THETAold)
    err_report <- paste0(':  rela_err_updating=', THETA_rela_err,';  R^2=', Rsquared)
    if (THETAfnorm) {
      THETA_rela_err_true <- fnorm(THETAnew-THETAtsr)/THETAfnorm
      err_report <- paste0(err_report,';  rela_err_true=',THETA_rela_err_true)
    }
    print(paste0('iteration ',iter,'/',iter_max,err_report))
    errmat <- rbind(errmat, c( (proc.time()-t0)[3], THETA_rela_err,
                               Rsquared, ifelse(THETAfnorm,THETA_rela_err_true,NA) ))
    if (THETA_rela_err < tol) {break}
  }
  hatMat <- ZL%*%solve(t(ZL)%*%ZL+n*LstarAL, t(ZL))
  colnames(errmat) <- c('running_time', 'relative_error_updating',
                        'R^2', 'relative_error_true')
  return(list(THETAest=THETAnew, hatMat=hatMat, res=res, error_matrix=errmat))
}


natspline <- function(breaks, halforder=2){
# spline(t) in span{1,t,t^2,...,t^(m-1)} + span{(t-breaks)^(m-1)} (m=halforder)
# boundary conditions s^{(m+j)}(0)=s^{(m+j)}(1)=0 for j=0,1,...,m-1
  pT <- length(breaks)
  V <- matrix(1, nrow=halforder, ncol=pT)
  if (halforder>1) for (k in 1:(halforder-1)) {
    V[k+1,] <- (1-breaks)^k
  }
  coef_mat <- rbind(diag(pT), cbind(matrix(0,halforder,halforder),
                           solve(V[,(pT-halforder+1):pT],-V[,1:(pT-halforder)])))
  nsbasis <- function(u){
    basis_mat <- matrix(NA, nrow=length(u), ncol=halforder+pT)
    for (k in 1:halforder) {
      basis_mat[,k] <- u^(k-1)
    }
    for (k in 1:pT) {
      basis_mat[,k+halforder] <- ((u>breaks[k])*(u-breaks[k]))^(2*halforder-1)
    }
    return(basis_mat%*%coef_mat)
  }
  ns_at_knots <- nsbasis(breaks)
  pen_mat <- matrix(0, nrow=halforder+pT, ncol=halforder+pT)
  for (i in 1:pT) for (j in 1:pT) {
    tmp <- 0
    for (k in 0:(halforder-1)) {
      tmp1 <- abs(breaks[i]-breaks[j])^(halforder-1-k)
      tmp2 <- (1-max(breaks[i],breaks[j]))^(halforder+k)
      tmp <- tmp + choose(halforder-1,k)*tmp1*tmp2/(halforder+k)
    }
    pen_mat[halforder+i,halforder+j] <- tmp
  }
  tmp <- (factorial(2*halforder-1)/factorial(halforder-1))^2
  pen_mat <- t(coef_mat) %*% (tmp*pen_mat) %*% coef_mat
  return(list(nsbasis=nsbasis, ns_at_knots=ns_at_knots, pen_mat=pen_mat))
}

ftreg <- function(y, Xtsr, breaks, r, M=0, rho=1e-6, halforder=2,
                  regcoef_mat=0, iter_max=80, tol=1e-8){
# y: response n-vector;
# Xtsr: covariate n*pT*p1*...*pD tensor;
# breaks: normalized time pT-vector with range (0,1);
# r: Tucker rank;
# M: covariate n*p0 matrix;
# rho: tuning parameter of roughness penalty;
# halforder: order of derivative to be regularized;
# regcoef_mat: true parameter pT*(p1*...*pD) matrix w.r.t. spline basis;
# iter_max: the maximum number of iterations;
# tol: tolerence for successful recovery /stoping criteria;
# Output estimated vectorized regression parameter, residuals, hat matrix, and running time
  t0 <- proc.time()
  p <- Xtsr@modes[-1]
  breaks.mid <- (c(0,breaks)+c(breaks,1))/2
  D <- diag(breaks.mid[-1]-breaks.mid[-(p[1]+1)])
  Ztsr <- ttm(Xtsr, D, 2)
  ns.results <- natspline(breaks, halforder)
  V <- matrix(1, nrow=p[1], ncol=halforder)
  if (halforder>1) for (k in 1:(halforder-1)) {
    V[,k+1] <- breaks^k
  }
  Amat <- ns.results$pen_mat + D%*%V%*% solve(t(V)%*%D%*%V, t(V)%*%D)
  print(paste("tuning parameter ρ =", rho))
  THETAtsr <- THETAtsr_ini <- as.tensor(0)
  if (sum(regcoef_mat^2)) {
    THETAtsr <- k_fold(ns.results$ns_at_knots%*%regcoef_mat, 1, p)
  }
  if (sum(M^2)) {
    gam_new <- solve(t(M)%*%M, t(M)%*%y)
    for (iter__ in 1:(iter_max/2)) {
      gam_old <- gam_new
      yc <- y - M %*% gam_old
      tsrreg.results <- tsrreg_RGN(yc, Ztsr, r, rho*Amat, THETAtsr, THETAtsr_ini, iter_max, tol)
      THETAtsr_ini <- tsrreg.results$THETAest
      gam_new <- gam_old + solve(t(M)%*%M, t(M)%*%tsrreg.results$res)
      if (sum((gam_new - gam_old)^2) < tol^2 * sum(gam_old^2)) {break}
    }
  } else {
    gam_old <- 0
    tsrreg.results <- tsrreg_RGN(y, Ztsr, r, rho*Amat, THETAtsr, THETAtsr_ini, iter_max, tol)
    THETAtsr_ini <- tsrreg.results$THETAest
  }
  regcoef_mat_est <- solve(t(ns.results$ns_at_knots)%*%ns.results$ns_at_knots,
                        t(ns.results$ns_at_knots)%*%k_unfold(THETAtsr_ini, 1)@data)
  return(list(res=tsrreg.results$res, hatMat=tsrreg.results$hatMat, regcoef_scal=gam_old,
              nsbasis=ns.results$nsbasis, regcoef_mat=regcoef_mat_est, time=(proc.time()-t0)[3]))
}

ftreg.gcv <- function(y, Xtsr, breaks, r, M=0, rho_list=10^(-8:-3), halforder=2,
                      regcoef_mat=0, nsbasis=NULL, iter_max=80, tol=1e-8){
# y: response n-vector;
# Xtsr: covariate n*pT*p1*...*pD tensor;
# breaks: normalized time pT-vector with range (0,1);
# r: Tucker rank;
# M: covariate n*p0 matrix;
# rho_list: tuning parameter of roughness penalty;
# halforder: order of derivative to be regularized;
# regcoef_mat: true parameter pT*(p1*...*pD) matrix w.r.t. spline basis;
# nsbasis: natural spline basis;
# iter_max: the maximum number of iterations;
# tol: tolerence for successful recovery /stoping criteria;
# Output GCV, AIC, BIC, RSS/n, and ISE if regcoef_mat and nsbasis are provided
  p <- Xtsr@modes[-1]
  breaks.mid <- (c(0,breaks)+c(breaks,1))/2
  D <- diag(breaks.mid[-1]-breaks.mid[-(p[1]+1)])
  Zmat <- k_unfold(ttm(Xtsr, D, 2), 1)@data
  n <- Xtsr@modes[1]
  trH <- RSS <- ISE <- NULL
  grid_point <- (0.5:(1e3-0.5))/1e3
  for (rho in rho_list) {
    print("  Set tuning parameter ...")
    ftreg_tmp <- ftreg(y, Xtsr, breaks, r, M, rho, halforder, regcoef_mat,
                       iter_max=iter_max, tol=tol)
    RSS <- c(RSS, sum(ftreg_tmp$res^2))
    trH <- c(trH, sum(diag(ftreg_tmp$hatMat)))
    err <- 0
    if (sum(regcoef_mat^2)) for (idx in 1:ngrid) {
      regcoef_at_grid <- c(nsbasis(grid_point[idx]) %*% regcoef_mat)
      regcoef_est_at_grid <- c(ftreg_tmp$nsbasis(grid_point[idx]) %*% ftreg_tmp$regcoef_mat)
      err <- err + sum((regcoef_at_grid - regcoef_est_at_grid)^2)/ngrid
    }
    ISE <- c(ISE, err)
  }
  GCV <- (RSS/n)/(1-trH/n)^2
  dof <- sum(r*(p-r)) + prod(r)
  AIC <- RSS + 2*dof
  BIC <- RSS + log(n)*dof
  return(list(GCV=GCV, AIC=AIC, BIC=BIC, RSSoverN=RSS/n, ISE=ISE))
}

# ftreg.cv <- function(y, Xtsr, breaks, r, M=0, rho_list=10^(-8:-3), halforder=2,
#                      cv_fold=10, iter_max=60, tol=1e-6){
#   # y: response n-vector;
#   # Xtsr: covariate n*pT*p1*...*pD tensor;
#   # breaks: normalized time pT-vector with range (0,1);
#   # r: Tucker rank;
#   # M: covariate n*p0 matrix;
#   # rho_list: tuning parameter of roughness penalty;
#   # halforder: order of derivative to be regularized;
#   # cv_fold: number of folds for cross-validation;
#   # iter_max: the maximum number of iterations;
#   # tol: tolerence for successful recovery /stoping criteria;
#   # Output CV and RSS
#   p <- Xtsr@modes[-1]
#   breaks.mid <- (c(0,breaks)+c(breaks,1))/2
#   D <- diag(breaks.mid[-1]-breaks.mid[-(p[1]+1)])
#   Zmat <- k_unfold(ttm(Xtsr, D, 2), 1)@data
#   n <- Xtsr@modes[1]
#   idx <- sample(1:n)
#   cv_list <- RSS <- NULL
#   for (rho in rho_list) {
#     print("  Set tuning parameter ...")
#     ftreg_tmp <- ftreg(y, Xtsr, breaks, r, M, rho, halforder, iter_max=1.5*iter_max, tol=tol/10)
#     RSS <- c(RSS, sum(ftreg_tmp$res^2))
#     n_sub <- n %/% cv_fold
#     cv <- tmp <- 0
#     for (k in 1:cv_fold) {
#       print(paste0('cv fold ',k,'/',cv_fold,' start'))
#       idx_k <- idx[(tmp+1):(tmp+n_sub)]
#       idx_minusk <- idx[-((tmp+1):(tmp+n_sub))]
#       Xtsr_minusk <- k_fold(k_unfold(Xtsr, 1)[idx_minusk,], 1, c(n - n_sub, p))
#       if (sum(M^2)) {
#         M_minusk <- M[idx_minusk,]
#         Mgam_k <- M[idx_k,] %*% solve(t(M_minusk)%*%M_minusk, t(M_minusk)%*%y[idx_minusk])
#       } else {
#         M_minusk <- Mgam_k <- 0
#       }
#       ftreg.results_k <- ftreg(y[idx_minusk], Xtsr_minusk, breaks, r, M_minusk,
#                                rho, halforder, iter_max=iter_max, tol=tol)
#       THETAtsr.est_vec <- c(ftreg.results_k$nsbasis(breaks)%*%ftreg.results_k$regcoef_mat)
#       cv <- cv + mean(( y[idx_k] - Mgam_k - Zmat[idx_k,] %*% THETAtsr.est_vec )^2)
#       print(paste0('cv fold ',k,'/',cv_fold,' end'))
#       tmp <- tmp + n_sub
#       if (k == cv_fold - n%%cv_fold) {
#         n_sub <- n_sub + 1
#       }
#     }
#     cv_list <- c(cv_list, cv/cv_fold)
#   }
#   o <- order(cv_list)
#   print(paste0("ρ=", rho_list[o], ", Cross-Validation=", round(cv_list[o], 14)))
#   o <- order(RSS)
#   print(paste0("ρ=", rho_list[o], ", RSS/n=", round(RSS[o]/n, 14)))
#   plot(rho_list, cv_list, log='x', xlab=expression(rho), ylab='', type='b')
#   lines(rho_list, RSS/n, type='b', lty=2, pch=2)
#   legend("topleft", c("CV","RSS/n"), lty=1:2, pch=1:2)
#   return(list(cv=cv_list, RSSoverN=RSS/n))
# }


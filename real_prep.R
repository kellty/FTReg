library(readr)
pheno <- read_csv("ADHD/allSubs_testSet_phenotypic_dx.csv")
pheno <- pheno[(pheno$Site==1)&(pheno$`ADHD Index`!='-999'),]
n_ADHD <- nrow(pheno)
y_ADHD <- as.numeric(pheno$`ADHD Index`) # Inattentive + Hyper/Impulsive
M_ADHD <- cbind(rep(1, n_ADHD), pheno$Gender, pheno$Age, as.numeric(pheno$`Med Status`),
               as.numeric(pheno$`Verbal IQ`), as.numeric(pheno$`Performance IQ`))
library(rTensor)
pT <- 50
Xtsr_ADHD <- k_fold(matrix(nrow=n_ADHD, ncol=64*64*33*pT), 1, c(n_ADHD,64,64,33,pT))
library(RNifti)
for (i in 1:n_ADHD) {
  tmp <- readNifti(paste0("ADHD/Peking_",as.character(pheno$ID[i]),
                          "_1/rest_1/NIfTI/rest.nii.gz"))
  Xtsr_ADHD[i,,,,] <- tmp[,,,101:(100+pT)]
}
ave_mat <- function(ntot, nsub){
  H <- matrix(0, nrow=ntot/nsub, ncol=ntot)
  for (i in 1:nrow(H)) {
    H[i,((i-1)*nsub+1):(i*nsub)] <- 1/nsub
  }
  return(H)
}
ave_z <- rbind(cbind(ave_mat(24,8),matrix(0,3,9)), c(rep(0,24),rep(1/9,9))) # ave_mat(33,4)
ave_mat_list <- list(ave_mat(64,8),ave_mat(64,8),ave_z)
Xtsr_ADHD <- ttl(Xtsr_ADHD, ave_mat_list, 2:4)
Xtsr_ADHD <- k_fold(k_unfold(Xtsr_ADHD, 5), 2, c(n_ADHD,pT,8,8,4))
breaks_ADHD <- (0.5:(pT-0.5))/pT
save(y_ADHD,M_ADHD,Xtsr_ADHD,breaks_ADHD, file="ADHD_prep.RData")

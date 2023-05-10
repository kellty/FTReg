# Implementation of Functional Tensor Regression

## Data
ADHD data from [http://fcon_1000.projects.nitrc.org/indi/adhd200/](http://fcon_1000.projects.nitrc.org/indi/adhd200/) was used for this work. The file "***ADHD_prep.RData***" was obtained after pre-processing by "***real_prep.R***".

## Code
- "***func.R***": Build functions to implement the Riemannian Gauss--Newton algorithm for functional tensor regression (FTReg).
- "***simu_begin.R***": Generate the simulation setting.
- "***simu_rho.R***": Perform FTReg on simulated data under different tuning parameters, with coefficient rank and sample size fixed.
- "***simu_r.R***": Perform FTReg on simulated data under different coefficient ranks, with tuning parameter selected to be optimal and sample size fixed.
- "***simu_n.R***": Perform FTReg on simulated data under different sample sizes, with tuning parameter selected to be optimal and coefficient rank fixed.
- "***real.R***": Perform FTReg on ADHD data, with coefficient rank and tuning parameter selected to be optimal.

## Workflow
First of all, make sure that "***func.R***" is in the working directory.
- Simulation: Run "***simu_begin.R***" and *then* any of "***simu_rho.R***", "***simu_r.R***" and "***simu_n.R***".
- Real data example: Make sure that "***ADHD_prep.RData***" is in the working directory. Run "***real.R***".

#' Get MLEs for coefficients and variance
#'
#' Calculate log likelihood for given parameters
#'
#' @param param a 4x1 vector of population mean, genetic effect, 
#' genetic variance and environmental variance
#' @param Phi a NxN kinship matrix for all phase 1 individuals
#' @param y a Nx1 vector of phenotype for all phase 1 individuals
#' @param X a Nx2 matrix, where the first column is 1 and the second 
#' column is the genotypes (0, 1 or 2) for all phase 1 individuals
#' @param pwt a nxn matrix of pairwise sampling probability of 
#' phase 2 individuals only 
#' @param R a Nx1 vector of sampling indicator (1 if sampled and 
#' 0 otherwise) for all phase 1 individuals
#'
#' @export
#' @return The deviance (-2*log-likelihood)
#'
#' @examples
#' data(sim)
#' loglik <- calcLL(c(50,-0.08,4.3,0.56), Phi, y, X, pwt, R)
calcLL <-
  function(param, Phi, y, X, pwt, R)
  {
    
    n <- nrow(Phi)
    if(!is.matrix(X)) X <- as.matrix(X)
    stopifnot(nrow(X) == n)
    if(!is.matrix(y)) y <- as.matrix(y)
    stopifnot(nrow(y) == n)
    
    # Variance-covariance matrix
    C <- param[3]*Phi+param[4]*diag(n)
    
    sub <- which(R==1)
    y_sub <- y[sub]
    X_sub <- as.matrix(X[sub,])
    C_sub_inv <- solve(C)[sub,sub]
    
    e <-  y_sub-X_sub%*%param[1:2]
    LL <- -0.5*(determinant(C, logarithm = TRUE)$modulus + t(e)%*%((1/pwt)*C_sub_inv)%*%e)
    
    if (!is.finite(LL)) stop("Oh NO!!!")
    
    -2*LL
    
  }

#' Fit a linear mixed model
#'
#' Fit a linear mixed model of the form y = Xb + e where e follows a
#' multivariate normal distribution with mean 0 and variance matrix
#' `sigma_g^2*Phi+sigma_e^2*I`, where `sigma_g^2` is the genetic
#' variance, `sigma_e^2` is the environmental variance,`Phi` is a 
#' known kniship matrix and `I` is the identity matrix.
#'
#' @param Phi a NxN kinship matrix for all phase 1 individuals
#' @param y a Nx1 vector of phenotype for all phase 1 individuals
#' @param X a Nx2 matrix, where the first column is 1 and the second 
#' column is the genotypes (0, 1 or 2) for all phase 1 individuals
#' @param pwt a nxn matrix of pairwise sampling probability of 
#' phase 2 individuals only 
#' @param R a Nx1 vector of sampling indicator (1 if sampled and 
#' 0 otherwise) for all phase 1 individuals
#' @param tol Tolerance for convergence
#'
#' @export
#' @return a 4x1 vector of population mean, genetic effect, 
#' genetic variance and environmental variance. Heritability is 
#' defined by genetic variance divided by the the sum of genetic 
#' and environmental variance
#'
#' @examples Fit LMM under outcome-dependent sampling
#' data(sim)
#' extreme <- which(y<=quantile(y,0.15)|y>=quantile(y,0.85))
#' middle <- which(y>quantile(y,0.15)&y<quantile(y,0.85))
#' N <- length(y)
#' n <- N/2
#' na <- length(extreme)
#' nb <- length(middle)
#' nn <- n-na
#' pwt11 <- matrix(1,na,na)
#' pwt12 <- matrix(nn/nb,nrow=na,ncol=nb)
#' pwt22 <- matrix((nn/nb)^2,nb,nb); diag(pwt22) <- nn/nb
#' pwt <- rbind(cbind(pwt11,pwt12),cbind(t(pwt12),pwt22))
#' ord <- order(c(extreme,middle))
#' pwt <- pwt[ord,ord]
#' sub <- sort(c(extreme,sample(middle,nn,replace = FALSE)))
#' R <- numeric(N)
#' R[sub] <- 1
#' result <- fitLMM(Phi, y, X, pwt[sub,sub], R)
#'
require(lme4qtl)
require(minqa)
fitLMM <-
  function(Phi, y, X, pwt, R, thres=1e-8)
  {
    n <- nrow(Phi)
    if(!is.matrix(X)) X <- as.matrix(X)
    stopifnot(nrow(X) == n)
    if(!is.matrix(y)) y <- as.matrix(y)
    stopifnot(nrow(y) == n)
    
    # Get initial parameters using unweighted lme4qtl
    sub <- which(R==1)
    df <- data.frame(ID=rownames(Phi[sub,sub]), trait=y[sub], dosage=X[sub,2])
    m1 <- lme4qtl::relmatLmer(trait ~ dosage + (1|ID), df, relmat = list(ID = Phi[sub,sub]), REML=FALSE)
    v1 <- lme4qtl::VarProp(m1)
    param <- c(m1@beta,v1$vcov)
    if (param[3]<thres) param[3]<-thres
    if (param[4]<thres) param[4]<-thres
    
    lower <- c(-Inf,-Inf,thres,thres)
    upper <- c(Inf,Inf,Inf,Inf)
    # maximize log likelihood
    out <- minqa::bobyqa(par=param, fn=calcLL, lower = lower, upper = upper,Phi=Phi, y=y, X=X, pwt=pwt, R=R)
    out
  }
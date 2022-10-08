#' Get MLEs for coefficients and variance
#'
#' Calculate log likelihood for given parameters
#'
#' @param param a 4x1 vector of population mean, genetic effect, 
#' total phenotypic variance and heritability
#' @param Phi a NxN kinship matrix for all phase 1 individuals
#' @param y a Nx1 vector of phenotype for all phase 1 individuals
#' @param X a Nx2 matrix, where the first column is 1 and the second 
#' column is the genotypes (0, 1 or 2) for all phase 1 individuals
#' @param p a nx1 vector of first-order sampling probability of 
#' phase 2 individuals only 
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
#' loglik <- calcLL(c(50,-0.08,4.3,0.56), Phi, y, X, p, pwt, R)
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
#' `sigma^2*(h^2*Phi+(1-h)*I)`, where `sigma^2` is the total phenotypic
#' variance, `h^2` is the heritability,`Phi` is a known kniship
#' matrix and `I` is the identity matrix.
#'
#' @param Phi a NxN kinship matrix for all phase 1 individuals
#' @param y a Nx1 vector of phenotype for all phase 1 individuals
#' @param X a Nx2 matrix, where the first column is 1 and the second 
#' column is the genotypes (0, 1 or 2) for all phase 1 individuals
#' @param p a nx1 vector of first-order sampling probability of 
#' phase 2 individuals only 
#' @param pwt a nxn matrix of pairwise sampling probability of 
#' phase 2 individuals only 
#' @param R a Nx1 vector of sampling indicator (1 if sampled and 
#' 0 otherwise) for all phase 1 individuals
#' @param tol Tolerance for convergence
#'
#' @export
#' @return a 4x1 vector of population mean, genetic effect, 
#' total phenotypic variance and heritability
#'
#' @examples
#' data(sim)
#' result <- fitLMM(Phi, y, X, p, pwt, R)
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
    
    # Get initial value for beta, sigmasq and heritability using unweighted lme4qtl
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
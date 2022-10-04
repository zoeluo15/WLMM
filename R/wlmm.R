#' Get MLEs for coefficients and variance
#'
#' For a fixed value for `hsq`, the heritability, calculate the
#' corresponding maximum likelihood estimates of `beta` and
#' `sigmasq`, with the latter being the total variance,
#' `sigmasq_g + sigmasq_e`.
#'
#' @param hsq heritability
#' @param Phi relatedness matrix 
#' @param y phenotypes 
#' @param X covariate matrix 
#' @param p first-order sampling probability for all individuals
#' @param pwt pairwise sampling probability for sampled individuals
#' @param R sampling indicator (1 if sampled and 0 otherwise) for all individuals
#'
#' @export
#' @return list containing `beta` and `sigmasq`, with residual
#' sum of squares as attributes.
#'
#' @examples
#' data(sim)
#' ml <- getMLsoln(0.5, Phi, y, X, p, pwt, R)
getMLsoln <-
  function(hsq, Phi, y, X, p, pwt, R)
  {
    N <- nrow(Phi)
    if(!is.matrix(X)) X <- as.matrix(X)
    stopifnot(nrow(X) == N)
    if(!is.matrix(y)) y <- as.matrix(y)
    stopifnot(nrow(y) == N)
    
    # Covariance matrix devide by sigma^2
    C <- hsq*Phi+(1-hsq)*diag(N)
    invsqrtC <- chol(solve(C))
    
    sub <- which(R==1)
    n <- length(sub)
    ij<-as.matrix(subset(expand.grid(i=1:n,j=1:n),i<j))
    ii<-ij[,1]
    jj<-ij[,2]
    
    invsqrtC_sub <- invsqrtC[sub,sub]
    
    X_sub <- as.matrix(X[sub,,drop=FALSE])
    y_sub <- y[sub]
    
    dKsi <- invsqrtC_sub[ii,]-invsqrtC_sub[jj,]
    t1 <- 2*sum(1/pwt[ij]*(dKsi%*%y_sub)*(dKsi%*%X_sub[,1]))
    t2 <- 2*sum(1/pwt[ij]*(dKsi%*%y_sub)*(dKsi%*%X_sub[,2]))
    t3 <- 2*sum(1/pwt[ij]*(dKsi%*%X_sub[,1])*(dKsi%*%X_sub[,2]))
    t4 <- 2*sum(1/pwt[ij]*(dKsi%*%X_sub[,1])^2)
    t5 <- 2*sum(1/pwt[ij]*(dKsi%*%X_sub[,2])^2)
    
    beta <- solve(matrix(data=c(t4,t3,t3,t5), nrow=2), c(t1,t2))
    beta[1] <- (t(y - beta[2]*X[,2])%*%(1/p))/sum(1/p)
    
    # calculate a bunch of matrices and RSS
    e <- invsqrtC_sub%*%(y[sub]-as.matrix(X[sub,])%*%beta)
    rss <- crossprod(e[ii]-e[jj],1/pwt[ij]*(e[ii]-e[jj]))/(N)
    
    sigmasq <- rss / N 
    
    # return value
    result <- list(beta=beta, sigmasq=sigmasq)
    attr(result, "rss") <- rss
    
    result
  }

#' Calculate log likelihood for a given heritability
#'
#' Calculate the log likelihood for a given value of the heritability, `hsq`.
#'
#' @param hsq heritability
#' @param Phi relatedness matrix 
#' @param y phenotypes 
#' @param X covariate matrix 
#' @param p first-order sampling probability for all individuals
#' @param pwt pairwise sampling probability for sampled individuals
#' @param R sampling indicator (1 if sampled and 0 otherwise) for all individuals
#'
#' @export
#' @return The log likelihood value, with the corresponding estimates
#' of `beta` and `sigmasq` included as attributes.
#'
#' @examples
#' data(sim)
#' loglik <- calcLL(0.5, Phi, y, X, p, pwt, R)
#' many_loglik <- calcLL(seq(0, 1, by=0.1), Phi, y, X, p, pwt, R)
calcLL <-
  function(hsq, Phi, y, X, p, pwt, R)
  {
    if(length(hsq) > 1)
      return(vapply(hsq, calcLL, 0, Phi, y, X, p, pwt, R))
    
    N <- nrow(Phi)
    
    # estimate beta and sigmasq
    MLsoln <- getMLsoln(hsq, Phi, y, X, p, pwt, R)
    beta <- MLsoln$beta
    sigmasq <- MLsoln$sigmasq
    
    # calculate log likelihood
    rss <- attr(MLsoln, "rss")
    
    LL <- -0.5*(log(det((hsq*Phi+(1-hsq)*diag(N)))) + N*log(rss))
    
    attr(LL, "beta") <- beta
    attr(LL, "sigmasq") <- sigmasq
    LL
  }

#' Fit a linear mixed model
#'
#' Fit a linear mixed model of the form y = Xb + e where e follows a
#' multivariate normal distribution with mean 0 and variance matrix
#' `sigmasq_g K + sigmasq_e I`, where `K` is a known kinship
#' matrix and `I` is the identity matrix.
#'
#' @param Phi relatedness matrix 
#' @param y phenotypes 
#' @param X covariate matrix 
#' @param p first-order sampling probability for all individuals
#' @param pwt pairwise sampling probability for sampled individuals
#' @param R sampling indicator (1 if sampled and 0 otherwise) for all individuals
#' @param check_boundary If TRUE, explicitly check log likelihood at 0 and 1.
#' @param tol Tolerance for convergence
#' @param compute_se if TRUE, return the standard error of the `hsq`
#' estimate using the Fisher Information matrix of the MLE estimate. The
#' standard error will be in an `attr` of  `hsq` in the output.
#' Currently requires `use_cpp = FALSE`, and so if `compute_se=TRUE`
#' we take `use_cpp=FALSE`.
#'
#' @importFrom stats optim
#'
#' @export
#' @return List containing estimates of `beta`, `sigmasq`,
#' `hsq`, `sigmasq_g`, and `sigmasq_e`, as well as the log
#' likelihood (`loglik`). If `compute_se=TRUE`, the output also
#' contains `hsq_se`.
#'
#' @examples
#' data(sim)
#' result <- fitLMM(Phi, y, X, p, pwt, R)
#'
#' # also compute SE
#' wSE <- fitLMM(Phi, y, X, p, pwt, R, compute_se = TRUE, use_cpp=FALSE)
#' c(hsq=wSE$hsq, SE=wSE$hsq_se)
#'
fitLMM <-
  function(Phi, y, X, p, pwt, R, check_boundary=TRUE, tol=1e-4, compute_se = FALSE)
  {
    N <- nrow(Phi)
    if(!is.matrix(X)) X <- as.matrix(X)
    stopifnot(nrow(X) == N)
    if(!is.matrix(y)) y <- as.matrix(y)
    stopifnot(nrow(y) == N)
    
    # maximize log likelihood
    out <- stats::optimize(calcLL, c(0, 1), Phi=Phi, y=y, X=X, p=p, pwt=pwt, R=R, maximum=TRUE, tol=tol)
    
    # Use the hessian to get the stanard errors; had to use `optim` here...
    if(compute_se){
      calcLL <- function(hsq, Phi, y, X, p, pwt, R){
        return(-1*calcLL(hsq, Phi, y, X, p, pwt, R))
      }
      opt2 <- optim(out$maximum, calcLL, Phi=Phi, y=y, X=X, p=p, pwt=pwt, R=R,
                    hessian=TRUE, lower = 0, upper = 1, method = "Brent")
      vc <- solve(opt2$hessian) # var-cov matrix
      se <- sqrt(diag(vc))        # standard errors
    }
    
    hsq <- out$maximum
    obj <- out$objective
    sigmasq <- attr(obj, "sigmasq")
    
    result <- list(beta=attr(obj, "beta"),
                   sigmasq=sigmasq, # total var
                   hsq=hsq,
                   sigmasq_g= hsq*sigmasq, # genetic variance
                   sigmasq_e = (1-hsq)*sigmasq, # residual variance
                   loglik = as.numeric(obj))
    if(compute_se) result$hsq_se <- se
    
    result
  }

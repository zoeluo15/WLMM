# Weighted maximumn likelihood estimation via Monte Carlo EM algorithm for generalized linear mixed models

rlatent<-function(mu, sigma2, observed){
  lower<-numeric(length(mu))
  lower[observed==0]<- -Inf
  upper<-numeric(length(mu))
  upper[observed==1]<- Inf
  truncnorm::qtruncnorm(runif(length(mu)), lower, upper, mu, sqrt(sigma2))
}

library(Matrix)
tr<-function(m) sum(diag(m))

sim_latent<-function(observed, latent_old, Xi_inv, xbeta){
  n<-length(observed)
  cond_var<-1/diag(Xi_inv)
  out <- numeric(n)
  for (i in 1:n){
    mu<-latent_old[i]-cond_var[i]*sum(Xi_inv[i,]*(latent_old-xbeta))  ##https://arxiv.org/pdf/1810.10559.pdf
    latent_old[i] <- out[i] <- rlatent(mu,cond_var[i],observed[i])
  }
  out
}

gibbs<-function(observed, Xi_inv, xbeta, LOTS=3000, BURNIN=100){
  n<-length(observed)
  latent_old<-xbeta #this is ystar
  sums<-numeric(n)
  sumsq<-matrix(0,n,n)
  for(t in 1:LOTS){
    latent_old<-sim_latent(observed,latent_old,Xi_inv,xbeta)
    if(t>BURNIN){
      sums<-sums+latent_old
      means<-sums/(t-BURNIN)
      sumsq<-sumsq+tcrossprod(latent_old-means)
    }
  }
  list(tmean=sums/(LOTS-BURNIN),tvar= sumsq/(LOTS-BURNIN))
}

condExpLL<-function(h, X, kinship, beta, G, sub, pwt){
  n <- nrow(kinship)
  Xi <- (h/(1-h))*kinship+Diagonal(n)
  Xi_sub_inv <- solve(Xi)[sub,sub]
  LL <- -0.5*as.numeric(determinant(Xi, logarithm = TRUE)$modulus +
                          tr((1/pwt*Xi_sub_inv)%*%G$var) +
                          crossprod(G$mean-X%*%as.matrix(beta),1/pwt*Xi_sub_inv)%*%(G$mean-X%*%as.matrix(beta)))
  LL
}

update_param<-function(X, y, FID, kinship, beta_old, h_old, sub, pwt, method, tol=1e-15){
  N <- nrow(kinship)
  Xi <- (h_old/(1-h_old))*kinship+Diagonal(N)
  Xi_sub <- Xi[sub,sub]
  Xi_sub_inv <- solve(Xi_sub)
  A<-list(); B<-list()
  fam<-unique(FID)
  for (i in 1:length(fam)) {
    idx<-which(FID==fam[i])
    X_f <- X[idx,]
    y_f <- y[idx]
    if (method=="Gibbs") {
      MV<-gibbs(y_f, Xi_sub_inv[idx,idx], X_f%*%as.matrix(beta_old))
    } else if (method=="Numerical"){ # Numerical integration
      nn<-length(idx)
      lower<-rep(0,nn); lower[y_f==0]<- -Inf
      upper<-rep(0,nn); upper[y_f==1]<- Inf
      MV<-tmvtnorm::mtmvnorm(mean=as.vector(X_f%*%as.matrix(beta_old)),sigma=Xi_sub[idx,idx],lower=lower,upper=upper,doComputeVariance=TRUE)
    } else {
      stop("Undefined method.")
    }
    A[[i]]<-MV$tmean
    B[[i]]<-MV$tvar
  }
  G<-list(mean=unlist(A),var=bdiag(B))
  out <- stats::optimize(condExpLL,c(0,1),X=X,kinship=kinship,beta=beta_old,G=G,sub=sub,pwt=pwt,maximum=TRUE,tol=tol)
  h_new <- out$maximum
  Xi <- (h_new/(1-h_new))*kinship+Diagonal(N)
  Xi_sub_inv <- solve(Xi)[sub,sub]
  XwV <- crossprod(X,1/pwt*Xi_sub_inv)
  beta_new <- solve(XwV%*%X,XwV%*%G$mean)
  list(h=h_new,beta=beta_new)
}

probit_em<-function(data, kinship, beta_init=NULL, h_init, R=NULL, pwt=NULL, maxit=150, method="Gibbs", VERBOSE=FALSE){
  N <- nrow(data)
  if (is.null(R)) R <- rep(1,N)
  if (is.null(pwt)) pwt <- matrix(1,N,N)
  sub <- which(R==1)
  X <- cbind(rep(1,N),data$genotype); X_sub <- as.matrix(X[sub,])
  y <- data$phenotype; y_sub <- y[sub]
  FID <- data$FID; FID_sub <- FID[sub]
  if (is.null(beta_init)) beta_init <- glm(y_sub~X_sub-1,family=binomial(link=probit))$coef
  new <- update_param(X=X_sub,y=y_sub,FID=FID_sub,kinship=kinship,beta_old=beta_init,h_old=h_init,sub=sub,pwt=pwt,method=method)
  beta <- new$beta[,1]
  h <- new$h
  sigmag <- sqrt(h/(1-h))
  out<-c(beta,sigmag,h)

  if (VERBOSE) print(data.frame(iter=1,beta0=out[1],beta1=out[2],sigma=out[3],h2=out[4]))
  
  for (m in 1:maxit) {
    new <- update_param(X=X_sub,y=y_sub,FID=FID_sub,kinship=kinship,beta_old=beta,h_old=h,sub=sub,pwt=pwt,method=method)
    beta <- new$beta[,1]
    h <- new$h
    sigmag <- sqrt(h/(1-h))
    out <- rbind(out,c(beta,sigmag,h))
    if (VERBOSE) print(data.frame(iter=m+1,beta0=out[m+1,1],beta1=out[m+1,2],sigma=out[m+1,3],h2=out[m+1,4]))
  }
  out
}

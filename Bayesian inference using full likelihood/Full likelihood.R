imputation <- nimbleFunction(
  run = function(geno = double(1), FID = double(2), indicator = double(1), 
                 rowx = double(), colx = double(), maf = double()) {
    returnType(double(1))
    p <- numeric(length = 3)
    if (colx==3|colx==4) {
      if (indicator[1]==1&indicator[2]==0) { # genotype of the first child is observed
        p[1] <- ((2-(geno[FID[rowx,1]]-1))/2)*(1-maf)
        p[2] <- ((2-(geno[FID[rowx,1]]-1))/2)*maf+((geno[FID[rowx,1]]-1)/2)*(1-maf)
        p[3] <- ((geno[FID[rowx,1]]-1)/2)*maf
        return(p)
      } else if (indicator[1]==0&indicator[2]==1) { # genotype of the second child is observed
        p[1] <- ((2-(geno[FID[rowx,2]]-1))/2)*(1-maf)
        p[2] <- ((2-(geno[FID[rowx,2]]-1))/2)*maf+((geno[FID[rowx,2]]-1)/2)*(1-maf)
        p[3] <- ((geno[FID[rowx,2]]-1)/2)*maf
        return(p)
      } else if (indicator[1]==1&indicator[2]==1) { # genotypes of both children are observed
        p[1] <- 0.5*((2-(geno[FID[rowx,1]]-1))/2)*(1-maf)+0.5*((2-(geno[FID[rowx,2]]-1))/2)*(1-maf)
        p[2] <- 0.5*(((2-(geno[FID[rowx,1]]-1))/2)*maf+((geno[FID[rowx,1]]-1)/2)*(1-maf))+
          0.5*(((2-(geno[FID[rowx,2]]-1))/2)*maf+((geno[FID[rowx,2]]-1)/2)*(1-maf))
        p[3] <- 0.5*((geno[FID[rowx,1]]-1)/2)*maf+0.5*((geno[FID[rowx,2]]-1)/2)*maf
        return(p)
      } else { # genotypes of both children are missing
        p[1] <- (1-maf)^2
        p[2] <- 2*maf*(1-maf)
        p[3] <- maf^2
        return(p)
      }
    } else {
      if (indicator[3]==1&indicator[4]==0) { # genotype of the father is observed
        p[1] <- ((2-(geno[FID[rowx,3]]-1))/2)*(1-maf)
        p[2] <- ((2-(geno[FID[rowx,3]]-1))/2)*maf+((geno[FID[rowx,3]]-1)/2)*(1-maf)
        p[3] <- ((geno[FID[rowx,3]]-1)/2)*maf
        return(p)
      } else if (indicator[3]==0&indicator[4]==1) { # genotype of the mother is observed
        p[1] <- ((2-(geno[FID[rowx,4]]-1))/2)*(1-maf)
        p[2] <- ((2-(geno[FID[rowx,4]]-1))/2)*maf+((geno[FID[rowx,4]]-1)/2)*(1-maf)
        p[3] <- ((geno[FID[rowx,4]]-1)/2)*maf
        return(p)
      } else if (indicator[3]==1&indicator[4]==1) { # genotypes of both parents are observed
        p[1] <- ((2-(geno[FID[rowx,3]]-1))/2)*((2-(geno[FID[rowx,4]]-1))/2)
        p[2] <- ((2-(geno[FID[rowx,3]]-1))/2)*((geno[FID[rowx,4]]-1)/2)+((geno[FID[rowx,3]]-1)/2)*((2-(geno[FID[rowx,4]]-1))/2)
        p[3] <- ((geno[FID[rowx,3]]-1)/2)*((geno[FID[rowx,4]]-1)/2)
        return(p)
      } else if ((colx==1&indicator[2]==1&indicator[3]==0&indicator[4]==0)|
                 (colx==2&indicator[1]==1&indicator[3]==0&indicator[4]==0)) { 
        # genotypes of both parents are missing
        # genotype of the sibling is observed
        p[1] <- 0.25
        p[2] <- 0.5
        p[3] <- 0.25
        return(p)
      } else { # genotypes of the all relatives are missing
        p[1] <- (1-maf)^2
        p[2] <- 2*maf*(1-maf)
        p[3] <- maf^2
        return(p)
      }
    }
  })

## define the model
code <- nimbleCode({
  beta_0 ~ dunif(-1000, 1000)
  beta_1 ~ dunif(-1000, 1000)
  sigma_g2 ~ dunif(0, 1000)
  sigma_e2 ~ dunif(0, 1000)
  # sampling
  for (i in 1:n_middle) {
    R[middle[i]] ~ dbinom(size=1, prob=frac) # sampling indicator of the non-extreme individuals
  }
  # impute missing genotype
  for (i in 1:n_rest) {
    indicator[i,1:4] <- c(R[FID[row_idx[i],1]],R[FID[row_idx[i],2]],R[FID[row_idx[i],3]],R[FID[row_idx[i],4]])
    out[i,1:3] <- imputation(geno=x_raw[],FID=FID[,],indicator=indicator[i,1:4],
                             rowx=row_idx[i],colx=col_idx[i],maf=maf)
    x[rest[i]] ~ dcat(prob=out[i,1:3])
  }
  # construct the covariance matrix
  for (j in 1:4) {
    for (k in 1:4) {
      Xi[j,k] <- sigma_g2*kinship[j,k]+sigma_e2*(j==k)
    }
  }
  # modeling
  for (i in 1:nf) {
    for (j in 1:4) {
      mu[i,j] <- beta_0 + beta_1*x[FID[i,j]]
    }
    y[i,1:4] ~ dmnorm(mu[i,1:4],cov=Xi[1:4,1:4])
  }
})

load("sim04(N1200)newb1.RData")
# Sample 50% the population 
extreme <- which(sim$phenotype<=quantile(sim$phenotype,0.166)|sim$phenotype>=quantile(sim$phenotype,0.834))
middle <- which(sim$phenotype>quantile(sim$phenotype,0.166)&sim$phenotype<quantile(sim$phenotype,0.834))

N <- nrow(sim)
FID <- matrix(1:N,ncol = 4); nf <- nrow(FID)
# Reorder y where rows correspond to families and cols correspond to family members
y <- matrix(NA,nf,4)
for (i in 1:nf) {
  y[i,] <- sim$phenotype[FID[i,]]
}
# rel matrix is the same for all families, so only keep one block
rel <- rel[FID[1,],FID[1,]]

# generate sample
sub <- sort(c(extreme,sample(middle,N*(0.5)-length(extreme))))
R <- numeric(N); R[sub] <- 1
# set x to NA for individuals that are not sampled
x_raw <- x <- sim$genotype+1
x[which(R==0)] <- NA

sample <- which(R==1); n_sample <- length(sample)
rest <- which(R==0); n_rest <- length(rest)
col_idx <- floor((rest-1)/nf)+1
row_idx <- rest-nf*(col_idx-1)
data <- list(y=y,x=x,R=R)
constants <- list(x_raw=x_raw,kinship=rel,middle=middle,sample=sample,rest=rest,FID=FID,
                  maf=0.2,nf=nf,frac=(N*(0.5)-length(extreme))/length(middle),n_middle=length(middle),
                  n_sample=n_sample,n_rest=n_rest,col_idx=col_idx,row_idx=row_idx)
inits <- list(beta_0=50, beta_1=0, sigma_g2=2, sigma_e2=2,x=sim$genotype+1)


lmmModel <- nimbleModel(code = code, constants = constants, data = data, inits = inits, check = FALSE)
spec <- configureMCMC(lmmModel, monitors = c("beta_0", "beta_1", "sigma_g2", "sigma_e2"), monitors2 = "x")
customLmmMCMC <- buildMCMC(spec)
ClmmModel <- compileNimble(lmmModel)
ClmmMCMC <- compileNimble(customLmmMCMC, project = lmmModel)
samples <- runMCMC(ClmmMCMC, niter = 10000)

param <- samples$samples[200:nrow(samples$samples),]
param[,1] <- param[,1]+param[,2] # Transformation of beta0
full <- colMeans(param)
est <- c(full[1],full[2],full[3]+full[4],full[4]/(full[3]+full[4]))
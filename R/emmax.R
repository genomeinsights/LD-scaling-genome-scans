library(data.table)
emmax <- function(Y,X,K,B=NULL,binary=FALSE,Covar=NULL,cores=1) {
  
  #The SNP based single-locus association test applied to wild data which accounts for relatedness structure.
  #Reference: Li Z, Kemppainen P, Rastas P, Merila, J. Linkage disequilibrium clustering-based approach for association mapping with tightly linked genome-wide data
  #Contact: Zitong Li, Email: lizitong1985@gmail.com
  #Input:
  #-Y is a vector consisting of the phenotype data
  #-X is a matrix of the genotype data, with row represents individuals and #column represents SNPs
  #-K is the pre-calculated kinship matrix (for example can be calculated by the 'A.mat' function in R package 'rrBLUP')
  #-B is the number of permutations 
  #Output:
  #-F: a vector of F statistic of the SNPs
  #-pval: a vector of original P-values of the SNPs
  #-pval.corr: the P-values adjusted by the permutation test of the SNPs
  
  
  n<-length(Y)
  m<-ncol(X)
  nbchunks = 2
  stopifnot(ncol(K) == n)
  stopifnot(nrow(K) == n)
  stopifnot(nrow(X) == n)
  stopifnot(nbchunks >= 2)
  
  #INTERCEPT
  
  if(!is.null(Covar) & length(unique(Covar))>1){
    Y <- resid(lm(Y~Covar))
  }   
  
  
  Xo<-rep(1,n)
  
  #K MATRIX NORMALISATION
  
  K_norm<-(n-1)/sum((diag(n)-matrix(1,n,n)/n)*K)*K
  
  #NULL MODEL
  
  null<-emma.REMLE(Y,as.matrix(Xo),K_norm)
  
  pseudoh<-null$vg/(null$vg+null$ve)
  
  
  ############################################################
  #EMMAX: to conduct association test and calculate p-values
  
  
  
  null<-emma.REMLE(Y,as.matrix(Xo),K_norm)
  
  M<-solve(chol(null$vg*K_norm+null$ve*diag(n)))
  Y_t<-crossprod(M,Y)
  Xo_t<-crossprod(M,Xo)
  
  RSS<-list()
  for (j in 1:(nbchunks-1)) {
    X_t<-crossprod(M,X[,((j-1)*round(m/nbchunks)+1):(j*round(m/nbchunks))])
    RSS[[j]]<-apply(X_t,2,function(x){sum(lsfit(cbind(Xo_t,x),Y_t,intercept = FALSE)$residuals^2)})
    rm(X_t)}
  X_t<-crossprod(M,X[,((j)*round(m/nbchunks)+1):(m)])
  RSS[[nbchunks]]<-apply(X_t,2,function(x){sum(lsfit(cbind(Xo_t,x),Y_t,intercept = FALSE)$residuals^2)})
  rm(X_t,j)
  
  RSSf<-unlist(RSS)
  RSS_H0<-rep(sum(lsfit(Xo_t,Y_t,intercept = FALSE)$residuals^2),m)
  df1<-1
  df2<-n-df1-1
  R2<-1-1/(RSS_H0/RSSf)
  F<-(RSS_H0/RSSf-1)*df2/df1
  pval<-pf(F,df1,df2,lower.tail=FALSE)
  
  
  ######################################################
  # B <- 10 n length(Covar)
  if(!is.null(B)){
    Fo <- sort(abs(F))
    
    Ord <- order(F)
    
    sigma <- K_norm*null$vg+diag(rep(1,n))*null$ve
    
    Ynew <- mvrnorm(B, rep(0, ncol(sigma)), sigma)
    
    if(binary){
      Min <- min(table(Y))
      
      Ynew <-  apply(Ynew, 1, function(x){
        Ord <- order(x)
        y <- rep(NA, length(Ord))
        y[Ord<=Min] <- 1
        y[Ord>Min] <- 0
        if(!is.null(Covar) & length(unique(Covar))>1){
          Y <- resid(lm(Y~Covar))
        }   
        
        return(y)
      })
      Ynew <- t(Ynew)
      
    }
    #FR <- matrix(0,nrow=B,ncol=m)
    #k <- 1
    FR <- do.call(rbind,mclapply(1:B, function(k){
      
      y <- Ynew[k,]
      
      null<-emma.REMLE(y,as.matrix(Xo),K_norm)
      
      
      Y_t<-crossprod(M,y)
      
      M<-solve(chol(null$vg*K_norm+null$ve*diag(n)))
      Y_t<-crossprod(M,y)
      Xo_t<-crossprod(M,Xo)
      
      RSS<-list()
      for (j in 1:(nbchunks-1)) {
        X_t<-crossprod(M,X[,((j-1)*round(m/nbchunks)+1):(j*round(m/nbchunks))])
        RSS[[j]]<-apply(X_t,2,function(x){sum(lsfit(cbind(Xo_t,x),Y_t,intercept = FALSE)$residuals^2)})
        rm(X_t)}
      X_t<-crossprod(M,X[,((j)*round(m/nbchunks)+1):(m)])
      RSS[[nbchunks]]<-apply(X_t,2,function(x){sum(lsfit(cbind(Xo_t,x),Y_t,intercept = FALSE)$residuals^2)})
      rm(X_t,j)
      
      RSSf<-unlist(RSS)
      RSS_H0<-rep(sum(lsfit(Xo_t,Y_t,intercept = FALSE)$residuals^2),m)
      df1<-1
      df2<-n-df1-1
      #R2<-1-1/(RSS_H0/RSSf)
      (RSS_H0/RSSf-1)*df2/df1
      
      
    },mc.cores = cores))
    
    Qmat <- t(apply(FR,1,cummax))
    
    Padj <- apply(t(matrix(rep(Fo,B),m)) < Qmat, 2, mean)
    
    o <- order(Ord)
    
    out <- Padj[o]
    names(out) <- names(pval)
    list('F'=F,'pval'=pval,'pval.corr'=out,'Rsq'=R2) 
  }else{
    list('F'=F,'pval'=pval,'Rsq'=R2) 
  }
  #Conduct permutation test
  
  
}

emma.REMLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10, esp=1e-10, eig.L = NULL, eig.R = NULL) {
  n <- length(y)
  t <- nrow(K)
  q <- ncol(X)
  
  #  stopifnot(nrow(K) == t)
  stopifnot(ncol(K) == t)
  stopifnot(nrow(X) == n)
  
  if ( det(crossprod(X,X)) == 0 ) {
    warning("X is singular")
    return (list(REML=0,delta=0,ve=0,vg=0))
  }
  
  if ( is.null(Z) ) {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.wo.Z(K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas <- matrix(eig.R$values,n-q,m) + matrix(delta,n-q,m,byrow=TRUE)
    Etasq <- matrix(etas*etas,n-q,m)
    LL <- 0.5*((n-q)*(log((n-q)/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Lambdas)))
    dLL <- 0.5*delta*((n-q)*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Lambdas))
    
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim,eig.R$values,etas))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim,eig.R$values,etas))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.REML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas)
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root,eig.R$values, etas))
      }
    }
  }
  else {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.w.Z(Z,K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    etas.1 <- etas[1:(t-q)]
    etas.2 <- etas[(t-q+1):(n-q)]
    etas.2.sq <- sum(etas.2*etas.2)
    
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas <- matrix(eig.R$values,t-q,m) + matrix(delta,t-q,m,byrow=TRUE)
    Etasq <- matrix(etas.1*etas.1,t-q,m)
    dLL <- 0.5*delta*((n-q)*(colSums(Etasq/(Lambdas*Lambdas))+etas.2.sq/(delta*delta))/(colSums(Etasq/Lambdas)+etas.2.sq/delta)-(colSums(1/Lambdas)+(n-t)/delta))
    
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z(llim,eig.R$values,etas.1,n,t,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z(ulim,eig.R$values,etas.1,n,t,etas.2.sq))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.REML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, n=n, t1=t, etas.2.sq = etas.2.sq )
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.REML.LL.w.Z(r$root,eig.R$values, etas.1, n, t, etas.2.sq ))
      }
    }
  }  
  
  maxdelta <- exp(optlogdelta[which.max(optLL)])
  maxLL <- max(optLL)
  if ( is.null(Z) ) {
    maxva <- sum(etas*etas/(eig.R$values+maxdelta))/(n-q)    
  }
  else {
    maxva <- (sum(etas.1*etas.1/(eig.R$values+maxdelta))+etas.2.sq/maxdelta)/(n-q)
  }
  maxve <- maxva*maxdelta
  
  return (list(REML=maxLL,delta=maxdelta,ve=maxve,vg=maxva))
}

emma.eigen.R.wo.Z <- function(K, X) {
  n <- nrow(X)
  q <- ncol(X)
  S <- diag(n)-X%*%solve(crossprod(X,X))%*%t(X)
  eig <- eigen(S%*%(K+diag(1,n))%*%S,symmetric=TRUE)
  stopifnot(!is.complex(eig$values))
  return(list(values=eig$values[1:(n-q)]-1,vectors=eig$vectors[,1:(n-q)]))
}

emma.eigen.R.w.Z <- function(Z, K, X, complete = TRUE) {
  if ( complete == FALSE ) {
    vids <-  colSums(Z) > 0
    Z <- Z[,vids]
    K <- K[vids,vids]
  }
  n <- nrow(Z)
  t <- ncol(Z)
  q <- ncol(X)
  
  SZ <- Z - X%*%solve(crossprod(X,X))%*%crossprod(X,Z)
  eig <- eigen(K%*%crossprod(Z,SZ),symmetric=FALSE,EISPACK=TRUE)
  if ( is.complex(eig$values) ) {
    eig$values <- Re(eig$values)
    eig$vectors <- Re(eig$vectors)    
  }
  qr.X <- qr.Q(qr(X))
  return(list(values=eig$values[1:(t-q)],
              vectors=qr.Q(qr(cbind(SZ%*%eig$vectors[,1:(t-q)],qr.X)),
                           complete=TRUE)[,c(1:(t-q),(t+1):n)]))   
}

emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <- exp(logdelta)
  etasq <- etas*etas
  ldelta <- lambda+delta
  return( 0.5*(nq*sum(etasq/(ldelta*ldelta))/sum(etasq/ldelta)-sum(1/ldelta)) )
}

emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <-  exp(logdelta)
  return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas*etas/(lambda+delta))))-sum(log(lambda+delta))) )
}
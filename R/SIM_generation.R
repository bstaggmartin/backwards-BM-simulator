#simulate trait and rate data under an EvoRates model
#get rid of intravar, set default intracov to 0

#' Simulate rate/trait data under an evorates model
#' 
#' 
#' This function simulates branchwise rates and corresponding trait data on a given phylogeny under an
#' "evorates"-style model potentially including rate variance ("\code{R_sig2}"), a trend ("\code{R_mu}"),
#' and/or tip error ("\code{Y_sig2}"). Capable of simulating multivariate trait evolution, but this
#' feature is not well-tested.
#' 
#' 
#' @param tree An object of class "\code{phylo}".
#' @param R0 The natural log of the starting rate at the root of the phylogeny.
#' @param Rsig2 The rate variance parameter, which controls the magnitude of random
#' changes in rates over time.
#' @param X0 The starting trait value(s) at the root of the phylogeny. Recycled if
#' \code{Xsig2} and/or \code{Ysig2} imply more traits.
#' @param Rmu The trend parameter, which controls whether rates tend to decrease
#' or increase over time.
#' @param Xsig2 The "base" rate of trait evolution, which scales both \code{R0}
#' and simulated branchwise rates. You should generally leave
#' this alone unless you want to simulate multivariate trait evolution, in which
#' case this should be a vector or "rate matrix" describing evolutionary rates
#' and covariances via its diagonal and off-diagonal elements, respectively. Vectors
#' are coerced diagonal matrices with 0 off-diagonal entries (i.e., no evolutionary
#' correlations). After simulation, \code{Xsig2} is normalized to have mean rate
#' of 1, and \code{R0}/branchwise rates are rescaled accordingly. Diagonal entries
#' are recycled if \code{X0} and/or \code{Ysig2} imply more traits. Generally speaking,
#' off-diagonal entries default to 0 if unprovided.
#' @param trait.names Name(s) of simulated trait(s). If a name is not provided
#' and/or \code{NA}, defaults to "\code{X<i>}" for the \code{i}th trait.
#' @param n.obs A vector of the number of observations per tip. Assumed to be
#' in the same order as \code{tree$tip.label} if unlabeled. Missing and/or
#' \code{NA} entries default to 1 if labeled.
#' 
#' 
#' @export
sim.evorates<-function(tree,R0=0,Rsig2=1,X0=0,Rmu=0,Xsig2=1,trait.names=NULL,
                       n.obs=rep(1,length(tree$tip.label)),Ysig2=0,percent=FALSE,
                       anc.states=FALSE,slow=FALSE,res=500){
  tree<-.coerce.tree(tree)$tree
  k<-max(sapply(list(X0,Xsig2,Ysig2),NROW))
  if(length(X0)<k){
    warning('X0 implies lower number of traits than Xsig2 and/or Ysig2: recycling X0 to match other inputs')
    X0<-rep(X0,length.out=k)
  }
  if(NROW(Xsig2)<k){
    warning('Xsig2 implies lower number of traits than X0 and/or Ysig2: recycled Xsig2 (assuming no covariance if undeclared) to match other inputs')
  }else if(length(dim(Xsig2))==0&k>1){
    warning('Xsig2 is a vector: turned Xsig2 into a matrix with no covariance')
  }
  Xsig2<-.coerce.to.cov.mat(Xsig2,k)
  chol.Xsig2<-t(chol(Xsig2))
  trait.names<-rep(c(trait.names,NA),length.out=k)
  def.trait.names<-paste('X',1:k,sep='')
  trait.names<-ifelse(is.na(trait.names),def.trait.names,trait.names)
  
  n<-length(tree$tip.label)
  edge<-tree$edge
  n_e<-nrow(edge)
  ndepths<-node.depth.edgelength(tree)
  elen<-tree$edge.length
  X<-matrix(NA,n+tree$Nnode,k)
  rownames(X)<-c(tree$tip.label,n+1:tree$Nnode)
  colnames(X)<-trait.names
  X[n+1,]<-X0
  n_R<-vector('numeric',max(edge))
  n_R[n+1]<-R0
  R<-vector('numeric',n_e)
  seed<-matrix(rnorm(n_e*(k+2)),nrow=n_e)
  if(slow){
    dt<-max(ndepths)/res
    nts<-floor(elen/dt)
    dts<-elen/(nts+1)
    additional.seed<-lapply(nts,function(ii) rnorm(ii+1))
  }
  for(i in 1:n_e){
    t<-elen[i]
    if(slow){
      tmp<-n_R[edge[i,1]]+cumsum(additional.seed[[i]]*sqrt(Rsig2*dts[i])+Rmu*dts[i])
      R[i]<-log(mean(exp(c(n_R[edge[i,1]],tmp))))
      n_R[edge[i,2]]<-tmp[length(tmp)]
    }else{
      n_R[edge[i,2]]<-n_R[edge[i,1]]+seed[i,1]*sqrt(Rsig2*t)
      R[i]<-mean(n_R[edge[i,]])+seed[i,2]*sqrt(Rsig2*t/12)
      if(Rmu!=0&t!=0){
        R[i]<-R[i]-log(abs(Rmu))-log(t)+log(abs(diff(exp(Rmu*ndepths[edge[i,]]))))
      }
    }
    X[edge[i,2],]<-X[edge[i,1],]+sqrt(t*exp(R[i]))*chol.Xsig2%*%seed[i,-(1:2)]
  }
  scalar<-1/mean(diag(Xsig2))
  Xsig2<-Xsig2*scalar
  R<-R-log(scalar);R0<-R0-log(scalar)
  out<-list('R_0'=R0,'R_sig2'=Rsig2,'X_0'=X0,'R'=R,'R_mu'=Rmu,'tree'=tree)
  if(anc.states){
    out$X<-X
  }else{
    out$X<-as.matrix(X[1:n,])
  }
  if(k>1){
    out$Xsig2<-Xsig2
    out$k<-k
  }
  out$trait.data<-as.matrix(X[1:n,])
  if(all(Ysig2==0)&any(n.obs>1)){
    warning('Ysig2 implies no intra-tip variation: set n.obs to 1 for each tip')
  }
  if(any(Ysig2>0)){
    if(NROW(Ysig2)<k){
      warning('Ysig2 implies lower number of traits than X0 and/or Xsig2: recycled Ysig2 (assuming no covariance if undeclared) to match other inputs')
    }else if(length(dim(Ysig2))==0&k>1){
      warning('Ysig2 is a vector: turned Ysig2 into a matrix with no covariance')
    }
    Ysig2<-.coerce.to.cov.mat(Ysig2,k)
    if(percent){
      vars<-apply(X,2,var)
      diag(Ysig2)<-diag(Ysig2^2)*vars
      Ycor<-Ysig2
      diag(Ycor)<-0
      if(any(Ycor<=-1)|any(Ycor>=1)){
        stop('If Ysig2 is expressed in terms of percentages, off-diagonal elements should represent correlations and be between -1 and 1')
      }
      diag(Ycor)<-1
      Ysig2<-diag(sqrt(diag(Ysig2)),nrow=k)%*%Ycor%*%diag(sqrt(diag(Ysig2)),nrow=k)
    }
    chol.Ysig2<-t(chol(Ysig2))
    if(is.null(names(n.obs))){
      n.obs<-rep(n.obs,length.out=n)
    }else{
      n.obs<-n.obs[tree$tip.label]
      n.obs[is.na(n.obs)]<-1
      names(n.obs)<-tree$tip.label
    }
    seed<-matrix(rnorm(k*sum(n.obs)),k,sum(n.obs))
    Y<-X[rep(1:n,n.obs),]+t(chol.Ysig2%*%seed)
    rownames(Y)<-rep(tree$tip.label,n.obs)
    out$trait.data<-Y
    out$Ysig2<-Ysig2
    out$n.obs<-n.obs
  }
  class(out)<-'evorates'
  out
}
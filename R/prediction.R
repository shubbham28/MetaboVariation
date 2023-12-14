#' Title
#' @title
#' Predict method for mixed models fitted with MCMCglmm in MetaboVariation.
#' @description
#' Compute posterior draws of the posterior predictive distribution. Can be performed for the data used to fit the model (posterior predictive checks) or for new data.
#' It is similar to the function \code{\link[MCMCglmm:predict.MCMCglmm]{predict}} from \code{\link[MCMCglmm]{MCMCglmm}} package.
#' @param object An object of class \code{\link{MetaboVariation}}.
#' @param seed The seed for random number generation to make results reproducible. A seed is default to 19205033.
#' @param nsim number of response vectors to simulate. Defaults to number of effective sample size of the model.
#' @param newdata An optional data.frame for which to evaluate predictions. If NULL (default), the original data of the model is used.
#' @param level A numeric scalar in the interval (0,1) giving the target probability content of the intervals.
#' @param verbose logical; if TRUE, warnings are issued with newdata when the original model has fixed effects that do not appear in newdata and/or newdata has random effects not present in the original model.
#'
#' @return A matrix of expectated values and Highest Posterior Distribution interval.
#' @export
#'
#'
prediction<-function(object,seed = NULL,nsim=nrow(object$Sol),newdata=NULL,level = 0.95,verbose = FALSE){
  if(!is.null(seed)){
    if(!is.integer(seed)){
      stop("seed must be an integer")
    }else{
      seed = 19205033
      set.seed(seed)
    }
  }
  if(!is.null(newdata)){
    suppressWarnings(object2<-MCMCglmm::MCMCglmm(fixed=object$Fixed$formula, random=object$Random$formula, rcov=object$Residual$formula, family=object$Residual$original.family, data=newdata, nitt=1, thin=1, burnin=0, ginverse=object$ginverse, verbose=FALSE, pr=TRUE, start=list(QUASI=FALSE)))
    find.fixed<-match(colnames(object2$Sol)[1:object2$Fixed$nfl], colnames(object$Sol))
    find.random<-match(colnames(object2$Sol)[-c(1:object2$Fixed$nfl)], colnames(object$Sol))
    if(any(is.na(find.fixed))){stop("model for newdata has fixed effects not present in original model")}
    if(verbose){
      if(any(is.na(find.random))){
        missing.random<-colnames(object2$Sol)[which(is.na(find.random))+object2$Fixed$nfl]
        warning(paste("model for newdata has random effects not present in original model:", paste(missing.random, collapse=", ")))
      }
      missing.fixed<-which(!colnames(object$Sol)[1:object$Fixed$nfl]%in%colnames(object2$Sol)[1:object2$Fixed$nfl])
      if(length(missing.fixed)>0){
        missing.fixed<-colnames(object$Sol)[1:object$Fixed$nfl][missing.fixed]
        warning(paste("original model has fixed effects not present in newdata:", paste(missing.fixed, collapse=", ")))
      }
    }
    object2$Sol<-object$Sol[,c(find.fixed, find.random), drop=FALSE]
    find.vcv<-match(colnames(object2$VCV), colnames(object$VCV))
    if(any(is.na(find.vcv))){stop("model for newdata has (co)variance terms not in original model")}
    object2$VCV<-object$VCV[,find.vcv, drop=FALSE]
    if(!is.null(object2$CP)){
      find.cp<-match(colnames(object2$CP), colnames(object$CP))
      if(any(is.na(find.cp))){stop("model for newdata has cutpoints not in original model")}
      object2$CP<-object$CP
    }
    object<-object2
    rm(object2)
  }
  ynew<-matrix(NA, nrow(object$X), nsim)
  unew<-matrix(NA,sum(object$Random$nfl*object$Random$nrl),1)
  enew<-matrix(NA,sum(object$Residual$nfl*object$Residual$nrl),1)

  it<-sample(1:nrow(object$Sol), nsim)

  for (i in 1:nsim) {
    cnt<-0
    cnt2<-0

    for(j in 1:length(object$Random$nfl)){

      nfl<-object$Random$nfl[j]
      nrl<-object$Random$nrl[j]
      nat<-object$Random$nat[j]
      set.seed((seed + i*nsim + j + i*j))
      Y<-matrix(rnorm(nrl*nfl),nrl,nfl)

      unew[1:(nrl*nfl)+cnt]<-as.vector(Y%*%chol(sqrt(matrix(object$VCV[it[i],cnt2+1:(nfl^2)],nfl,nfl))))
      cnt<-cnt+(nrl*nfl)
      cnt2<-cnt2+nfl^2
    }
    cnt<-0
    for(j in 1:length(object$Residual$nfl)){

      nfl<-object$Residual$nfl[j]
      nrl<-object$Residual$nrl[j]

      Y<-matrix(rnorm(nrl*nfl),nrl,nfl)

      enew[1:(nrl*nfl)+cnt]<-as.vector(Y%*%chol(matrix(object$VCV[it[i],cnt2+1:(nfl^2)],nfl,nfl)))

      cnt<-cnt+(nrl*nfl)
      cnt2<-cnt2+nfl^2
    }
    ynew[,i]<-as.vector(object$X%*%object$Sol[it[i],1:ncol(object$X)]+object$ZR%*%enew+object$Z%*%object$Sol[it[i],(ncol(object$X) + 1):ncol(object$Sol)])

    if(!is.null(object$Random$nfl)){

      #    ynew[,i]<-ynew[,i]+as.vector(object$Z%*%unew)


    }
  }
  if(nsim==1){
    ynew<-as.vector(ynew)
  }
  else{
    ynew = t(ynew)
  }
  pred<-matrix(colMeans(ynew), dim(ynew)[2],1)
  pred<-cbind(pred, coda::HPDinterval(coda::mcmc(ynew), prob=level))
  colnames(pred)<-c("fit", "lwr", "upr")
  rownames(pred)<-1:dim(pred)[1]
  result = list("posterior" = ynew,"summary"=pred)
  return(result)
}

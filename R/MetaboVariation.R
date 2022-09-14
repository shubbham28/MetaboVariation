`%dopar%` <- foreach::`%dopar%`
#' @title
#' Fits a Bayesian generalised linear model to repeated measurements of metabolites.
#' @description
#' Implements a Bayesian generalised linear model (BGLM) to flag individuals with notable variation in their metabolite levels.
#' @details
#' The Bayesian generalised linear model (BGLM) to identify individuals with variations in their
#' metabolite levels is fitted. The BGLM can include covariates and learns the intra-individual variation in metabolite levels using repeated measures.
#' Posterior predictive distributions of metabolite levels at the individual level are provided. These are used to identify individuals with observed metabolite level outside the
#' \code{cutoff}% Highest Posterior Distribution (HPD) interval at one time point; such individuals with notable variation are flagged.
#' The model builds on the Bayesian framework using the function \code{\link[MCMCglmm:MCMCglmm]{MCMCglmm}} from the \code{\link[MCMCglmm]{MCMCglmm}} package.
#'
#' The BGLM models the metabolite value \eqn{y^m_{it}} for the individual \eqn{i} at timepoint \eqn{j},
#' \deqn{ y^m_{it} = \beta_{0} + \sum_{l=1}^L\beta_{1,l}x_{il} + S_{it} + \epsilon_{it}}
#'
#' where \eqn{y^m_{it}} is the metabolite level for individual \eqn{i} at time point \eqn{t}, the term \eqn{X_{il}} is the covariate \eqn{l \in 1, 2, ..., L} for individual \eqn{i}. The term \eqn{\beta_0} is the mean intercept, while \eqn{\beta_{1,l}} is the regression coefficient for the covariate \eqn{l}. The random effect for the \eqn{i^{th}} individual is denoted by \eqn{S_{it}} for timepoint \eqn{t} and \eqn{\epsilon_{it}} is the random error for individual \eqn{i} at timpoint \eqn{t} assumed \eqn{\epsilon_{ij} \sim N(0,\sigma^2)}.
#'
#' @param data A data frame containing data of all variables to be used in the BGLM.
#' @param individual_ids A character string detailing the name of the column in the data that contains the individual ids.
#' @param metabolite A string or vector of string containing the name of metabolites to model.
#' @param covariates A list of string containing the name of columns in the data that contains the covariates.
#' @param iter Number of total iterations per chain (including warmup; defaults to 2000).
#' @param warmup A positive integer specifying number of warmup iterations. This also specifies the number of iterations used for stepsize
#'   adaptation, so warmup draws should not be used for inference. The number
#'   of warmup should not be larger than \code{iter} and the default is
#'   \code{iter/2}.
#' @param thin Thinning rate. Must be a positive integer. Set thin > 1 to save memory and computation time if \code{iter} is large.
#' @param cutoff Cutoff for the width of Highest Posterior Distribution interval of the posterior predictive distribution. The value should be less than 100. Default is set to 95.
#' @param cores Number of cores to use when executing the chains in parallel, which defaults to half of the total cores you have in the system.
#' @param seed The seed for random number generation to make results reproducible. A seed is default to 19205033.
#' @return The function returns an object of \code{\link{MetaboVariation}}. The Following values are returned for each metabolite object.
#' \itemize{
#' \item metabolite - contains the name of the metabolite.
#' \item significant_covariates - a list of covariates that have significant relationship (at the 90% level) with the metabolite levels in individuals.
#' \item result - a summary of the posterior predictive distributions of the individuals in the cohort. The result contains the mean, \code{cutoff} % HPD interval width, tails of the HPD interval and the observed value of the metabolite for the individual for that timepoint. The result also has a binary flag that shows whether the observed value lies within the HPD interval or not where 1 denotes the observed value lies in the interval and 0 denotes the observed value lies outside the interval.
#' \item warmup - number of total warmup iterations per chain.
#' \item iter - number of total iterations per chain used for modelling.
#' \item Rhat - the potential scale reduction statistic, also known as the Gelman-Rubin statistic which measures the extent to which chains are converging. The further the value of the statistic from 1, the poorer the convergence of the chains. The MCMC chains are considered converged if the value lies between 0.9 and 1.05.
#' }
#' @export
#' @seealso \code{\link{flagged.Individuals}}, \code{\link{plot.MetaboVariation}}
#' @examples
#' \dontrun{
#' data(metabol.data)
#' metabolite_list = colnames(metabol.data)[5:length(colnames(metabol.data))]
#' metabolites = get.metabolites(list = metabolite_list)
#' covariates = c("SexM.1F.2","Age","BMI")
#' individual_id = "Individual_id"
#' model = MetaboVariation(data = metabol.data,individual_ids = individual_id,
#' metabolite = metabolites[1], covariates = covariates)
#' model = MetaboVariation(data = metabol.data,individual_ids = individual_id,
#' metabolite = metabolites[1:3])
#' }
#'
#' @import parallel foreach stats

MetaboVariation <- function ( data,individual_ids,metabolite,covariates=NULL,iter = 5000,
          warmup = floor(iter/2),thin = 2,cutoff=0.95,cores = NULL,seed = NULL){
  if(!is.null(seed)){
    if(!is.integer(seed)){
      stop("seed must be an integer")
    }else{
      seed = 19205033
      set.seed(seed)
    }
  }
  future = getOption("future", TRUE)
  file = paste(metabolite,"model",sep="_")
  i <- NULL
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    cores <- getOption("mc.cores",2)
  } else {
    # use half of the cores in devtools::test()
    if(is.null(cores)){
      cores <- getOption("mc.cores",parallel::detectCores()/2)
    }
  }

  col_list = colnames(data)
  if(!(all(covariates %in% col_list)) & !(is.null(covariates))){
    stop("Covariates are not present in dataset")
  }
  if(length(individual_ids)!=1){
    stop("Multiple columns for individual id are passed")
  }
  if(!(individual_ids %in% col_list)){
    stop("Individual id column is not present in dataset")
  }
  if(length(covariates)==0){
    formula = paste("measurement ~ 1")
  }
  else{
    formula = paste("measurement ~",paste(covariates,collapse = " + "))
  }
  if(length(metabolite)==1){
    var_list = na.omit(stringr::str_extract(col_list,stringr::regex(paste(".*",metabolite,".*",sep = ""))))
    if(length(var_list)==0){
      stop("Metabolite passed is not present in dataset")
    }
    new_data = reshape2::melt(data,id.vars = c(individual_ids,covariates),measure.vars = var_list,variable.name = "occurrence",value.name = "measurement")
    new_data = na.omit(new_data)
    print(paste("Model building for",metabolite))
    set.seed(seed)
    model <- parallel::mclapply(1:cores, function(i) {
      MCMCglmm::MCMCglmm(as.formula(formula),random = as.formula(paste("~us(1):",individual_ids)),data = new_data,pr=TRUE,nitt = iter,thin = thin,burnin = warmup)
    }, mc.cores=cores)
    model <- lapply(model, function(m) m$Sol)
    model <- do.call(coda::mcmc.list, model)
    t = summary(model)
    rf = coda::gelman.diag(model)
    print(paste("Calculating posterior predictive distribution for",metabolite))
    div = nrow(new_data)%/% 25
    obs = nrow(new_data)%/% div - 1
    doParallel::registerDoParallel(cores)
    result<-foreach::foreach(i=0:obs,.combine=rbind,.packages = c('MCMCglmm'))%dopar%{
      # for (i in 0:obs) {
      if(i==obs){
        drop = c((i*div + 1):nrow(new_data))
      }
      else{
        drop = i*div + c(1:div)
      }
      set.seed(seed)
      model = MCMCglmm::MCMCglmm(as.formula(formula),random = as.formula(paste("~us(1):",individual_ids)),data = new_data[-drop,],pr=TRUE,
                                 nitt = iter,thin = thin,burnin = warmup)
      prediction(model,newdata = new_data[drop,])$summary
    }

    brief_result = cbind(result,result[,3] - result[,2],new_data$measurement)
    brief_result = cbind(brief_result,(brief_result[,2] < brief_result[,5] & brief_result[,5] < brief_result[,3]))
    colnames(brief_result) = c("mean",paste0((.5-(cutoff)/2)*100," cutoff"),paste0((.5+(cutoff)/2)*100," cutoff"),"interval","original","identifier")
    rownames(brief_result) = paste(new_data[,individual_ids],new_data$occurrence)
    co_list = cbind(t$statistics[,1:2],t$quantiles[,c(1,5)],apply(t$quantiles, 1, function(x){ifelse(x[1]*x[5] > 0,TRUE,FALSE) }))
    colnames(co_list) = c("Mean","SD","2.5%","97.5%","Significance")
    co_list = co_list[-1,]
    return_list = list("significant_covariates" = co_list,"result" = brief_result,"warmup" = warmup,"iter"=iter,
                       "Rhat" =  mean(c(rf$psrf[,1],rf$mpsrf)),"metabolite" = metabolite)

    remove(model)
    gc()
    class(return_list) = c("MetaboVariation","meta.single_model")
    return(return_list)
  }
    else{
      final_list = list()
      for (met in metabolite) {
        var_list = na.omit(stringr::str_extract(col_list,stringr::regex(paste(".*",met,".*",sep = ""))))
        if(length(var_list)==0){
          stop("Metabolite passed is not present in dataset")
        }
        new_data = reshape2::melt(data,id.vars = c(individual_ids,covariates),measure.vars = var_list,variable.name = "occurrence",value.name = "measurement")
        new_data = na.omit(new_data)
        print(paste("Model building for",met))
        set.seed(seed)
        model <- parallel::mclapply(1:cores, function(i) {
          MCMCglmm::MCMCglmm(as.formula(formula),random = as.formula(paste("~us(1):",individual_ids)),data = new_data,pr=TRUE,nitt = iter,thin = thin,burnin = warmup)
        }, mc.cores=cores)
        model <- lapply(model, function(m) m$Sol)
        model <- do.call(coda::mcmc.list, model)
        t = summary(model)
        rf = coda::gelman.diag(model)
        print(paste("Calculating posterior predictive distribution for",met))
        div = nrow(new_data)%/% 25
        obs = nrow(new_data)%/% div - 1
        doParallel::registerDoParallel(cores)
        result<-foreach::foreach(i=0:obs,.combine=rbind,.packages = c('MCMCglmm'))%dopar%{
          # for (i in 0:obs) {
          if(i==obs){
            drop = c((i*div + 1):nrow(new_data))
          }
          else{
            drop = i*div + c(1:div)
          }
          set.seed(seed)
          model = MCMCglmm::MCMCglmm(as.formula(formula),random = as.formula(paste("~us(1):",individual_ids)),data = new_data[-drop,],pr=TRUE,
                                     nitt = iter,thin = thin,burnin = warmup)
          prediction(model,newdata = new_data[drop,])$summary
        }

        brief_result = cbind(result,result[,3] - result[,2],new_data$measurement)
        brief_result = cbind(brief_result,(brief_result[,2] < brief_result[,5] & brief_result[,5] < brief_result[,3]))
        colnames(brief_result) = c("mean",paste0((.5-(cutoff)/2)*100," cutoff"),paste0((.5+(cutoff)/2)*100," cutoff"),"interval","original","identifier")
        rownames(brief_result) = paste(new_data[,individual_ids],new_data$occurrence)
        co_list = cbind(t$statistics[,1:2],t$quantiles[,c(1,5)],apply(t$quantiles, 1, function(x){ifelse(x[1]*x[5] > 0,TRUE,FALSE) }))
        colnames(co_list) = c("Mean","SD","2.5%","97.5%","Significance")
        co_list = co_list[-1,]
        final_list[[met]] = list("metabolite" = met,"significant_covariates" = co_list,"result" = brief_result,"warmup" = warmup,"iter"=iter,
                                 "Rhat" = mean(c(rf$psrf[,1],rf$mpsrf)))

        class(final_list[[met]]) = c("MetaboVariation","meta.single_model")
        remove(model)
        gc()

      }
      class(final_list) = c("MetaboVariation","meta.multi_model")
      return(final_list)

    }
}




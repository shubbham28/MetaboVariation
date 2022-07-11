`%dopar%` <- foreach::`%dopar%`
#' @title
#' Fits a Bayesian generalised linear model to repeated measurements of metabolites.
#' @description
#' Implements a Bayesian generalised linear model (BGLM) to identify the flagged individuals with notable variations.
#' @details
#' The Bayesian generalised linear model (BGLM) to identify individuals with variations in their
#' metabolite levels is fitted by employing STAN. The BGLM can include covariates and learns the intra-individual variation in metabolite levels using repeated measures.
#' Posterior predictive distributions of metabolite levels at the individual level are provided. These are used to identify individuals with observed metabolite level outside the
#' \code{cutoff}% central credible interval at one time point; such individuals with notable variation are flagged.
#' The model builds on the Bayesian framework using the function \code{\link[brms:brm]{brm}} from \code{\link[brms:brm]{brms}} package.
#'
#' The BGLM models the metabolite value \eqn{M_{ij}} for the individual \eqn{i} at timepoint \eqn{j},
#' \deqn{ y^m_{it} = \beta_{0} + \sum_{l=1}^L\beta_{1,l}x_{il} + S_{it} + \epsilon_{it}}
#'
#' where \eqn{y^m_{it}} is the metabolite level for individual \eqn{i} at time point \eqn{t}, the term \eqn{X_{il}} is the covariate \eqn{l \in 1, 2, ..., L} for individual \eqn{i}. The notation \eqn{\beta_0} is the mean intercept, while \eqn{\beta_{1,l}} is the regression coefficient for the covariate \eqn{l}. The random effect for the \eqn{i^{th}} individual is denoted by \eqn{S_{it}} for timepoint \eqn{t} and \eqn{\epsilon_{it}} is the random error for individual \eqn{i} at timpoint \eqn{t} assumed \eqn{\epsilon_{ij} \sim N(0,\sigma^2)}.
#'
#' @param data A data frame containing data of all variables to be used in the BGLM.
#' @param individual_ids A character string detailing the name of the column in the data that contains the individual ids.
#' @param metabolite A string or vector of string containing the name of metabolites to model.
#' @param covariates A list of string containing the name of columns in the data that contains the covariates.
#' @param save_brms_model A logical value indicating whether the model should be saved or not. By default, it is set as FALSE.
#' @param full_posterior A logical value indicating whether the full posterior distribution should be returned or not. By default, it is set as FALSE.
#' @param iter Number of total iterations per chain (including warmup; defaults to 2000).
#' @param warmup A positive integer specifying number of warmup (aka burnin)
#'   iterations. This also specifies the number of iterations used for stepsize
#'   adaptation, so warmup draws should not be used for inference. The number
#'   of warmup should not be larger than \code{iter} and the default is
#'   \code{iter/2}.
#' @param thin Thinning rate. Must be a positive integer. Set thin > 1 to save memory and computation time if iter is large.
#' @param cutoff Cutoff for the credible interval of posterior predictive distribution. The value should be less than 100. Default is set to 95.
#' @param ... Further arguments are passed to \code{\link[brms]{brm}} for backend processes.
#'
#' @return function returns an object of "MetaboVariation". Following values are returned for each metabolite object.
#' \itemize{
#' \item metabolite - contains the name of metabolite.
#' \item significant_covariates - a list of covariates that have significant relationship (at the 90% level) with the metabolite levels in individuals.
#' \item result - a summary of the posterior predictive distribution of the cohort. The result contains the mean, \code{cutoff} % credible interval width, tails of the credible interval and the observed value of the metabolite for the individual for that timepoint. The result also has a binary flag that shows whether the observed value lies within the credible interval or not.
#' \item warmup - number of total warmup iterations per chain.
#' \item iter - number of total iterations per chain used for modelling.
#' \item Rhat - the potential scale reduction statistic, also known as the Gelman-Rubin statistic which measures the extent to which chains are converging. The further the value of the statistic from 1, the poorer the convergence of the chains. The MCMC chains are considered converged if the value lies between 0.9 and 1.05.
#' \item full_posterior - it contains the full posterior predictive distribution of the data.
#' }
#' @export
#' @seealso \code{\link{flaggedIndividuals}}, \code{\link{plot.MetaboVariation}}
#' @examples
#' \dontrun{
#' data(metabol.data)
#' metabolite_list = colnames(metabol.data)[5:length(colnames(metabol.data))]
#' metabolites = getmetabolites(list = metabolite_list)
#' covariates = c("SexM.1F.2","Age","BMI")
#' individual_id = "Individual_id"
#' model = MetaboVariation(data = metabol.data,individual_ids = individual_id,
#' metabolite = metabolites[1], covariates = covariates)
#' model = MetaboVariation(data = metabol.data,individual_ids = individual_id,
#' metabolite = metabolites[1:3], full_posterior = TRUE)
#' }
#'
#' @import parallel foreach stats

MetaboVariation <- function ( data,individual_ids,metabolite,covariates=NULL, save_brms_model = FALSE,full_posterior = FALSE,iter = 2000,
          warmup = floor(iter/2),thin = 1,cutoff=95,...){

  future = getOption("future", TRUE)
  file = paste(metabolite,"model",sep="_")
  seed = 19205033
  i <- NULL
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    cores <- getOption("mc.cores",2)
  } else {
    # use half of the cores in devtools::test()
    cores <- getOption("mc.cores",parallel::detectCores()/2)
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
    formula = paste("measurement ~ (1|",individual_ids,")")
  }
  else{
    formula = paste("measurement ~",paste(covariates,collapse = " + "), "+ (1|",individual_ids,")")
  }
  lower_cutoff = round((100 - cutoff)/2,1)
  upper_cutoff = round(100 - lower_cutoff,1)
  if(length(metabolite)==1){
    var_list = na.omit(stringr::str_extract(col_list,stringr::regex(paste(".*",metabolite,".*",sep = ""))))
    if(length(var_list)==0){
      stop("Metabolite passed is not present in dataset")
    }
    new_data = reshape2::melt(data,id.vars = c(individual_ids,covariates),measure.vars = var_list,variable.name = "occurrence",value.name = "measurement")
    new_data = na.omit(new_data)
    print(paste("Model building for",metabolite))
    model = brms::brm(formula,data = new_data,iter = iter,
                      warmup = warmup, cores = cores, seed = seed, file = file,thin = thin,backend = "cmdstanr")
    rows <- nrow(new_data)
    doParallel::registerDoParallel(cores)
    print(paste("Calculating posterior predictive distribution for",metabolite))
    result<-foreach::foreach(i=1:rows,.combine=cbind,.packages = c('brms','cmdstanr','future'))%dopar%{
      new_model = update(model,file=file,newdata = new_data[-i,],silent=2,refresh=0)
      rbind.data.frame(paste(new_data[i,individual_ids], new_data$occurrence[i]),brms::posterior_predict(new_model,newdata = new_data[i,]))
    }
    if(!save_brms_model){
      file.remove(paste(file,".rds",sep=""))
    }
    colnames(result) = result[1,]
    result = result[-1,]
    result = apply(result,2,as.numeric)
    brief_result = t(apply(result, 2, function(x) {c(mean(x),quantile(x,probs = c(lower_cutoff/100,upper_cutoff/100)),quantile(x,upper_cutoff/100)-quantile(x,lower_cutoff/100))}))
    brief_result = cbind(brief_result,new_data$measurement)
    brief_result = cbind(brief_result,(brief_result[,2] < brief_result[,5] & brief_result[,5] < brief_result[,3]))
    colnames(brief_result) = c("mean",paste0(lower_cutoff,"%"),paste0(upper_cutoff,"%"),"QR","original","identifier")
    co_list = summary(model)$fixed[names(na.omit(apply(summary(model)$fixed[2:nrow(summary(model)$fixed),3:4],1,function(x) {ifelse(prod(x[1],x[2])>=0,1,NA )}))),3:4]
    if(nrow(co_list)==0){
      co_list = NULL
    }

    return_list = list("significant_covariates" = co_list,"result" = brief_result,"warmup" = warmup,"iter"=iter,
                       "Rhat" = summary(model)$spec_pars$Rhat)
    if(full_posterior){
      return_list[["full_posterior"]] = result
    }
    remove(model)
    gc()
    class(return_list) = c("MetaboVariation","meta.single_model")
    return(return_list)
  }
    else{
      final_list = list()
      for (met in metabolite) {
        file = paste(met,"model",sep="_")
        var_list = na.omit(stringr::str_extract(col_list,stringr::regex(paste(".*",met,".*",sep = ""))))
        if(length(var_list)==0){
          stop("Metabolite passed is not present in dataset")
        }
        new_data = reshape2::melt(data,id.vars = c(individual_ids,covariates),measure.vars = var_list,variable.name = "occurrence",value.name = "measurement")
        new_data = na.omit(new_data)
        print(paste("Model building for",met))
        model = brms::brm(formula,data = new_data,iter = iter,
                          warmup = warmup, cores = cores, seed = seed, file = file,thin = thin,backend = "cmdstanr")

        rows <- nrow(new_data)
        doParallel::registerDoParallel(cores)
        print(paste("Calculating posterior predictive distribution for",met))
        result<-foreach::foreach(i=1:rows,.combine=cbind,.packages = c('brms','cmdstanr','future'))%dopar%{
          new_model = update(model,file=file,newdata = new_data[-i,],silent=2,refresh=0)
          rbind.data.frame(paste(new_data[i,individual_ids], new_data$occurrence[i]),brms::posterior_predict(new_model,newdata = new_data[i,]))
        }
        if(!save_brms_model){
          file.remove(paste(file,".rds",sep=""))
        }
        colnames(result) = result[1,]
        result = result[-1,]
        result = apply(result,2,as.numeric)
        brief_result = t(apply(result, 2, function(x) {c(mean(x),quantile(x,probs = c(lower_cutoff/100,upper_cutoff/100)),quantile(x,lower_cutoff/100)-quantile(x,upper_cutoff/100))}))
        brief_result = cbind(brief_result,new_data$measurement)
        brief_result = cbind(brief_result,(brief_result[,2] < brief_result[,5] & brief_result[,5] < brief_result[,3]))
        colnames(brief_result) = c("mean",paste0(lower_cutoff,"%"),paste0(upper_cutoff,"%"),"QR","original","identifier")
        co_list = summary(model)$fixed[names(na.omit(apply(summary(model)$fixed[2:nrow(summary(model)$fixed),3:4],1,function(x) {ifelse(prod(x[1],x[2])>=0,1,NA )}))),3:4]
        if(nrow(co_list)==0){
          co_list = NULL
        }
        final_list[[met]] = list("metabolite" = met,"significant_covariates" = co_list,"result" = brief_result,"warmup" = warmup,"iter"=iter,
                                 "Rhat" = summary(model)$spec_pars$Rhat)
        if(full_posterior){
          final_list[[met]][["full_posterior"]] = result
        }
        class(final_list[[met]]) = c("MetaboVariation","meta.single_model")
        remove(model)
        gc()

      }
      class(final_list) = c("MetaboVariation","meta.multi_model")
      return(final_list)

    }
}


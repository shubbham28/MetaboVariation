`%dopar%` <- foreach::`%dopar%`
#' @title
#' Fits a Bayesian generalised linear model to repeated measurements of metabolites.
#' @description
#' Implements a Bayesian generalised linear model (BGLM) to flag individuals with notable variation in their metabolite levels.
#' @details
#' The Bayesian generalised linear model (BGLM) is used to identify individuals with variations in their
#' metabolite levels. The BGLM can include covariates and learns the intra-individual variation in metabolite levels from the repeated measurements.
#' Posterior predictive distributions of metabolite levels at the individual level are provided. These are used to identify individuals with observed metabolite level outside the
#' \code{cutoff}% Highest Posterior Distribution (HPD) interval at one time point; such individuals with notable variation are flagged.
#' The model builds on the Bayesian framework using the function \code{\link[MCMCglmm:MCMCglmm]{MCMCglmm}} from the \code{\link[MCMCglmm]{MCMCglmm}} package.
#'
#' The BGLM models the metabolite value \eqn{y^m_{it}} for the individual \eqn{i} at timepoint \eqn{t},
#' \deqn{ y^m_{it} = \beta_{0} + \sum_{l=1}^L\beta_{1,l}x_{il} + S_{it} + \epsilon_{it}}
#'
#' where \eqn{y^m_{it}} is the metabolite level for individual \eqn{i} at time point \eqn{t}, the term \eqn{x_{il}} is the covariate \eqn{l \in 1, 2, ..., L} for individual \eqn{i}. The term \eqn{\beta_0} is the mean intercept, while \eqn{\beta_{1,l}} is the regression coefficient for the covariate \eqn{l}. The random effect for the \eqn{i^{th}} individual is denoted by \eqn{S_{it}} for timepoint \eqn{t} and \eqn{\epsilon_{it}} is the random error for individual \eqn{i} at timpoint \eqn{t} assumed \eqn{\epsilon_{ij} \sim N(0,\sigma^2)}.
#'
#' @param data A data frame containing data of all variables to be used in the BGLM. You can refer \code{\link{metabol.data}} for structure of the data.
#' @param individual_ids A character string detailing the name of the column in the data that contains the individual ids.
#' @param metabolites A string or vector of string containing the name of metabolites to model.
#' @param covariates A list of string containing the name of columns in the data that contains the covariates.
#' @param iter Number of total iterations per chain (including warmup; defaults to 2000).
#' @param warmup A positive integer specifying number of warmup iterations. This also specifies the number of iterations used for stepsize
#'   adaptation, so warmup draws should not be used for inference. The number
#'   of warmup should not be larger than \code{iter} and the default is
#'   \code{iter/2}.
#' @param thin Thinning rate. Must be a positive integer. Set thin > 1 to save memory and computation time if \code{iter} is large.
#' @param cutoff A list of cutoff for the width of Highest Posterior Distribution interval of the posterior predictive distribution. The value should be less than 1. Default cutoff are set to a list of (0.91,0.93,0.95,0.97,0.99)
#' @param cores Number of cores to use when executing the chains in parallel, which defaults to half of the total cores you have in the system.
#' @param seed The seed for random number generation to make results reproducible. A seed is default to 19205033.
#' @param type Determines whether to model with or without considering the correlation among metabolites. Options: "dependent" or "independent".
#' @param n_chains Number of chains to use in the model.
#' @param prior list of prior specifications having 4 possible elements: R (random effects of individuals) G (random errors of the model), B (fixed effects). B is a list containing the expected value (mu) and a (co)variance matrix (V) representing the strength of belief. The priors for the variance structures (R and G) are lists with the expected (co)variances (V) and degree of belief parameter (nu) for the inverse-Wishart where nu should be greater than number of metabolites.
#' @return The function returns an object of \code{\link{MetaboVariation}}. The Following values are returned for each metabolite object.
#' \itemize{
#' \item model - contains intermediate statistical models used to see the convergence.
#' \item type - shows which modelling is done. It can be either "dependent" or "independent"
#' \item significant_covariates - a list of covariates that have significant relationship (at the 95% level) with the metabolites.
#' \item result - a summary of the posterior predictive distributions of the individuals in the cohort. The result contains the mean, \code{cutoff} % HPD interval width, tails of the HPD interval and the observed value of the metabolite for the individual for that timepoint. The result also has a binary flag that shows whether the observed value lies within the HPD interval or not where 1 denotes the observed value lies in the interval and 0 denotes the observed value lies outside the interval.
#' \item chain_convergence - the potential scale reduction statistic, also known as the Gelman-Rubin statistic which measures the extent to which chains are converging. The further the value of the statistic from 1, the poorer the convergence of the chains. The MCMC chains are considered converged if the value lies between 0.9 and 1.05.
#' }
#' @export
#' @references Hadfield, J. D. (2010). MCMC Methods for Multi-Response Generalized Linear Mixed Models: The MCMCglmm R Package. Journal of Statistical Software, 33(2), 1–22. https://doi.org/10.18637/jss.v033.i02
#' @seealso \code{\link{radar.plot}}, \code{\link{circos.plot}}, \code{\link{metabolite.heatmap}}
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
#' @import parallel foreach

MetaboVariation <- function ( data,individual_ids,metabolites,covariates=NULL,iter = 5000,
                                warmup = 2000,thin = 2,cutoff=c(0.95,0.975,0.99),cores = NULL,seed = NULL,type="dependent",n_chains=4,prior=NULL){
    if(!is.null(seed)){
      if(!is.numeric(seed)){
        stop("seed must be an integer")
      }
    }else{
      seed = 19205033
    }
    set.seed(seed)
    #file = paste(metabolites,"model",sep="_")
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

    var_list = stats::na.omit(stringr::str_extract(col_list,stringr::regex(paste(".*",metabolites[1],".*",sep = ""))))
    final_data = reshape2::melt(data,id.vars = c(individual_ids,covariates),measure.vars = var_list,variable.name = "occurrence",value.name = metabolites[1])
    s = gsub(metabolites[1],"",var_list)
    levels(final_data$occurrence) = gsub("\\D","",s)
    for (i in 2:length(metabolites)) {
      var_list = stats::na.omit(stringr::str_extract(col_list,stringr::regex(paste(".*",metabolites[i],".*",sep = ""))))
      new_data = reshape2::melt(data,id.vars = c(individual_ids,covariates),measure.vars = var_list,variable.name = "occurrence",value.name = metabolites[i])
      s = gsub(metabolites[i],"",var_list)
      levels(new_data$occurrence) = gsub("\\D","",s)
      final_data <- merge(final_data,new_data,by=c(c(individual_ids,covariates),"occurrence"))
    }
    final_data = final_data[order(final_data$occurrence),]
    final_data = stats::na.omit(final_data)
    formula = paste("cbind(",paste(metabolites,collapse = ", "),")~ trait - 1 + ",paste("trait:",covariates,collapse = " + "))
    if(type == "dependent"){
      random = stats::as.formula(paste("~us(trait):",individual_ids))
    }else if(type == "independent"){
      random = stats::as.formula(paste("~idh(trait):",individual_ids))
    }
    if(type == "dependent"){
      rcov = ~us(trait):units
    }else if(type == "independent"){
      rcov = ~idh(trait):units
    }

    if(is.null(prior)){
      print("Calculating prior")
      old_model = list()
      for (metabolite in metabolites) {
        sub_model <- parallel::mclapply(1:n_chains, function(i) {
          set.seed(seed + i)
          MCMCglmm::MCMCglmm(stats::as.formula(paste(metabolite, "~",paste(covariates,collapse = " + "))),random = stats::as.formula(paste("~us(1):",individual_ids)),data = final_data,pr=TRUE,nitt = iter/2,thin = thin,burnin = warmup/2)
        }, mc.cores=cores)
        model_summary = list()
        model_summary$Sol = c()
        model_summary$VCV = c()
        for (i in 1:length(sub_model)) {
          model_summary$Sol = rbind(model_summary$Sol,sub_model[[i]]$Sol[,1:(1+length(covariates))])
          model_summary$VCV = rbind(model_summary$VCV,sub_model[[i]]$VCV)
        }
        old_model[[metabolite]] = model_summary
        remove(sub_model)
      }

      prior_mean = c()
      prior_sd = c()
      prior_G = c()
      prior_R = c()
      for (metabolite in metabolites) {
        sub = apply(old_model[[metabolite]]$Sol,2,mean)
        sub = sapply(sub,function(x){ifelse(abs(x) > 1, 10 * ifelse(x/10>1,round(x/10),round(x/10,1)), 1)})
        names(sub) = paste0(metabolite,"_",names(sub))
        prior_mean = c(prior_mean,sub)
        sub = apply(old_model[[metabolite]]$Sol,2,stats::sd)
        sub = sapply(sub,function(x){ifelse(abs(x) > 1, 10 * ifelse(x/10>1,round(x/10),round(x/10,1)), 1)})
        names(sub) = paste0(metabolite,"_",names(sub))
        prior_sd = c(prior_sd,sub)
        sub = apply(old_model[[metabolite]]$VCV,2,mean)[[1]]**0.5
        sub = ifelse(abs(sub) > 1, 10 * ifelse(sub/10>1,round(sub/10),round(sub/10,1)), 1)
        names(sub) = metabolite
        prior_G = c(prior_G,sub)
        sub = apply(old_model[[metabolite]]$VCV,2,mean)[[2]]**0.5
        sub = ifelse(abs(sub) > 1, 10 * ifelse(sub/10>1,round(sub/10),round(sub/10,1)), 1)
        names(sub) = metabolite
        prior_R = c(prior_R,sub)
      }

      prior_mean = prior_mean[order(
        !grepl("Intercept", names(prior_mean)),
        !grepl(covariates[1], names(prior_mean)),# Sort "Intercept" first
        charmatch(gsub(".+_", "", names(prior_mean)), covariates),
        match(gsub("_.+", "", names(prior_mean)), metabolites)  # Sort metabolites last
      )]

      prior_sd = prior_sd[order(
        !grepl("Intercept", names(prior_sd)),
        !grepl(covariates[1], names(prior_sd)),# Sort "Intercept" first
        charmatch(gsub(".+_", "", names(prior_sd)), covariates),
        match(gsub("_.+", "", names(prior_sd)), metabolites)  # Sort metabolites last
      )]

      prior_R <- diag(prior_R)
      off_diag_values <- sample(c(0.1, -0.1), (nrow(prior_R)*(nrow(prior_R)-1))/2, replace = TRUE)
      prior_R[upper.tri(prior_R)] <- off_diag_values
      prior_R <- prior_R + t(prior_R) - diag(diag(prior_R))

      prior_G <- diag(prior_G)
      off_diag_values <- sample(c(0.1, -0.1), (nrow(prior_G)*(nrow(prior_G)-1))/2, replace = TRUE)
      prior_G[upper.tri(prior_G)] <- off_diag_values
      prior_G <- prior_G + t(prior_G) - diag(diag(prior_G))
      prior <- list(B=list(mu = prior_mean,V=diag(prior_sd**2)),R = list(V = prior_R, nu = length(metabolites)*1.1), G = list( G1 = list(V = prior_G, nu = length(metabolites)*1.1)))

    }
    print("Building model")
    model <- parallel::mclapply(1:n_chains, function(i) {
      sub_seed = seed+i-1
      set.seed(sub_seed)
      MCMCglmm::MCMCglmm(stats::as.formula(formula),random = random,rcov = rcov,data = final_data,pr=TRUE,family = rep("gaussian",length(metabolites)),   nitt = iter,thin = thin,burnin = warmup,prior=prior,singular.ok=TRUE)
    }, mc.cores=n_chains)

    summarised <- lapply(model, function(m) m$Sol[,1:((1+length(covariates))*length(metabolites))])
    summarised <- do.call(coda::mcmc.list, summarised)
    t = summary(summarised)
    rf = coda::gelman.diag(summarised)
    co_list = cbind(t$statistics[,1:2],t$quantiles[,c(1,5)],apply(t$quantiles, 1, function(x){ifelse(x[1]*x[5] > 0,TRUE,FALSE) }))
    colnames(co_list) = c("Mean","SD","2.5%","97.5%","Significance")
    co_list = co_list[-c(1:length(metabolites)),]
    rownames(co_list) = gsub("trait","",rownames(co_list))
    rownames(co_list) = gsub(":"," : ",rownames(co_list))
    print("Sampling Posterior predicive distribution")

    doParallel::registerDoParallel(cores)
    div = nrow(final_data)%/% 25
    obs = nrow(final_data)%/% div - 1
    dependent_result<-foreach::foreach(i=0:obs,.combine=rbind,.packages = c('MCMCglmm'))%dopar%{
      if(i==obs){
        drop = (i*div + 1):nrow(final_data)
      }
      else{
        drop = i*div + c(1:div)
      }
      sub_model = update_MCMCglmm(initial_model = model,new_data = final_data[-drop,], iter_reduced = as.numeric(iter)/2,thin = thin,burnin_reduced = as.numeric(warmup)/2,prior = prior,seed = seed,n_chains = n_chains)


      result = prediction(sub_model,newdata = final_data[drop,],levels = cutoff,prior=prior,seed = seed)$summary
      brief_result = cbind(result,unlist(final_data[drop,(length(c(individual_ids,covariates))+2):ncol(final_data)]))
      #         brief_result = cbind(brief_result,(brief_result[,2] < brief_result[,4] & brief_result[,4] < brief_result[,3]))
      brief_result = cbind(brief_result,rep(metabolites,each=length(drop)))
      rownames(brief_result) = rep(paste(final_data[drop,individual_ids],final_data[drop,]$occurrence),length(metabolites))
      brief_result
    }

    rownames(dependent_result) = paste(rownames(dependent_result),dependent_result[,ncol(dependent_result)])
    #test = t(apply(dependent_result, 1, function (x){as.numeric(x[1:(length(x)-1)])+metabolite_means[x[length(x)]]}))
    test = t(apply(dependent_result, 1, function (x){as.numeric(x[1:(length(x)-1)])}))
    colnames(test) = c("fit",as.vector(outer(c("lwr","upr"),cutoff,paste0)),"original")
    for (level in cutoff) {
      test = cbind(test,!(test[,paste0("lwr",level)] < test[,"original"] & test[,"original"] < test[,paste0("upr",level)]))
      colnames(test)=c(colnames(test)[1:ncol(test)-1],paste0("flag",level))
    }
    final = list("model" = model,"type" = type,"significant_covariates" = co_list,"result" = test,"chain_convergence" = rf)
    class(final) = "MetaboVariation"
    return(final)
}

update_MCMCglmm <- function(initial_model, new_data, iter_reduced, thin, burnin_reduced, prior,seed,n_chains) {
  # Extract parameter estimates from the initial model
  initial_matrix = c()
  for (i in 1:n_chains) {
    initial_matrix = rbind(initial_matrix,initial_model[[i]]$VCV)
  }
  initial_matrix = apply(initial_matrix,2,mean)
  initial_G = initial_matrix[1:(length(initial_matrix)/2)]
  initial_R = initial_matrix[(length(initial_matrix)/2+1):length(initial_matrix)]
  n_traits = length(initial_model[[1]]$Residual$original.family)

  if(length(initial_G)==n_traits){
    initial_G <- list(G1 = diag(initial_G))
    initial_R <- list(V = diag(initial_R))
  }
  if(length(initial_G)==n_traits**2){
    initial_G <- list(G1 = matrix(initial_G,nrow = n_traits))
    initial_R <- list(V = matrix(initial_R,nrow = n_traits))
  }
  # Create the list for the start argument
  start_values <- list(G = initial_G, R = initial_R)
  # Update model using new data and initial parameter estimates
  MCMCglmm::MCMCglmm(fixed = initial_model[[1]]$Fixed$formula, random = initial_model[[1]]$Random$formula, rcov = initial_model[[1]]$Residual$formula, data = new_data,
                     pr=TRUE,family = initial_model[[1]]$Residual$original.family, prior = prior, burnin = burnin_reduced, nitt = iter_reduced, thin = thin, start = start_values,singular.ok = TRUE)

}

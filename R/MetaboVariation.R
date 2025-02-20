`%dopar%` <- foreach::`%dopar%`
#' @title
#' Fit a Bayesian generalised linear model to repeated measurements of metabolites.
#' @description
#' Fits a Bayesian generalised linear model (BGLM) to flag individuals with notable variation in their metabolite levels.
#' @details
#' The multivariate BGLM can include covariates and flag the intra-individual variation in metabolite levels from repeated measurements.
#' Posterior predictive distributions of metabolite levels at the individual level are used to flag individuals with observed metabolite levels outside the
#' \code{interval.width}% highest posterior density (HPD) interval at one time point.
#' The function uses the \code{\link[MCMCglmm:MCMCglmm]{MCMCglmm}} function from the \code{\link[MCMCglmm]{MCMCglmm}} package.
#'
#' For the model, we assume the metabolite values of \eqn{M} metabolites from \eqn{N} individuals at \eqn{T} time points, organised into a matrix \eqn{\boldsymbol{Y}} of dimension \eqn{(N \times T) \times M}  are modelled as:
#'  \deqn{\boldsymbol{Y} = \boldsymbol{X} \boldsymbol{\beta} + \boldsymbol{S} + \boldsymbol{\epsilon}}
#'  where the matrix \eqn{\boldsymbol{X}} contains a column of \eqn{1}s for the intercept and the \eqn{L} covariates, the \eqn{(L+1) \times M} matrix \eqn{\boldsymbol{\beta}} represents the regression coefficients and the matrix \eqn{\boldsymbol{S}} (of dimension \eqn{(N \times T) \times M)} captures the random effects of all individuals at each time point and each metabolite. Here, \eqn{\boldsymbol{S} = \boldsymbol{Z}\boldsymbol{u}} where \eqn{\boldsymbol{Z}} is the \eqn{(N \times T) \times (N \times M)} binary design matrix and \eqn{\boldsymbol{u}} contains the random effects for the metabolites. The term \eqn{\boldsymbol{\epsilon}} represents the random error associated with each measurement.
#'  The fixed effects (\eqn{\boldsymbol{\beta}}) follow a multivariate normal prior distribution with mean \eqn{\boldsymbol{\beta}_0} and prior covariance matrix \eqn{\boldsymbol{B}}. For the random effects (\eqn{\boldsymbol{u}}) and residuals (\eqn{\boldsymbol{\epsilon}}), a multivariate normal prior is assumed with means of \eqn{0} along with dense covariance matrices \eqn{\boldsymbol{G}} and \eqn{\boldsymbol{R}}, respectively; \eqn{\boldsymbol{G}} and \eqn{\boldsymbol{R}} have inverse Wishart prior distributions.
#'  The inverse Wishart distribution is characterised by two hyperparameters: the scale matrices \eqn{\boldsymbol{\Sigma}^2} and \eqn{\boldsymbol{\Sigma}^2_\epsilon} for random effects \eqn{\boldsymbol{G}} and \eqn{\boldsymbol{R}}, respectively, and the degrees of freedom \eqn{\nu} for both \eqn{\boldsymbol{G}} and \eqn{\boldsymbol{R}}.
#'  Under dependent settings, the scale matrices \eqn{\boldsymbol{\Sigma}^2} and \eqn{\boldsymbol{\Sigma}^2_\epsilon} are assumed to be dense matrices while under independent setting, they are assumed to be diagonal matrices; details for the specification for \eqn{\boldsymbol{\Sigma}^2}, \eqn{\boldsymbol{\Sigma}^2_\epsilon} and \eqn{\nu} are detailed in the prior argument.
#'
#' @param data A data frame containing data on all variables to be used in the BGLM. Refer to \code{\link{metabol.data}} for the structure of the data frame.
#' @param individual_ids A character string detailing the name of the column in the data that contains the individual IDs.
#' @param metabolites A string or vector of strings containing the names of metabolites to model.
#' @param covariates A vector of strings containing the names of columns in the data that contain the covariates.
#' @param type Determines whether to model with ("dependent") or without ("independent") considering the correlation among metabolites.
#' @param iter The number of iterations per chain (including warmup; defaults to 5000).
#' @param warmup A positive integer specifying the number of warmup iterations. This also specifies the number of iterations used for step size adaptation, so warmup draws should not be used for inference. The value of warmup should not be larger than \code{iter} and the default is
#'   \code{iter/2}.
#' @param thin Thinning rate. Must be a positive integer. Set thin > 1 to save memory and computation time if \code{iter} is large.
#' @param interval.width A numeric vector of desired widths of the highest posterior density interval of the posterior predictive distribution. Values should be less than 1. Default widths are 0.95, 0.975, and 0.99.
#' @param cores The number of cores to use when fitting the models in parallel, which defaults to half of the total cores you have in the system.
# #' @param seed The seed for random number generation to make results reproducible. A seed is default to 19205033.
#' @param n_chains Number of chains to use when fitting the model. Default is 4.
#' @param prior A list of prior hyperparameters with three entries: one for B (fixed effects), one for G (random effects of individuals), and one for R (residual errors of the model).
#'
#' * **B**: This sub-list specifies the hyperparameters of the multivariate normal prior for the fixed effects (regression coefficients). It should contain two elements:
#'     - `mu`: The expected value (mean) of the regression coefficients, assumed to follow a multivariate normal distribution (MVN). By default, `mu` is set to the posterior mean of the regression coefficients obtained from BGLM fitted to each metabolite univaraitely.
#'     - `V`: The covariance matrix representing the strength of belief in `mu`. By default, `V` is set to a diagonal matrix where the diagonal elements correspond to the variances of posterior distributions of the regression coefficients  obtained from BGLM fitted to each metabolite univaraitely.
#'
#' * **G**: This sub-list specifies the hyperparameters of the inverse Wishart prior for the covariance matrix for the random effects of individuals. It should contain two elements:
#'     - `V`: The scale matrix of the inverse-Wishart distribution, which defines the covariance structure of the random effects. By default, the diagonal elements of `V` are set to the variances of the posterior distributions of the random effects as obtained from fitting BGLM to each metabolite univariately. The off-diagonal entries of the `V` are presumed to have an absolute value of `0.1`, and the signs of these off-diagonal terms are fixed to match the signs from the correlation matrix of the observed metabolite data.
#'     - `nu`: The degrees of freedom of the inverse-Wishart distribution. By default, `nu` is set to 150% of the number of metabolites (`M`).
#'
#' * **R**: This sub-list specifies the hyperparameters of the inverse Wishart prior for the covariance matrix for the residual errors of the model. It should contain two elements:
#'     - `V`: The scale matrix of the inverse-Wishart distribution for the residual errors. By default, the diagonal elements of `V` are set to the variances of the posterior distributions of the residuals as obtained from fitting BGLM to each metabolite univariately. The off-diagonal entries of the `V` are presumed to have an absolute value of `0.1`, and the signs of these off-diagonal terms are fixed to match the signs from the correlation matrix of the observed metabolite data
#'     - `nu`: The degrees of freedom of the inverse-Wishart distribution. Like G, `nu` is set to 150% of the number of metabolites (`M`).
#'
#' @return The function returns an object of class \code{\link{MetaboVariation}}. The following values are returned.
#' \itemize{
#' \item result - a data frame containing the summary of the HPDs of the posterior predictive distributions of the individuals in the cohort and a binary flag for each HPD showing whether the individual is flagged under each HPD interval width or not. The data frame contains following columns: the mean, the lower and upper bounds of each HPD interval, the observed value of the metabolite for the individual for that timepoint and the binary flag for each HPD interval. The binary flag indicates whether the observed value lies within the HPD interval or not where 1 denotes the observed value lies outside the interval and 0 denotes the observed value lies in the interval. The row names specify the individual, time point, and metabolite for each observation.
#' \item significant_covariates - a data frame containing the estimates and 95% credible intervals for each regression coefficient for each covariate across all metabolites along with a binary flag indicating whether each covariate is significant in the model for the corresponding metabolite.
#' \item chain_convergence - the potential scale reduction factor, also known as the Gelman-Rubin statistic which measures the extent to which chains are converged. The further the value of the statistic from 1, the poorer the convergence of the chains. The MCMC chains are considered converged if the value lies between 0.9 and 1.05.
#' \item model - contains BGLM model outputs.
#' \item type - shows which model is fitted. It can be either "dependent" or "independent".
#' }
#' @export
#' @references Hadfield, J. D. (2010). MCMC Methods for Multi-Response Generalized Linear Mixed Models: The MCMCglmm R Package. *Journal of Statistical Software*, 33(2), 1–22.
#' @references Gelman, A., & Rubin, D. B. (1992). Inference from Iterative Simulation Using Multiple Sequences. *Statistical Science*, 7(4), 457–472.
#' @seealso \code{\link{radar.plot}}, \code{\link{circos.plot}}, \code{\link{metabolite.heatmap}}, \code{\link{metpair.heatmap}}
#' @examples
#' \dontrun{
#' # Load the simulated data and extract the metabolites names.
#' data(metabol.data)
#' metabolite_list = colnames(metabol.data)[5:length(colnames(metabol.data))]
#' metabolites = get.metabolites(list = metabolite_list)
#' covariates = c("SexM.1F.2","Age","BMI")
#' individual_id = "Individual_id"
#'
#' # Run MetaboVariation on first three metabolites under dependent setting.
#' model = MetaboVariation(data = metabol.data,individual_ids = individual_id,
#' metabolite = metabolites[1:3], covariates = covariates,type="dependent")
#'
#' #' # Run MetaboVariation on first three metabolite under independent setting.
#' model = MetaboVariation(data = metabol.data,individual_ids = individual_id,
#' metabolite = metabolites[1:3], covariates = covariates,type="independent")
#' }
#'
#' @import parallel foreach

MetaboVariation <- function ( data,individual_ids,metabolites,covariates=NULL,type="dependent",iter = 5000,
                                warmup = 2000,thin = 2,interval.width=c(0.95,0.975,0.99),cores = NULL,n_chains=4,prior=NULL){
    # if(!is.null(seed)){
    #   if(!is.numeric(seed)){
    #     stop("seed must be an integer")
    #   }
    # }else{
    #   seed = 19205033
    # }
    # set.seed(seed)
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
    sign_matrix = sign(stats::cor(final_data[,metabolites]))
    if(is.null(covariates)){
      formula = paste("cbind(",paste(metabolites,collapse = ", "),")~ trait - 1 ")

    }else{
      formula = paste("cbind(",paste(metabolites,collapse = ", "),")~ trait - 1 + ",paste("trait:",covariates,collapse = " + "))

    }
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
          # set.seed(seed + i)
          if(is.null(covariates)){
            sub_formula = stats::as.formula(paste(metabolite, "~ 1"))
          }else{
            sub_formula = stats::as.formula(paste(metabolite, "~",paste(covariates,collapse = " + ")))
          }
          MCMCglmm::MCMCglmm(sub_formula,random = stats::as.formula(paste("~us(1):",individual_ids)),data = final_data,pr=TRUE,nitt = iter/2,thin = thin,burnin = warmup/2)
        }, mc.cores=cores)
        model_summary = list()
        model_summary$Sol = c()
        model_summary$VCV = c()
        for (i in 1:length(sub_model)) {
          if(is.null(covariates)){
            model_summary$Sol = c(model_summary$Sol,sub_model[[i]]$Sol[,1])
          }else{
            model_summary$Sol = rbind(model_summary$Sol,sub_model[[i]]$Sol[,1:(1+length(covariates))])
          }
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
        if(is.null(covariates)){
          sub = mean(old_model[[metabolite]]$Sol)
          sub = sapply(sub,function(x){ifelse(abs(x) > 1, 10 * ifelse(x/10>1,round(x/10),round(x/10,1)), 1)})
          names(sub) = paste0(metabolite,"_Intercept")
          prior_mean = c(prior_mean,sub)
          sub = sd(old_model[[metabolite]]$Sol)
          sub = sapply(sub,function(x){ifelse(abs(x) > 1, 10 * ifelse(x/10>1,round(x/10),round(x/10,1)), 1)})
          names(sub) = paste0(metabolite,"_Intercept")
          prior_sd = c(prior_sd,sub)
        }else{
          sub = apply(old_model[[metabolite]]$Sol,2,mean)
          sub = sapply(sub,function(x){ifelse(abs(x) > 1, 10 * ifelse(x/10>1,round(x/10),round(x/10,1)), 1)})
          names(sub) = paste0(metabolite,"_",names(sub))
          prior_mean = c(prior_mean,sub)
          sub = apply(old_model[[metabolite]]$Sol,2,stats::sd)
          sub = sapply(sub,function(x){ifelse(abs(x) > 1, 10 * ifelse(x/10>1,round(x/10),round(x/10,1)), 1)})
          names(sub) = paste0(metabolite,"_",names(sub))
          prior_sd = c(prior_sd,sub)
        }
        sub = apply(old_model[[metabolite]]$VCV,2,mean)[[1]]**0.5
        sub = ifelse(abs(sub) > 1, 10 * ifelse(sub/10>1,round(sub/10),round(sub/10,1)), 1)
        names(sub) = metabolite
        prior_G = c(prior_G,sub)
        sub = apply(old_model[[metabolite]]$VCV,2,mean)[[2]]**0.5
        sub = ifelse(abs(sub) > 1, 10 * ifelse(sub/10>1,round(sub/10),round(sub/10,1)), 1)
        names(sub) = metabolite
        prior_R = c(prior_R,sub)
      }
      if(is.null(covariates)){
        prior_mean = prior_mean[order(names(prior_mean))]

        prior_sd = prior_sd[order(names(prior_mean))]
      }else{
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
      }


      prior_R <- diag(prior_R)
      off_diag_values <- sample(c(0.1), (nrow(prior_R)*(nrow(prior_R)-1))/2, replace = TRUE)
      prior_R[upper.tri(prior_R)] <- off_diag_values
      prior_R <- prior_R + t(prior_R) - diag(diag(prior_R))
      prior_R <- prior_R * sign_matrix


      prior_G <- diag(prior_G)
      off_diag_values <- sample(c(0.1), (nrow(prior_G)*(nrow(prior_G)-1))/2, replace = TRUE)
      prior_G[upper.tri(prior_G)] <- off_diag_values
      prior_G <- prior_G + t(prior_G) - diag(diag(prior_G))
      prior_G <- prior_G * sign_matrix
      prior <- list(B=list(mu = prior_mean,V=diag(prior_sd**2)),R = list(V = prior_R, nu = length(metabolites)*1.5), G = list( G1 = list(V = prior_G, nu = length(metabolites)*1.5)))

    }
    print("Building model")
    model <- parallel::mclapply(1:n_chains, function(i) {
      # sub_seed = seed+i-1
      # set.seed(sub_seed)
      MCMCglmm::MCMCglmm(stats::as.formula(formula),random = random,rcov = rcov,data = final_data,pr=TRUE,family = rep("gaussian",length(metabolites)),   nitt = iter,thin = thin,burnin = warmup,prior=prior,singular.ok=TRUE)
    }, mc.cores=n_chains)

    summarised <- lapply(model, function(m) m$Sol[,1:((1+length(covariates))*length(metabolites))])
    summarised <- do.call(coda::mcmc.list, summarised)
    t = summary(summarised)
    rf = coda::gelman.diag(summarised)
    if(!is.null(covariates)){
      co_list = cbind(t$statistics[,1:2],t$quantiles[,c(1,5)],apply(t$quantiles, 1, function(x){ifelse(x[1]*x[5] > 0,TRUE,FALSE) }))
      colnames(co_list) = c("Mean","SD","2.5%","97.5%","Significance")
      co_list = co_list[-c(1:length(metabolites)),]
      rownames(co_list) = gsub("trait","",rownames(co_list))
      rownames(co_list) = gsub(":"," : ",rownames(co_list))}
    print("Sampling Posterior predicive distribution")

    doParallel::registerDoParallel(cores)
    div = nrow(final_data)%/% 25
    obs = nrow(final_data)%/% div - 1
    dependent_result<-foreach::foreach(i=0:obs,.combine=rbind,.packages = c('MCMCglmm'))%dopar%{
      if(i==obs){
        drop = (i*div + 1):nrow(final_data)
      }else{
        drop = i*div + c(1:div)
      }
      sub_model = update_MCMCglmm(initial_model = model,new_data = final_data[-drop,], iter_reduced = as.numeric(iter)/2,thin = thin,burnin_reduced = as.numeric(warmup)/2,prior = prior,n_chains = n_chains)


      result = prediction(sub_model,newdata = final_data[drop,],levels = interval.width,prior=prior)$summary
      brief_result = cbind(result,unlist(final_data[drop,(length(c(individual_ids,covariates))+2):ncol(final_data)]))
      #         brief_result = cbind(brief_result,(brief_result[,2] < brief_result[,4] & brief_result[,4] < brief_result[,3]))
      brief_result = cbind(brief_result,rep(metabolites,each=length(drop)))
      rownames(brief_result) = rep(paste(final_data[drop,individual_ids],final_data[drop,]$occurrence),length(metabolites))
      brief_result
    }

    rownames(dependent_result) = paste(rownames(dependent_result),dependent_result[,ncol(dependent_result)])
    #test = t(apply(dependent_result, 1, function (x){as.numeric(x[1:(length(x)-1)])+metabolite_means[x[length(x)]]}))
    test = t(apply(dependent_result, 1, function (x){as.numeric(x[1:(length(x)-1)])}))
    colnames(test) = c("mean",as.vector(outer(c("lwr","upr"),interval.width,paste0)),"original")
    for (level in interval.width) {
      test = cbind(test,!(test[,paste0("lwr",level)] < test[,"original"] & test[,"original"] < test[,paste0("upr",level)]))
      colnames(test)=c(colnames(test)[1:ncol(test)-1],paste0("flag",level))
    }
    if(is.null(covariates)){
      final = list("result" = test,"chain_convergence" = rf,"model" = model,"type" = type)

    }else{
      final = list("result" = test,"significant_covariates" = co_list,"chain_convergence" = rf,"model" = model,"type" = type)

    }
    class(final) = "MetaboVariation"
    return(final)
}

update_MCMCglmm <- function(initial_model, new_data, iter_reduced, thin, burnin_reduced, prior,n_chains) {
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

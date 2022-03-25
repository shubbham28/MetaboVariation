`%dopar%` <- foreach::`%dopar%`
#' @title
#' Models the metabolomic profiles to get individual traits
#'@description
#' Fits a generalised linear model on the passed metabolites and provide results for single or multiple metabolites. It works on bayesian framework using the function \code{\link[brms]{brm}} from BRMS package.
#'
#' @param data An object of class data.frame (or one that can be coerced to that class) containing data of all variables used for modelling.
#' @param individual_ids A character string that denotes the column name of individual ids present in the dataset passed.
#' @param metabolite A string or vector containing metabolites name that needs to be modeled.
#' @param covariates A list of column names present in the data that denotes the covariates data for the individuals in the data.
#' @param save_brms_model A logical value indicating whether the model should be saved or not. By default, it is set as FALSE
#' @param full_posterior A logical value indicating whether the full posterior distribution should be included or not. By default, it is set as FALSE and the summarized result is passed by the function.
#' @param iter Number of total iterations per chain (including warmup; defaults to 2000).
#' @param warmup A positive integer specifying number of warmup (aka burnin)
#'   iterations. This also specifies the number of iterations used for stepsize
#'   adaptation, so warmup draws should not be used for inference. The number
#'   of warmup should not be larger than \code{iter} and the default is
#'   \code{iter/2}.
#' @param ... Further arguments are passed to \code{\link[brms]{brm}} for backend processes.
#'
#' @return modelling returns an object of "MetaboVariation" class with "meta.single_model" and "meta.multi_model" subclass. Sublass "meta.multi_model" is a collection of subclass "meta.single_model" objects.
#' And subclass "meta.single_model" contains following values.
#' \itemize{
#' \item significant_covariates - a list of covariates that influences the metabolite value significantly for the individuals. Based on the posterior distribution of the covariates, if zero lies in the 90% credible interval, it is assumed that particular covariate has no significant realtion with the metabolite value.
#' \item result - It stores a summarized result for the posterior predictive distribution of the cohort. It contains the mean, 95% credible interval range, tails of the credible interval and the original value of the individual for that timepoint. It also has a flag that shows whether the original values lies within the credible interval or not.
#' \item iterations - Number of total iterations per chain used for modelling (including warmup; defaults to 2000).
#' \item warmup - A positive integer specifying number of warmup (aka burnin) iterations. This also specifies the number of iterations used for stepsize adaptation, so warmup draws should not be used for inference. The number of warmup should not be larger than iter and the default is iter/2.
#' \item Rhat - Rhat refers to the potential scale reduction statistic, also known as the Gelman-Rubin statistic.
#' Thus, it measures the extent to which chains are converging. The further the value of the statistic from 1, the worse. Model is considered good if the value lies between 0.9 and 1.05.
#' \item metabolite - Contains the name of metabolite on which modelling is done.
#' \item full_posterior - It contains the full posterior predictive distiribution of the data. It can be used for further analysis as per the user requirements.
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' data(metabol.data)
#' metabolite_list = colnames(metabol.data)[5:length(colnames(metabol.data))]
#' metabolites = get.metabolites(list = metabolite_list)
#' covariates = c("SexM.1F.2","Age","BMI")
#' individual_id = "Sample.Id"
#' model = MetaboVariation(data = metabol.data,individual_ids = individual_id,
#' metabolite = metabolites[1], covariates = covariates)
#' model = MetaboVariation(data = metabol.data,individual_ids = individual_id,
#' metabolite = metabolites[1:3], full_posterior = TRUE)
#' }
#'
#' @import parallel foreach stats

MetaboVariation <- function ( data,individual_ids,metabolite,covariates=NULL, save_brms_model = FALSE,full_posterior = FALSE,iter = 2000,
          warmup = floor(iter/2),...){

  future = getOption("future", TRUE)
  file = paste(metabolite,"model",sep="_")
  seed = 19205033
  i <- NULL
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    cores <- getOption("mc.cores",2)
  } else {
    # use all cores in devtools::test()
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
  if(length(metabolite)==1){
    var_list = na.omit(stringr::str_extract(col_list,stringr::regex(paste(".*",metabolite,".*",sep = ""))))
    if(length(var_list)==0){
      stop("Metabolite passed is not present in dataset")
    }
    new_data = reshape2::melt(data,id.vars = c(individual_ids,covariates),measure.vars = var_list,variable.name = "occurrence",value.name = "measurement")
    new_data = na.omit(new_data)
    print(paste("Model building for",metabolite))
    model = brms::brm(formula,data = new_data,iter = iter,
                      warmup = warmup, future = future, seed = seed, file = file)
    rows <- nrow(new_data)
    doParallel::registerDoParallel(cores)
    print(paste("Calculating posterior predictive distribution for",metabolite))
    result<-foreach::foreach(i=1:rows,.combine=cbind,.packages = c('brms','rstan','future'))%dopar%{
      new_model = update(model,file=file,newdata = new_data[-i,],silent=2,refresh=0)
      rbind.data.frame(paste(new_data[i,individual_ids], new_data$occurrence[i]),brms::posterior_predict(new_model,newdata = new_data[i,]))
    }
    if(!save_brms_model){
      file.remove(paste(file,".rds",sep=""))
    }
    colnames(result) = result[1,]
    result = result[-1,]
    result = apply(result,2,as.numeric)
    brief_result = t(apply(result, 2, function(x) {c(mean(x),quantile(x,probs = c(0.025,0.975)),quantile(x,0.975)-quantile(x,0.025))}))
    brief_result = cbind(brief_result,new_data$measurement)
    brief_result = cbind(brief_result,(brief_result[,2] < brief_result[,5] & brief_result[,5] < brief_result[,3]))
    colnames(brief_result) = c("mean","2.5%","97.5%","QR","original","identifier")
    co_list = names(na.omit(apply(summary(model)$fixed[2:nrow(summary(model)$fixed),3:4],1,function(x) {ifelse(prod(x[1],x[2])>=0,1,NA )})))
    if(identical(co_list,character(0))){
      co_list = NULL
    }

    return_list = list("significant_covariates" = co_list,"result" = brief_result,"warmup" = warmup,"iterations"=iter,
                       "Rhat" = summary(model)$spec_pars$Rhat)
    if(full_posterior){
      return_list[["full_posterior"]] = result
    }
    gc()
    remove(model)
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
                          warmup = warmup, future = future, seed = seed, file = file)

        rows <- nrow(new_data)
        doParallel::registerDoParallel(cores)
        print(paste("Calculating posterior predictive distribution for",met))
        result<-foreach::foreach(i=1:rows,.combine=cbind,.packages = c('brms','rstan','future'))%dopar%{
          new_model = update(model,file=file,newdata = new_data[-i,],silent=2,refresh=0)
          rbind.data.frame(paste(new_data[i,individual_ids], new_data$occurrence[i]),brms::posterior_predict(new_model,newdata = new_data[i,]))
        }
        if(!save_brms_model){
          file.remove(paste(file,".rds",sep=""))
        }
        colnames(result) = result[1,]
        result = result[-1,]
        result = apply(result,2,as.numeric)
        brief_result = t(apply(result, 2, function(x) {c(mean(x),quantile(x,probs = c(0.025,0.975)),quantile(x,0.975)-quantile(x,0.025))}))
        brief_result = cbind(brief_result,new_data$measurement)
        brief_result = cbind(brief_result,(brief_result[,2] < brief_result[,5] & brief_result[,5] < brief_result[,3]))
        colnames(brief_result) = c("mean","2.5%","97.5%","QR","original","identifier")
        co_list = names(na.omit(apply(summary(model)$fixed[2:nrow(summary(model)$fixed),3:4],1,function(x) {ifelse(prod(x[1],x[2])>=0,1,NA )})))
        if(identical(co_list,character(0))){
          co_list = NULL
        }
        final_list[[met]] = list("metabolite" = met,"significant_covariates" = co_list,"result" = brief_result,"warmup" = warmup,"iterations"=iter,
                                 "Rhat" = summary(model)$spec_pars$Rhat)
        if(full_posterior){
          final_list[[met]][["full_posterior"]] = result
        }
        class(final_list[[met]]) = c("MetaboVariation","meta.single_model")
        gc()
        remove(model)

      }
      class(final_list) = c("MetaboVariation","meta.multi_model")
      return(final_list)

    }
}


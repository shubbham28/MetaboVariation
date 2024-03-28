#' `%>%` <- magrittr::`%>%`
#' `%like%` <- data.table::`%like%`
#' globalVariables(c("metabolite","flag","original","timepoint","lower","upper","nonzero_count","n","total_count","percentage","value"))
#'
#' #' @title Plot Flagged Individuals for Metabolites
#' #'
#' #' @description
#' #' This function generates plots to visualize the results obtained after fitting a MetaboVariation Model to metabolomics data.
#' #' Two types of plots can be produced, each offering a different perspective on the results.
#' #'
#' #' @param x An object of class \code{\link{MetaboVariation}} containing the fitted model results.
#' #' @param type A character specifying the type of plot to generate: either \code{radarplot} or \code{heatmap}.
#' #' @param interval The interval for the Highest Posterior Distribution (HPD) that you want to visualize.
#' #' @param cutoff The cutoff for the heatmap, based on timepoints. This parameter determines the threshold for including individuals based on the percentage of timepoints they are flagged in.
#' #' It ranges from 0 to 100.
#' #' @param threshold A positive integer indicating the minimum number of metabolites an individual must be flagged across different metabolites to be included in the heatmap.
#' #' @param individual Used with \code{type = "radarplot"} to visualize data for a specific individual.
#' #' @param ... Additional arguments to be passed to the plotting functions.
#' #'
#' #' @return This function will show the results of modelling in the following two ways:
#' #' \itemize{
#' #'   \item Radar plot - The radar plot provides an overview of an individual's metabolic profile across different timepoints. Individual's metabolic values depicted as a colored line while the remaining cohort are presented as gray lines. The shaded area indicates the posterior predictive HPD interval for the individual.
#' #'   The metabolite where individual is flagged is marked by asterisk.
#' #'   \item Heatmap plot: The heatmap illustrates the individuals flagged in the metabolites, with the shades of color indicating the frequency of flags across timepoints in the \code{\link{MetaboVariation}} object.
#' #' }
#' #' @seealso \code{\link{MetaboVariation}}
#' #'
#' #' @examples
#' #' \dontrun{
#' #' data(metabol.data)
#' #' metabolite_list = colnames(metabol.data)[5:length(colnames(metabol.data))]
#' #' metabolites = get.metabolites(list = metabolite_list)
#' #' covariates = c("SexM.1F.2","Age","BMI")
#' #' individual_id = "Individual_id"
#' #' model = MetaboVariation(data = metabol.data,individual_ids = individual_id,
#' #' metabolite = metabolites[1:3], covariates = covariates,full_posterior = TRUE)
#' #' plot(x = model, type = "circos")
#' #' plot(x = model, type = "metabolites_count",timepoints = c(1,2,4))
#' #'}
#' #'
#' #' @export
#' #' @aliases plot
#' #' @method plot MetaboVariation
#' #'
#' #'
#' plot.MetaboVariation <- function(x,type = NULL,interval = NULL, cutoff = NULL, threshold = 1, individual = NULL, ...){
#'   if(!inherits(x,"MetaboVariation")){
#'     stop("Model passed is not of class 'Metabovariation'")
#'   }
#'   if(type == "radarplot"){
#'     test=x$result
#'     plot_df = test[rownames(test) %like% individual,]
#'     plot_df = plot_df[,c(paste0(c("lwr","upr"),interval),"original",paste0("flag",interval))]
#'     plot_df = data.frame(plot_df)
#'     colnames(plot_df) = c("lower","upper","value","flag")
#'     plot_df$timepoint = gsub(".* (.*) .*","\\1",rownames(plot_df))
#'     plot_df$metabolite = gsub(".* .* (.*)","\\1",rownames(plot_df))
#'     original_df = test[!rownames(test) %like% individual,"original"]
#'     original_df = data.frame(original_df)
#'     colnames(original_df) = c("original")
#'     original_df$timepoint = gsub(".* (.*) .*","\\1",rownames(original_df))
#'     original_df$individual = gsub("(.*) .* .*","\\1",rownames(original_df))
#'     original_df$metabolite = gsub(".* .* (.*)","\\1",rownames(original_df))
#'     original_df = original_df[original_df$timepoint %in% unique(plot_df$timepoint),]
#'     for (met in unique(plot_df$metabolite)) {
#'       subdf = plot_df[plot_df$metabolite == met,]
#'       s = c(unlist(subdf[,1:3]),original_df[original_df$metabolite == met,1])
#'       # s = sapply(unlist(subdf[,1:3]), function(x){(x-min(unlist(subdf[,1:3])))/(max(unlist(subdf[,1:3]))-min(unlist(subdf[,1:3])))})
#'       # s = sapply(s, function(x){x = x-mean(s)})
#'       # subdf[,1:3] = matrix(sapply(s, function(x){(1-0)/(1+1)*(x+1)+0}),nrow = 4)
#'       s = sapply(s, function(x){(x-min(s))/(max(s)-min(s))})
#'       s = sapply(s, function(x){x = x-mean(s)})
#'       s = sapply(s, function(x){(1-0)/(1+1)*(x+1)+0})
#'       rows = nrow(plot_df[plot_df$metabolite == met,])
#'       subdf[,1:3] = matrix(s[1:(3*rows)],nrow = rows)
#'       plot_df[plot_df$metabolite == met,] = subdf
#'       original_df[original_df$metabolite == met,1] = s[-(1:(3*rows))]
#'     }
#'
#'     flagged_metabolites <- plot_df %>%
#'       dplyr::group_by(metabolite) %>%
#'       dplyr::filter(any(flag == 1)) %>%
#'       dplyr::distinct(metabolite)
#'     flagged_metabolites = unique(c(unlist(flagged_metabolites),metabolites))
#'     #plot_df = plot_df[plot_df$metabolite %in% flagged_metabolites,]
#'     #original_df = original_df[original_df$metabolite %in% flagged_metabolites,]
#'
#'     metabolites = unique(plot_df$metabolite)
#'     rows = nrow(plot_df)
#'     plot_df = rbind(plot_df,plot_df[plot_df$metabolite == metabolites[1],])
#'     plot_df[(rows+1):nrow(plot_df),"metabolite"] = "Z"
#'     rows = nrow(original_df)
#'     original_df = rbind(original_df,original_df[original_df$metabolite == metabolites[1],])
#'     original_df[(rows+1):nrow(original_df),"metabolite"] = "Z"
#'
#'     ggplot2::ggplot()+
#'       ggplot2::geom_line(data = original_df,ggplot2::aes(x=metabolite,y=original,group=individual),color="lightgray",alpha = 0.75,show.legend = FALSE)+
#'       ggplot2::geom_line(data = plot_df,ggplot2::aes(x=metabolite,y=value,group=timepoint,color=timepoint),show.legend = FALSE)+
#'       ggplot2::geom_line(data = plot_df,ggplot2::aes(x=metabolite,y=lower,group=timepoint,color=timepoint),linetype = "dashed",show.legend = FALSE)+
#'       ggplot2::geom_line(data = plot_df,ggplot2::aes(x=metabolite,y=upper,group=timepoint,color=timepoint),linetype = "dashed",show.legend = FALSE)+
#'       ggplot2::geom_ribbon(data = plot_df,ggplot2::aes(x=metabolite,y=value,group=timepoint,ymin=lower,ymax=upper,fill=timepoint),alpha=0.25,show.legend = FALSE)+
#'       ggplot2::geom_text(data = plot_df,ggplot2::aes(x=metabolite,y=1,label=ifelse(metabolite!="Z",ifelse(flag==1,paste0(metabolite,"*"),metabolite),""),color = timepoint,fontface = ifelse(flag==1,"bold.italic","plain"),family = ifelse(flag==1,"serif","mono")),show.legend = FALSE)+
#'       ggplot2::theme_light()+
#'       ggplot2::coord_polar()+
#'       ggplot2::labs(x="",y="")+
#'       ggplot2::scale_y_continuous(breaks = NULL)+
#'       ggplot2::facet_wrap(~timepoint)+ggplot2::scale_x_discrete(breaks = metabolites,labels=rep("",length(metabolites)),expand = c(0,0))+
#'       ggplot2::theme(text = ggplot2::element_text(size = 20),panel.grid.major = ggplot2::element_line(linewidth = 0))
#'   }
#'
#'   if(type == "heatmap"){
#'     test= x$result
#'     cutoff = as.numeric(cutoff)/100
#'     subdf = test[,paste0("flag",interval)]
#'     subdf = data.frame(subdf)
#'     colnames(subdf) = c("flag")
#'     subdf$individual = gsub("(.*) .* .*","\\1",rownames(subdf))
#'     subdf$timepoint = gsub(".* (.*) .*","\\1",rownames(subdf))
#'     subdf$metabolite = gsub(".* .* (.*)","\\1",rownames(subdf))
#'     #subdf$flag = sapply(subdf$flag,function(x){ifelse(x==0,1,0)})
#'     timepoints=length(unique(subdf$timepoint))
#'     metabolite_flagged <- stats::aggregate(flag ~ individual + metabolite, data = subdf, FUN = sum)
#'     #metabolite_flagged$flag = sapply(metabolite_flagged$flag,function(x){ifelse(x>=(timepoints*as.numeric(type)),"Yes","No")})
#'     individuals = unique(metabolite_flagged[metabolite_flagged$flag>=(timepoints*as.numeric(cutoff)),"individual"])
#'     metabolite_flagged$flag = factor(metabolite_flagged$flag)
#'     metabolite_flagged = metabolite_flagged[metabolite_flagged$individual %in% individuals,]
#'     flagged_metabolites <- metabolite_flagged %>%
#'       dplyr::group_by(metabolite) %>%
#'       dplyr::summarize(nonzero_count = sum(flag != 0), total_count = dplyr::n()) %>%
#'       dplyr::mutate(percentage = nonzero_count / total_count) %>%
#'       dplyr::filter(percentage > 0.1) %>%
#'       dplyr::distinct(metabolite)
#'     metabolite_flagged = metabolite_flagged[metabolite_flagged$metabolite %in% flagged_metabolites,]
#'     flagged_individuals <- metabolite_flagged %>%
#'       dplyr::group_by(individual) %>%
#'       dplyr::summarize(nonzero_count = sum(flag != 0), total_count = dplyr::n()) %>%
#'       dplyr:: mutate(percentage = nonzero_count / total_count) %>%
#'       dplyr::filter(percentage > threshold/100) %>%
#'       dplyr::distinct(individual)
#'     metabolite_flagged = metabolite_flagged[metabolite_flagged$individual %in% unlist(flagged_individuals),]
#'
#'
#'     color <- grDevices::colorRampPalette(c("white", "#2B7CE9"))(length(levels(metabolite_flagged$flag))+1)
#'     names(color) = c("",levels(metabolite_flagged$flag))
#'     ggplot2::ggplot(metabolite_flagged,ggplot2::aes(x=metabolite,y=individual,fill=flag))+ggplot2::geom_tile(color = "white")+
#'       ggplot2::scale_fill_manual(values = color)+ggplot2::theme_minimal()+ggplot2::theme(text = ggplot2::element_text(size = 20),axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))  + ggplot2::labs(x = "Metabolites",y="Individuals",fill = "Flag")
#'
#'   }
#'
#' }

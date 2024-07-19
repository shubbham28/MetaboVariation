`%>%` <- magrittr::`%>%`
`%like%` <- data.table::`%like%`
globalVariables(c("metabolite","flag","original","timepoint","lower","upper","value","x"))


#' @title Radar plot for an individual.
#' @description
#' Visualise a radar plot showing the metabolic profile of an individual across all time points and for all metabolites.
#'
#' @param model An object of class \code{\link{MetaboVariation}} containing the fitted model results.
#' @param interval The interval for the Highest Posterior Distribution (HPD) that you want to visualise.
#' @param individual The ID for individual you want to visualise.
#' @param title Title to plot
#'
#' @return Returns a radar plot for a specific individual at the specified interval.
#' @export
#'
#' @examples
#' \dontrun{
#' # Load the simulated data and extract the metabolites names.
#' data(metabol.data)
#' metabolite_list = colnames(metabol.data)[5:length(colnames(metabol.data))]
#' metabolites = get.metabolites(list = metabolite_list)
#' covariates = c("SexM.1F.2","Age","BMI")
#' individual_id = "Individual_id"
#'
#' # Run the MetaboVariation.
#' model = MetaboVariation(data = metabol.data,individual_ids = individual_id,
#' metabolite = metabolites[1:3], covariates = covariates,cutoff=c(0.95,0.975,0.99))
#'
#' # Visualise the radar plot.
#' radarplor(model,interval=0.975,individual = 3)
#'}

radar.plot <- function(model,interval,individual,title=""){
  if(!inherits(model,"MetaboVariation")){
    stop("Model passed is not of class 'Metabovariation'")
  }
  # if(model$type!="dependent"){
  #   stop("Function is only used with dependent model")
  # }
  test=model$result
  plot_df = test[rownames(test) %like% paste0("^",individual," "),]
  plot_df = plot_df[,c(paste0(c("lwr","upr"),interval),"original",paste0("flag",interval))]
  plot_df = data.frame(plot_df)
  colnames(plot_df) = c("lower","upper","value","flag")
  plot_df$timepoint = gsub(".* (.*) .*","\\1",rownames(plot_df))
  plot_df$metabolite = gsub(".* .* (.*)","\\1",rownames(plot_df))
  original_df = test[!rownames(test) %like% individual,"original"]
  original_df = data.frame(original_df)
  colnames(original_df) = c("original")
  original_df$timepoint = gsub(".* (.*) .*","\\1",rownames(original_df))
  original_df$individual = gsub("(.*) .* .*","\\1",rownames(original_df))
  original_df$metabolite = gsub(".* .* (.*)","\\1",rownames(original_df))
  original_df$timepoint = factor(original_df$timepoint)
  plot_df$timepoint = factor(plot_df$timepoint,levels = levels(original_df$timepoint))
  original_df = original_df[original_df$timepoint %in% unique(plot_df$timepoint),]
  for (met in unique(plot_df$metabolite)) {
    subdf = plot_df[plot_df$metabolite == met,]
    s = c(unlist(subdf[,1:3]),original_df[original_df$metabolite == met,1])
    # s = sapply(unlist(subdf[,1:3]), function(x){(x-min(unlist(subdf[,1:3])))/(max(unlist(subdf[,1:3]))-min(unlist(subdf[,1:3])))})
    # s = sapply(s, function(x){x = x-mean(s)})
    # subdf[,1:3] = matrix(sapply(s, function(x){(1-0)/(1+1)*(x+1)+0}),nrow = 4)
    s = sapply(s, function(x){(x-min(s))/(max(s)-min(s))})
    s = sapply(s, function(x){x = x-mean(s)})
    s = sapply(s, function(x){(1-0)/(1+1)*(x+1)+0})
    rows = nrow(plot_df[plot_df$metabolite == met,])
    subdf[,1:3] = matrix(s[1:(3*rows)],nrow = rows)
    plot_df[plot_df$metabolite == met,] = subdf
    original_df[original_df$metabolite == met,1] = s[-(1:(3*rows))]
  }

  # flagged_metabolites <- plot_df %>%
  #   group_by(metabolite) %>%
  #   filter(any(flag == 1)) %>%
  #   distinct(metabolite)
  #plot_df = plot_df[plot_df$metabolite %in% flagged_metabolites,]
  #original_df = original_df[original_df$metabolite %in% flagged_metabolites,]

  metabolites = unique(plot_df$metabolite)
  rows = nrow(plot_df)
  plot_df = rbind(plot_df,plot_df[plot_df$metabolite == metabolites[1],])
  plot_df[(rows+1):nrow(plot_df),"metabolite"] = "Z"
  rows = nrow(original_df)
  original_df = rbind(original_df,original_df[original_df$metabolite == metabolites[1],])
  original_df[(rows+1):nrow(original_df),"metabolite"] = "Z"

  ggplot2::ggplot()+
    ggplot2::geom_line(data = original_df,ggplot2::aes(x=metabolite,y=original,group=individual),color="lightgray",alpha = 0.75,show.legend = FALSE)+
    ggplot2::geom_line(data = plot_df,ggplot2::aes(x=metabolite,y=value,group=timepoint,color=timepoint),show.legend = FALSE)+
    ggplot2::geom_line(data = plot_df,ggplot2::aes(x=metabolite,y=lower,group=timepoint,color=timepoint),linetype = "dashed",show.legend = FALSE)+
    ggplot2::geom_line(data = plot_df,ggplot2::aes(x=metabolite,y=upper,group=timepoint,color=timepoint),linetype = "dashed",show.legend = FALSE)+
    ggplot2::geom_ribbon(data = plot_df,ggplot2::aes(x=metabolite,y=value,group=timepoint,ymin=lower,ymax=upper,fill=timepoint),alpha=0.25,show.legend = FALSE)+
    ggplot2::geom_text(data = plot_df,ggplot2::aes(x=metabolite,y=1,label=ifelse(metabolite!="Z",ifelse(flag==1,paste0(metabolite,"*"),metabolite),""),color = timepoint,fontface = ifelse(flag==1,"bold.italic","plain"),family = ifelse(flag==1,"serif","mono")),show.legend = FALSE)+
    ggplot2::theme_light()+
    ggplot2::coord_polar()+
    ggplot2::labs(x="",y="",title = title)+
    ggplot2::scale_y_continuous(breaks = NULL)+
    ggplot2::facet_wrap(~timepoint,drop=FALSE)+ggplot2::scale_x_discrete(breaks = metabolites,labels=rep("",length(metabolites)),expand = c(0,0))+ggplot2::theme(text = ggplot2::element_text(size = 20),panel.grid.major = ggplot2::element_line(linewidth = 0))
  #labels = function(x){ifelse(plot_df$flag==1,as.expression(x),as.expression(bquote(bold(.(x)))))}
}

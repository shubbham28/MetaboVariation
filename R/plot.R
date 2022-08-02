`%>%` <- magrittr::`%>%`
globalVariables(c("timepoint","value"))
unit = grid::unit

#' @title
#' Plots the flagged individuals for the metabolites.
#' @description
#' A function to plot the results after fitting a Bayesian generalised linear model (BGLM) model to data.
#' Two types of plots can be produced depicting results from different perspectives.
#'
#' @param x An object of class \code{\link{MetaboVariation}}
#' @param type A character specifying the plot type (eg: "circos", "metabolites_count"). By default, "circos" is shown.
#' @param timepoints A numeric vector containing timepoints to be shown in the plots. By default, it is NULL and will show all the timepoints.
#' @param ... Further arguments are ignored.
#'
#' @return plot will show the results of modelling in the following two ways:
#' \itemize{
#'   \item Circos plot - A plot showing the posterior predictive credible interval for all the individuals across timepoints. If the observed value
#'   is outside the credible interval, that bar will be flagged as black and the individual's id will be shown in __bold__.
#'   \item Metabolites count plot - A bar plot that shows the number of individuals that are flagged in each timepoint for all metabolites in the \code{\link{MetaboVariation}} object.
#' }
#' @seealso \code{\link{MetaboVariation}}

#'
#' @examples
#' \dontrun{
#' data(metabol.data)
#' metabolite_list = colnames(metabol.data)[5:length(colnames(metabol.data))]
#' metabolites = get.metabolites(list = metabolite_list)
#' covariates = c("SexM.1F.2","Age","BMI")
#' individual_id = "Individual_id"
#' model = MetaboVariation(data = metabol.data,individual_ids = individual_id,
#' metabolite = metabolites[1:3], covariates = covariates,full_posterior = TRUE)
#' plot(x = model, type = "circos")
#' plot(x = model, type = "metabolites_count",timepoints = c(1,2,4))
#'}
#'@aliases plot
#' @method plot MetaboVariation
#' @export
#'
plot.MetaboVariation <- function(x,type = "circos",timepoints = NULL,...){
  model = x
  if(!inherits(model,"MetaboVariation")){
    stop("Model passed is not of class 'Metabovariation'")
  }
  colors = grDevices::topo.colors(50)[seq(1, 50, 5)]
  if(!(type %in% c("circos","metabolites_count"))){
    stop("Invalid plot type")
  }

   #Single metabolite model
  if(inherits(model,"meta.single_model")){

    #Circos plot
    if(type == "circos"){
    values = cbind.data.frame(stringr::str_split_fixed(rownames(model$result),pattern=" ",n=2),model$result[,4])
    status = cbind.data.frame(stringr::str_split_fixed(rownames(model$result),pattern=" ",n=2),model$result[,6])
    colnames(values) = c("id","timepoint","value")
    colnames(status) = c("id","timepoint","status")
    values$timepoint = sapply(values$timepoint,function(x){as.numeric(gsub("[^0-9]", "",unique(x)))})
    time_data = as.numeric(gsub("[^0-9]", "",unique(values$timepoint)))
    if(!is.null(timepoints)){
      if(!all(timepoints %in% time_data)){
        stop("Timepoints passed are not present in the data")
      }
    }
    values = tidyr::spread(values,timepoint,value)
    status = tidyr::spread(status,timepoint,status)
    size = cbind.data.frame(stringr::str_split_fixed(rownames(model$result),pattern=" ",n=2),model$result[,6])
    colnames(size) = c("id","timepoint","status")
    size = tidyr::spread(size,timepoint,status)
    size[is.na(size)] = 1
    size = apply(size[,-1],1,function(x){sum(x)%/%length(time_data)})
    outliers = cbind.data.frame(stringr::str_split_fixed(rownames(model$result),pattern=" ",n=2),"identifier" = model$result[,6])

    circlize::circos.clear()
    circlize::circos.par(cell.padding=c(0,0,0,0),gap.degree=1,points.overflow.warning=F)
    circlize::circos.initialize(letters[1], xlim = c(0, nrow(values)))
    circlize::circos.track(ylim=c(0,1),bg.border=NA,track.height=.2,track.margin=c(.01,0),
                           panel.fun=function(x,y)for(i in 1:nrow(values))
                             circlize::circos.text(i,0,values$id[i],adj=c(0,.5),facing="clockwise",niceFacing=T,cex = ifelse(size[i] == 0, 0.75,0.5),font=ifelse(size[i] == 0, 4,1)))
    if(is.null(timepoints)){
      height = 0.75/ncol(values)
      for (i in 2:ncol(values)) {
        status[,i] = sapply(status[,i], function(x){ifelse(x==0,1,colors[i])})
        circlize::circos.track(ylim = c(0,1),bg.border=NA,track.margin=c(0,0),track.height=height, panel.fun = function(x, y) {
          value = scales::rescale(values[,i])
          circlize::circos.barplot(value, 1:nrow(values),border = "white",bar_width = 0.75, col = status[,i])
        })
        lgd = ComplexHeatmap::Legend(at = paste0("Timepoint",colnames(values)[2:ncol(values)]), type = "lines", legend_gp = grid::gpar(col = colors[2:ncol(values)], lwd = 3), title_position = "topleft", title = "Timepoints", ncol = 4)
        ComplexHeatmap::draw(lgd, x = unit(0.2, "npc"), y = unit(0.05, "npc"), just = c("left", "bottom"))
        ComplexHeatmap::draw(ComplexHeatmap::Legend(at = paste(model$metabolite,"plot")), x = unit(0.425, "npc"), y = unit(0.9, "npc"), just = c("left", "bottom"))
      }
    }
    else if(is.numeric(timepoints) & max(timepoints) <= ncol(values)-1){

      height = 0.6/length(timepoints)
      for (i in timepoints) {
        i = i+1
        status[,i] = sapply(status[,i], function(x){ifelse(x==0,1,colors[i])})
        circlize::circos.track(ylim = c(0,1),bg.border=NA,track.margin=c(0,0),track.height=height, panel.fun = function(x, y) {
          value = scales::rescale(values[,i])
          circlize::circos.barplot(value, 1:nrow(values),border = "white",bar_width = 0.75, col = status[,i])
        })
        lgd = ComplexHeatmap::Legend(at = paste0("Timepoint",timepoints), type = "lines", legend_gp = grid::gpar(col = colors, lwd = 3), title_position = "topleft", title = "Timepoints", ncol = 4)
        ComplexHeatmap::draw(lgd, x = unit(0.2, "npc"), y = unit(0.05, "npc"), just = c("left", "bottom"))
        ComplexHeatmap::draw(ComplexHeatmap::Legend(at = paste(model$metabolite,"plot")), x = unit(0.425, "npc"), y = unit(0.9, "npc"), just = c("left", "bottom"))
      }


    }
    }#End of Circos plot

     # Metabolite count plot
    if(type == "metabolites_count"){
      individuals = c()
      status = cbind.data.frame(stringr::str_split_fixed(rownames(model$result),pattern=" ",n=2),model$result[,6])
      colnames(status) = c("id","timepoint","status")
      time_data = as.numeric(gsub("[^0-9]", "",unique(status$timepoint)))
      status = tidyr::spread(status,timepoint,status)
      individuals = cbind(individuals,apply(!status[,-1],2,function(x){ifelse(x==FALSE | is.na(x),0,1)}))
      sum_m = apply(individuals,2,sum)
      if(!is.null(timepoints)){
        if(!all(timepoints %in% time_data)){
          stop("Timepoints passed are not present in the data")
        }
      }
      if(is.null(timepoints)){
      sum_plot = data.frame(sum_m)
      rownames(sum_plot) = paste0("Timepoint ",1:length(sum_m))
      }
      else if(is.numeric(timepoints) & max(timepoints) <= ncol(status)-1){
        sum_plot = data.frame(sum_m[timepoints])
        sum_plot$Timepoint = paste0("Timepoint ",timepoints)
        colnames(sum_plot) = c("sum_m","Timepoint")
      }
      sum_plot$Timepoint = rownames(sum_plot)
      plotly::plot_ly(sum_plot,x = ~Timepoint, y = ~sum_m, type = "bar",marker = list(color=colors[1],width = 0.2)) %>%
        plotly::layout(yaxis = list(title = 'No. of individuals flagged',showgrid = FALSE,showline = TRUE))
    }    # End Metabolite count plot

  }

  # Multiple metabolite model
  else if(inherits(model,"meta.multi_model")) {

    #Circos plot
    if(type == "circos"){
      for (met in names(model)) {
      values = cbind.data.frame(stringr::str_split_fixed(rownames(model[[met]]$result),pattern=" ",n=2),model[[met]]$result[,4])
      status = cbind.data.frame(stringr::str_split_fixed(rownames(model[[met]]$result),pattern=" ",n=2),model[[met]]$result[,6])
      colnames(values) = c("id","timepoint","value")
      colnames(status) = c("id","timepoint","status")
      values$timepoint = sapply(values$timepoint,function(x){as.numeric(gsub("[^0-9]", "",unique(x)))})
      time_data = as.numeric(gsub("[^0-9]", "",unique(values$timepoint)))
      if(!is.null(timepoints)){
        if(!all(timepoints %in% time_data)){
          stop("Timepoints passed are not present in the data")
        }
      }

      values = tidyr::spread(values,timepoint,value)
      status = tidyr::spread(status,timepoint,status)
      size = cbind.data.frame(stringr::str_split_fixed(rownames(model[[met]]$result),pattern=" ",n=2),model[[met]]$result[,6])
      colnames(size) = c("id","timepoint","status")
      size = tidyr::spread(size,timepoint,status)
      size[is.na(size)] = 1
      size = apply(size[,-1],1,function(x){sum(x)%/%length(time_data)})
      outliers = cbind.data.frame(stringr::str_split_fixed(rownames(model[[met]]$result),pattern=" ",n=2),"identifier" = model[[met]]$result[,6])

      circlize::circos.clear()
      circlize::circos.par(cell.padding=c(0,0,0,0),gap.degree=1,points.overflow.warning=F)
      circlize::circos.initialize(letters[1], xlim = c(0, nrow(values)))
      circlize::circos.track(ylim=c(0,1),bg.border=NA,track.height=.2,track.margin=c(.01,0),
                             panel.fun=function(x,y)for(i in 1:nrow(values))
                               circlize::circos.text(i,0,values$id[i],adj=c(0,.5),facing="clockwise",niceFacing=T,cex = ifelse(size[i] == 0, 0.75,0.5),font=ifelse(size[i] == 0, 4,1)))
      if(is.null(timepoints)){
        height = 0.75/ncol(values)
        for (i in 2:ncol(values)) {
          status[,i] = sapply(status[,i], function(x){ifelse(x==0,1,colors[i])})
          circlize::circos.track(ylim = c(0,1),bg.border=NA,track.margin=c(0,0),track.height=height, panel.fun = function(x, y) {
            value = scales::rescale(values[,i])
            circlize::circos.barplot(value, 1:nrow(values),border = "white",bar_width = 0.75, col = status[,i])
          })
          lgd = ComplexHeatmap::Legend(at = paste0("Timepoint",colnames(values)[2:ncol(values)]), type = "lines", legend_gp = grid::gpar(col = colors[2:ncol(values)], lwd = 3), title_position = "topleft", title = "Timepoints", ncol = 4)
          ComplexHeatmap::draw(lgd, x = unit(0.2, "npc"), y = unit(0.05, "npc"), just = c("left", "bottom"))
          ComplexHeatmap::draw(ComplexHeatmap::Legend(at = paste(met,"plot")), x = unit(0.425, "npc"), y = unit(0.9, "npc"), just = c("left", "bottom"))
        }
      }
      else if(is.numeric(timepoints) & max(timepoints) <= ncol(values)-1){

        height = 0.6/length(timepoints)
        for (i in timepoints) {
          i = i+1
          status[,i] = sapply(status[,i], function(x){ifelse(x==0,1,colors[i])})
          circlize::circos.track(ylim = c(0,1),bg.border=NA,track.margin=c(0,0),track.height=height, panel.fun = function(x, y) {
            value = scales::rescale(values[,i])
            circlize::circos.barplot(value, 1:nrow(values),border = "white",bar_width = 0.75, col = status[,i])
          })
          lgd = ComplexHeatmap::Legend(at = paste0("Timepoint",timepoints), type = "lines", legend_gp = grid::gpar(col = colors, lwd = 3), title_position = "topleft", title = "Timepoints", ncol = 4)
          ComplexHeatmap::draw(lgd, x = unit(0.2, "npc"), y = unit(0.05, "npc"), just = c("left", "bottom"))
          ComplexHeatmap::draw(ComplexHeatmap::Legend(at = paste(met,"plot")), x = unit(0.425, "npc"), y = unit(0.9, "npc"), just = c("left", "bottom"))
        }

      }


    }#End of for loop
    }#End of Circos plot

    # Metabolite count plot
    if(type == "metabolites_count"){
      individuals = c()
      for (met in names(model)) {
        status = cbind.data.frame(stringr::str_split_fixed(rownames(model[[met]]$result),pattern=" ",n=2),model[[met]]$result[,6])
        colnames(status) = c("id","timepoint","status")
        time_data = as.numeric(gsub("[^0-9]", "",unique(status$timepoint)))
        if(!is.null(timepoints)){
          if(!all(timepoints %in% time_data)){
            stop("Timepoints passed are not present in the data")
          }
        }
        status = tidyr::spread(status,timepoint,status)
        individuals = cbind(individuals,apply(!status[,-1],2,function(x){ifelse(x==FALSE | is.na(x),0,1)}))
      }

      sum_m = apply(individuals,2,sum)
      sum_plot = c()
      if(is.null(timepoints)){
        timepoints = length(sum_m)%/%length(names(model))
        timepoints = 1:timepoints
      }
      for (met in names(model)) {
        sum_plot = rbind.data.frame(sum_plot,c(model[[met]]$metabolite,sum_m[paste(model[[met]]$metabolite,"_",timepoints,sep = "")]))

      }
      colnames(sum_plot) = c("Metabolites",paste0("Timepoint.",timepoints))
      for (i in 2:ncol(sum_plot)) {
        sum_plot[,i] = as.numeric(sum_plot[,i])
      }
    #  sum_plot = reshape2::melt(sum_plot,id.vars = "Metabolites", variable.name = "Timepoints",value.name = "count")
      #return(sum_plot)
     # figure = ggplot2::ggplot(sum_plot, ggplot2::aes(fill=Timepoints, y=count, x=Metabolites)) +  ggplot2::geom_bar(position="dodge", stat="identity")
      #print(figure)
      #plotly::ggplotly(figure)
     figure = plotly::plot_ly(x = sum_plot[["Metabolites"]],y = sum_plot[[paste0("Timepoint.",timepoints[1])]], type = "bar",marker = list(color=colors[timepoints[1]]),width = 0.2,name = paste0("Timepoint ",timepoints[1])) %>%
      plotly::layout(yaxis = list(title = 'No. of individuals flagged',showgrid = FALSE,showline = TRUE,range = c(min(sum_plot[,-1]) - 1,max(sum_plot[,-1]) + 1)), barmode = 'group')
      for (i in 2:length(timepoints)) {
      figure = figure %>% plotly::add_trace(y = sum_plot[[paste0("Timepoint.",timepoints[i])]],marker = list(color=colors[i]),width = 0.2,name = paste0("Timepoint ",timepoints[i]))
      }
      return(plotly::config(plotly::as_widget(figure)))
    }    # End Metabolite count plot

  }


}

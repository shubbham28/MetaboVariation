`%>%` <- magrittr::`%>%`
`%like%` <- data.table::`%like%`
globalVariables(c("metabolite","flag","original","timepoint","lower","upper","value","flag"))
unit = grid::unit

#' @title Plots the metabolite-specific circos plot.
#' @description
#' Visualise a circos plot showing a highest posterior density (HPD) interval for posterior predictive distributions for all individuals across all time points for a specific metabolite.
#'
#' @param model An object of class \code{\link{MetaboVariation}} containing the fitted model results.
#' @param metabolite A string containing the name of the metabolite you want to visualise.
#' @param interval.width The width of the highest posterior density (HPD) interval considered. Must be a numeric between 0 and 1 with default value of 0.95.
#' @param title Text for the title of the plot.
#'
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
#' # Run MetaboVariation on first three metabolites.
#' model = MetaboVariation(data = metabol.data,individual_ids = individual_id,
#' metabolite = metabolites[1:3], covariates = covariates,interval.width=c(0.95,0.975,0.99))
#'
#' # Plots a circos plot for first metabolite with an interval width of 0.975.
#' circos.plot(model,metabolite = metabolites[1],interval.width=0.975)
#'}

circos.plot <- function(model,metabolite,interval.width = 0.95,title = ""){
  if(!inherits(model,"MetaboVariation")){
    stop("Model passed is not of class 'Metabovariation'")
  }
  colors = grDevices::topo.colors(50)[seq(1, 50, 5)]

    data = cbind.data.frame(stringr::str_split_fixed(rownames(model$result),pattern=" ",n=3),model$result[,colnames(model$result)[colnames(model$result) %like% interval.width]])
    colnames(data) = c("id","timepoint","metabolite","lower","upper","flag")
    data = data[data$metabolite %like% metabolite,]
    values = cbind.data.frame(data$id,data$timepoint,data$upper - data$lower)
    status = cbind.data.frame(data$id,data$timepoint,data$flag)
    colnames(values) = c("id","timepoint","value")
    colnames(status) = c("id","timepoint","status")
    values$timepoint = sapply(values$timepoint,function(x){as.numeric(gsub("[^0-9]", "",unique(x)))})
    time_data = as.numeric(gsub("[^0-9]", "",unique(values$timepoint)))

    values = tidyr::spread(values,timepoint,value)
    status = tidyr::spread(status,timepoint,status)
    size = status
    size[is.na(size)] = 0
    size = apply(size[,-1],1,function(x){(length(time_data)-sum(x))%/%length(time_data)})
    outliers = status

    circlize::circos.clear()
    circlize::circos.par(cell.padding=c(0,0,0,0),gap.degree=1,points.overflow.warning=F,start.degree = 90)
    circlize::circos.initialize(letters[1], xlim = c(0, nrow(values)))
    circlize::circos.track(ylim=c(0,1),bg.border=NA,track.height=.2,track.margin=c(.01,0),
                           panel.fun=function(x,y)for(i in 1:nrow(values))
                             circlize::circos.text(i,0,values$id[i],adj=c(0,.5),facing="clockwise",niceFacing=T,cex = ifelse(size[i] == 0, 0.75,0.5),font=ifelse(size[i] == 0, 4,1)))
      height = 0.75/ncol(values)
      for (i in 2:ncol(values)) {
        status[,i] = sapply(status[,i], function(x){ifelse(x==1,1,colors[i])})
        circlize::circos.track(ylim = c(0,1),bg.border=NA,track.margin=c(0,0),track.height=height, panel.fun = function(x, y) {
          limit = stats::quantile(values[,i],c(0.25,0.75),na.rm = TRUE) + c(-1.5,+1.5)*stats::IQR(values[,i],na.rm = TRUE)
          value = values[,i]
          # value[which(value < limit[1])] = limit[1]
          # value[which(value > limit[2])] = limit[2]
          value = scales::rescale(value)
          circlize::circos.barplot(value, 1:nrow(values),border = "white",bar_width = 0.75, col = status[,i])
        })
      }
        lgd = ComplexHeatmap::Legend(at = paste0("",colnames(values)[2:ncol(values)]), type = "lines", legend_gp = grid::gpar(col = colors[2:ncol(values)], lwd = 3), title_position = "topleft", title = "Timepoints", ncol = 4)
        ComplexHeatmap::draw(lgd, x = unit(0.2, "npc"), y = unit(0.05, "npc"), just = c("left", "bottom"))
        if(title!=""){
          ComplexHeatmap::draw(ComplexHeatmap::Legend(at = title), x = unit(0.425, "npc"), y = unit(0.9, "npc"), just = c("left", "bottom"))
        }


    }#End of Circos plot

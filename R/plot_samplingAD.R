

#' plot_samplingAD
#'
#' Plots the set of feasible combinations of individuals and observations per individual needed to obtain a statistical power given an effect size, alpha level and metric of interest. Based on Pirla, Taquet and Quoidbach (2021).
#' @param power Minimum required power. Input as a number from 0.01 to 0.99. Defaults to 0.8 (80\% power).
#' @param metric Metric of interest. Must be a single character from "Average", "Rel.SD", "SD", "RMSSD", "TKEO", "PAC" or "Autocorrelation".
#' @param r Number from 0.01 to 0.99 indicating the expected effect size (Pearson correlation) of interest.
#' @param p.value Alpha level. Must be one of the following numbers: 0.01, 0.05, 0.001, 0.005, 0.001. Defaults to 0.05.
#' @keywords Affect Dynamics, Power Analysis
#' @export
#' @examples
#' plot_samplingAD(power=0.8, metric="SD", r=0.1, p.value=0.01)
#' @references Pirla, Taquet and Quoidbach (2021). ADD REFERENCE



plot_samplingAD<-function(power=0.8, metric, r, p.value=0.05 ){

  #Warning Messages/Errors:
  if(missing(metric) | sum(data.powerAD$Metric ==metric)==0 | length(metric)>1 |is.character(metric)==F ){stop("Please enter a single valid metric (as string character) from the following ones: Average, Rel.SD, SD, RMSSD, TKEO, PAC, Autocorrelation")}
  if(r<0.01 |r>0.99 |is.numeric(r)==F){stop("Please, provide a valid effect size (r) between 0.01 and 0.99.")}
  if(sum(p.value==c(0.1,0.05,0.01,0.005,0.001))==0){stop("Please, provide a valid alpha level (p.value). Valid alpha levels are 0.1, 0.05, 0.01, 0.005, 0.001.")}
  if(power<0 |power>1 |is.numeric(power)==F){stop("Please, provide a valid required power between 0.01 and 0.99.")}

  #Obtain Valid Observations from df
  r<-round(r, digits = 2)
  selected<-data.powerAD[data.powerAD$Metric==metric & data.powerAD$r==r & data.powerAD$P.value==p.value,]
  selected<-selected[selected$Power>=power,]


  ind.selected<-unique(selected$Individuals)
  selected2<-as.data.frame(matrix(ncol = 6))
  names(selected2)<-c("Power", "Individuals", "Samples", "r", "P.value", "Metric" )
  for (i in ind.selected) {
    row.sel<-selected[selected$Individuals==i & selected$Samples==min(selected$Samples[selected$Individuals==i]),]
    selected2<-rbind(selected2,row.sel)
  }
  selected3<-selected2[-1,]
  rownames(selected3) <- NULL


  if(nrow(selected3)==0){warning("No plot created. No feasible combinations of individuals and samples for this combination of effect size, metric, power and significance level." )
  }else{
    #Create Plot
    graph<-ggplot2::ggplot(selected3, ggplot2::aes(x=Individuals, y=Samples)) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::ylim(0,50) +
      ggplot2::ggtitle("Feasible Combinations") +
      ggplot2::scale_x_continuous(trans = "log10", limits = c(9,5200), breaks = c(10,20,40,80,160,320,640,1280,2560,5120))

    print(graph)
    invisible(graph)
  }


}

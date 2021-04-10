
#' samplingAD
#'
#' Estimates the set of feasible combinations of individuals and observations per individual needed to obtain a statistical power given an effect size, alpha level and metric of interest. Based on Pirla, Taquet and Quoidbach (2021).
#' @param power Minimum required power. Input as a number from 0.01 to 0.99. Defaults to 0.8 (80\% power).
#' @param metric Metric of interest. Must be a single character from "Average", "Rel.SD", "SD", "RMSSD", "TKEO", "PAC" or "Autocorrelation".
#' @param r Number from 0.01 to 0.99 indicating the expected effect size (Pearson correlation) of interest.
#' @param p.value Alpha level. Must be one of the following numbers: 0.01, 0.05, 0.001, 0.005, 0.001. Defaults to 0.05.
#' @keywords Affect Dynamics, Power Analysis
#' @export
#' @examples
#' samplingAD(power=0.8, metric="SD", r=0.1, p.value=0.01)
#' @references Pirla, Taquet and Quoidbach (2021). ADD REFERENCE

samplingAD<-function(power=0.8,  metric, r, p.value=0.05){

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
  selected2<-selected2[-1,]
  sm<-unique(selected2$Samples)
  selected3<-as.data.frame(matrix(ncol = 6))
  names(selected3)<-c("Power", "Individuals", "Samples", "r", "P.value", "Metric" )
  for (i in sm) {
    row.sel<-selected2[selected2$Samples==i & selected2$Individuals==min(selected2$Individuals[selected2$Samples==i]),]
    selected3<-rbind(selected3,row.sel)
  }

  selected3<-selected3[-1,]
  rownames(selected3) <- NULL

  if(nrow(selected3)==0){warning("No feasible combinations of individuals and samples for this combination of effect size, metric, power and significance level.")
  }else{

    #Print Output:
    if(metric=="Average"){metric.me<-"trait affect."}
    if(metric=="Rel.SD"){metric.me<-"the Relative Standard Deviation (Rel. SD) in affect."}
    if(metric=="SD"){metric.me<-"the Standard Deviation (SD) in affect."}
    if(metric=="RMSSD"){metric.me<-"the Root Mean Square of Successive Differences (RMSSD) in affect."}
    if(metric=="TKEO"){metric.me<-"the Teager-Kaiser Energy Operator (TKEO) of affect."}
    if(metric=="PAC"){metric.me<-"the Probability of Acute Change (PAC) in affect."}
    if(metric=="Autocorrelation"){metric.me<-"the Autocorrelation coefficient of affect."}


    cat("\nFeasible combinations of individuals and samples per individual to obtain a power of", power*100, "%  to detect an effect of size r=", r, "using an alpha of", p.value,
        "when interested in",metric.me , "\n\nFeasible sampling combinations:\n")
    cat("\n")
    print(selected3[,1:3], row.names=F)
    cat("\n -----------------------\n")
  }
  #Invisible output
  invisible(selected3)
}

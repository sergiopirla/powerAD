
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

  if(sum(selected$Power>=power)>0){
    sel.b<-data.frame(Individuals=c(NA),Samples=c(NA),Power=c(NA) )
    for (i in 1:10) {
      s<-5*i
      if(nrow(selected[selected$Power>=power & selected$Samples==s,])>0){
        ind.high<-min(selected$Individuals[selected$Power>=power & selected$Samples==s])
        ind.low<-max(selected$Individuals[selected$Power<=power& selected$Samples==s])

        selected.b<-selected[selected$Individuals %in% c(ind.low, ind.high) & selected$Samples==s,]
        selected.b<-unique(selected.b)
        rownames(selected.b) <- NULL

        power.gain.ind<- (mean(selected.b$Power[selected.b$Individuals==ind.high]) - mean(selected.b$Power[selected.b$Individuals==ind.low]))/(ind.high - ind.low)
        if(is.nan(power.gain.ind)){power.gain.ind<-0}
        Individuals<-seq(from=ind.low, to=ind.high, by=10)
        da<-data.frame(Individuals)
        da$Samples<-s
        Power.baseline<-selected.b$Power[selected.b$Individuals==ind.low & selected.b$Samples==s]
        da$Power<- Power.baseline + (da$Individuals - ind.low)*power.gain.ind
        da<-da[da$Individuals==min(da$Individuals[da$Power>=power]),]
        sel.b[i,]<-da
      }else{
        da<-data.frame(Individuals=NA)
        da$Samples<-s
        da$Power<- NA
        sel.b[i,]<-da
      }}

    sel.b<-na.omit(sel.b)

    row.n<-nrow(sel.b)
    for(i in 2:row.n){
      a<-as.numeric(sel.b$Individuals[(row.n - i +2)])
      b<-as.numeric(sel.b$Individuals[1:((row.n - i +2)-1)])
      if(sum(a>b)>0 | sum(a==b)>0){sel.b<-sel.b[-(row.n - i +2),]}
    }

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

    selected3$Power<-round(selected3$Power,2)
    selected3$Samples<-round(selected3$Samples,0)
    selected3$Individuals<-round(selected3$Individuals,0)
    sel.b$Power<-round(sel.b$Power,2)
    sel.b<-sel.b[,c(3,1,2)]

    if(metric=="Average"){metric.me<-"trait affect."}
    if(metric=="Rel.SD"){metric.me<-"the Relative Standard Deviation (Rel. SD) in affect."}
    if(metric=="SD"){metric.me<-"the Standard Deviation (SD) in affect."}
    if(metric=="RMSSD"){metric.me<-"the Root Mean Square of Successive Differences (RMSSD) in affect."}
    if(metric=="TKEO"){metric.me<-"the Teager-Kaiser Energy Operator (TKEO) of affect."}
    if(metric=="PAC"){metric.me<-"the Probability of Acute Change (PAC) in affect."}
    if(metric=="Autocorrelation"){metric.me<-"the Autocorrelation coefficient of affect."}

    cat("\nNumber of individuals and affect reports per individual (samples) needed to achieve a power of", power*100, "% or more to detect a Pearson correlation of size r =", r, "between", metric.me, "and a given variable using a two-tailed t-test and an alpha of", p.value,".\n")
    cat("\n")
    print(sel.b, row.names=F)
    cat("\n -----------------------\n")
    cat("Power is estimated through a linear interpolation using the sample combinations included in our main analyses.")
    cat("We refrain from making power extrapolations and therefore, only consider sampling approaches that range between 10 and 5120 participants and from 5 to 50 affect reports per participant (samples). The following table presents the minimal sampling combinations included in our main analyses that yielded the specified power: \n")
    cat("\n")
    print(selected3[,1:3], row.names=F)
    cat("\n -----------------------\nHow to report:\n")
    cat("\n")
    cat("Power analysis for affect dynamic studies (Pirla et al., 2021) suggests that our sampling strategy achieved a statistical power of", power*100, "% to detect a Pearson correlation of size r =",r,"using a two-tailed t-test with an alpha of",p.value,".")
    cat("\n -----------------------\nReference:\n")
    cat("\n")
    cat("Pirla, S., Taquet, M., & Quoidbach, J. (2021). Measuring Affect Dynamics: An Empirical Framework. https://doi.org/10.31219/osf.io/x2ywa")



    #Invisible output
    invisible(sel.b)

  }else{
    warning("No feasible combinations of individuals and samples for this combination of effect size, metric, power and significance level.")}
}

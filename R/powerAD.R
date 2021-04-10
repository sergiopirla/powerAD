
#' powerAD
#'
#' Estimates the statistical power of an affect dynamics time series study for a given sampling, effect size, alpha level and metric of interest. Based on Pirla, Taquet and Quoidbach (2021).
#' @param individuals Number of individuals sampled.
#' @param samples Number of affect observations per individual.
#' @param metric Metric of interest. Must be a single character from "Average", "Rel.SD", "SD", "RMSSD", "TKEO", "PAC" or "Autocorrelation".
#' @param r Number from 0.01 to 0.99 indicating the expected effect size (Pearson correlation) of interest.
#' @param p.value Alpha level. Must be one of the following numbers: 0.01, 0.05, 0.001, 0.005, 0.001. Defaults to 0.05.
#' @keywords Affect Dynamics, Power Analysis
#' @export
#' @examples
#' powerAD(individuals=500, samples=10, metric="SD", r=0.1, p.value=0.01)
#' @references Pirla, Taquet and Quoidbach (2021). ADD REFERENCE

powerAD<-function(individuals, samples, metric, r, p.value=0.05){

  #Warning Messages/Errors:
  if(missing(individuals) | individuals<0 | length(individuals)>1 | is.numeric(individuals)==F){stop("Please enter a valid sample size (number of individuals)")}
  if(missing(samples) | samples<0 | length(samples)>1 | is.numeric(samples)==F){stop("Please enter a valid number of observations per participant")}
  if(missing(metric) | sum(data.powerAD$Metric ==metric)==0 | length(metric)>1 |is.character(metric)==F ){stop("Please enter a single valid metric (as string character) from the following ones: Average, Rel.SD, SD, RMSSD, TKEO, PAC, Autocorrelation")}
  if(individuals>5120){warning("Number of individuals entered exceeds the maximum. Power computed for maximum number of individuals (max=5120)")}
  if(samples>50){warning("Number of observations per participant exceeds the maximum. Power computed for maximum number of samples per individual (max=50)")}
  if(individuals<10){warning("Number of individuals lower than the minimum. Power computed for minimum number of individuals (min=10)")}
  if(samples<5){warning("Number of observations per participant lower than the minimum. Power computed for minimum number of samples per individual (min=5)")}
  if(r<0.01 |r>0.99|is.numeric(r)==F){stop("Please, provide a valid effect size (r) between 0.01 and 0.99.")}
  if(sum(p.value==c(0.1,0.05,0.01,0.005,0.001))==0){stop("Please, provide a valid alpha level (p.value). Valid alpha levels are 0.1, 0.05, 0.01, 0.005, 0.001.")}

  #Obtain Valid Observations from df

  if(individuals>5120){ind.high<-Inf}else{ind.high<-min(data.powerAD$Individuals[data.powerAD$Individuals>=individuals])}
  if(individuals<10){ind.low<-Inf}else{ind.low<-max(data.powerAD$Individuals[data.powerAD$Individuals<=individuals])}
  if(samples<5){s.low<-Inf}else{s.low<-max(data.powerAD$Samples[data.powerAD$Samples<=samples])}
  if(samples>50){s.high<-Inf}else{s.high<-min(data.powerAD$Samples[data.powerAD$Samples>=samples])}
  r<-round(r, digits = 2)
  selected<-data.powerAD[data.powerAD$Individuals %in% c(ind.low, ind.high) & data.powerAD$Samples %in% c(s.low, s.high) & data.powerAD$Metric==metric & data.powerAD$r==r & data.powerAD$P.value==p.value,]
  selected<-unique(selected)
  rownames(selected) <- NULL

  #Compute linear aproximation of Power:
  power.gain.sample<- (mean(selected$Power[selected$Samples==s.high]) - mean(selected$Power[selected$Samples==s.low]))/5
  power.gain.ind<- (mean(selected$Power[selected$Individuals==ind.high]) - mean(selected$Power[selected$Individuals==ind.low]))/(ind.high - ind.low)
  power.gain.sample<-power.gain.sample * (samples - min(selected$Samples))
  power.gain.ind<-power.gain.ind * (individuals - min(selected$Individuals))
  if(is.nan(power.gain.sample)){power.gain.sample<-0}
  if(is.nan(power.gain.ind)){power.gain.ind<-0}
  aprox.Power<- (selected$Power[selected$Individuals==min(selected$Individuals) & selected$Samples==min(selected$Samples)] + power.gain.ind + power.gain.sample )

  #Print Output:
  if(metric=="Average"){metric.me<-"trait affect."}
  if(metric=="Rel.SD"){metric.me<-"the Relative Standard Deviation (Rel. SD) in affect."}
  if(metric=="SD"){metric.me<-"the Standard Deviation (SD) in affect."}
  if(metric=="RMSSD"){metric.me<-"the Root Mean Square of Successive Differences (RMSSD) in affect."}
  if(metric=="TKEO"){metric.me<-"the Teager-Kaiser Energy Operator (TKEO) of affect."}
  if(metric=="PAC"){metric.me<-"the Probability of Acute Change (PAC) in affect."}
  if(metric=="Autocorrelation"){metric.me<-"the Autocorrelation coefficient of affect."}


  cat("\nPower to detect an effect of size r=", r, "using an alpha of", p.value,
      "when interested in",metric.me , "\n\nEmpirical Power (closest sampling approaches):\n")
  cat("\n")
  print(selected[,1:3], row.names=F)
  cat("\n -----------------------\nLinear aproximation of Power under the sampling approach specified (",
      individuals, "individuals and", samples, "observations per individual):\n")
  cat("Aprox. Power=", aprox.Power*100, "%")

  #Invisible output
  inv<-list(selected,aprox.Power)
  names(inv)<-c("Empirical.Power", "Linear.Aproximation")
  invisible(inv)
}

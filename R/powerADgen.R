
#' powerADgen
#'
#' Estimates the statistical power of an affect dynamics time series study for a given sampling, effect size, alpha level and metric of interest.
#' @param individuals Number of individuals sampled.
#' @param samples Number of affect observations per individual.
#' @param metric Metric of interest. Must be a single character from "Average", "Rel.SD", "SD", "RMSSD", "TKEO", "PAC" or "Autocorrelation".
#' @param r Number from 0.01 to 0.99 indicating the expected effect size (Pearson correlation) of interest.
#' @param p.value Alpha level. Must be one of the following numbers: 0.01, 0.05, 0.001, 0.005, 0.001. Defaults to 0.05.
#' @param min Affect scale lower bound.
#' @param max Affect scale upper bound.
#' @param data Data frame containing affect observations used perfrom power analysis.
#' @param id Name of variable containing id observations.
#' @param affect Name of variable containing affect observations.
#' @keywords Affect Dynamics, Power Analysis
#' @export
#' @examples
#' powerADgen(individuals=500, samples=10, metric="SD", r=0.1, p.value=0.01)
#' @references Pirla, Taquet and Quoidbach (2021). ADD REFERENCE

powerADgen<-function(individuals, samples, metric, r, p.value=0.05, min, max, data,id,affect, time, perm=100 ){

  #Error Messages
  if(missing(individuals) | individuals<0 | length(individuals)>1 | is.numeric(individuals)==F){stop("Please enter a valid sample size (number of individuals)")}
  if(missing(samples) | samples<0 | length(samples)>1 | is.numeric(samples)==F){stop("Please enter a valid number of observations per participant")}
  if(missing(metric) | length(metric)>1 |is.character(metric)==F ){stop("Please enter a single valid metric (as string character) from the following ones: Average, Rel.SD, SD, RMSSD, TKEO, PAC, Autocorrelation")}
  if(missing(id) | id %in% names(data) ==F ){stop("Please enter a valid id variable name.")}
  if(missing(time) | id %in% names(data) ==F ){stop("Please enter a valid time variable name.")}
  if(missing(affect) | affect %in% names(data) ==F ){stop("Please enter a valid affect variable name.")}
  if(metric=="Rel.SD" & (missing(min) | missing(max))){stop("Please enter the affect scale boundaries.")}
  if(r<0.01 |r>0.99|is.numeric(r)==F){stop("Please, provide a valid effect size (r) between 0.01 and 0.99.")}
  if(sum(p.value==c(0.1,0.05,0.01,0.005,0.001))==0){stop("Please, provide a valid alpha level (p.value). Valid alpha levels are 0.1, 0.05, 0.01, 0.005, 0.001.")}
  names(data)[which(names(data)==id)]<-"id"
  names(data)[which(names(data)==time)]<-"time"
  names(data)[which(names(data)==affect)]<-"Happiness"
  data<-data[,c("id", "time", "Happiness")]
  data<-na.omit(data)
  data<-dplyr::group_by(data, id)
  data<-dplyr::mutate(data, n=dplyr::n())
  data<-data[data$n>samples,]
  if(individuals>=length(unique(data$id))){stop("Combination of individuals and samples not valid. Not enough observations in data to simulate power analysis.")}


  #Note at this point, the dataframe names of interest are fixed (id, time, Happiness)



  #Average
  ###############################################################################

  if(metric=="Average"){

    #Estimate Measure all datapoints
    da <- dplyr::group_by(data, id)
    da <- dplyr::mutate(da, av.happ=mean(Happiness))
    da <- dplyr::slice_head(da)

    res <- c()

    for (i in 1:perm) {

      #Simulate variable
      da$sim <- faux::rnorm_pre(da$av.happ, mu = 0, sd = 1, r = r)

      #Resample individuals
      rand.id <- sample(da$id, size=individuals)
      data.res <- data[data$id %in% rand.id,]

      #Resample observations per individual
      data.res <- dplyr::group_by(data.res, id)
      data.res <- dplyr::slice_sample(data.res, n=samples)

      #Re-estimate measure
      data.res <- dplyr::mutate(data.res, av.happ2=mean(Happiness))
      data.res <- dplyr::slice_head(data.res)
      data.res <- data.res[,c("id", "av.happ2")]

      #Merge
      data.res <- merge(data.res, da, by="id")

      #Check if significant
      a <- stats::cor.test(data.res$av.happ2, data.res$sim)
      res[i]<-ifelse(a$estimate>0 & a$p.value<p.value,1,0)

    }}



  #SD
  ###############################################################################

  if(metric=="SD"){

    #Estimate Measure all datapoints
    da <- dplyr::group_by(data, id)
    da <- dplyr::mutate(da, sd.happ=sd(Happiness))
    da <- dplyr::slice_head(da)

    res <- c()

    for (i in 1:perm) {

      #Simulate variable
      da$sim <- faux::rnorm_pre(da$sd.happ, mu = 0, sd = 1, r = r)

      #Resample individuals
      rand.id <- sample(da$id, size=individuals)
      data.res <- data[data$id %in% rand.id,]

      #Resample observations per individual
      data.res <- dplyr::group_by(data.res, id)
      data.res <- dplyr::slice_sample(data.res, n=samples)

      #Re-estimate measure
      data.res <- dplyr::mutate(data.res, sd.happ2=sd(Happiness))
      data.res <- dplyr::slice_head(data.res)
      data.res <- data.res[,c("id", "sd.happ2")]

      #Merge
      data.res <- merge(data.res, da, by="id")

      #Check if significant
      a <- stats::cor.test(data.res$sd.happ2, data.res$sim)
      res[i]<-ifelse(a$estimate>0 & a$p.value<p.value,1,0)

    }}


  #Rel.SD
  ###############################################################################

  if(metric=="Rel.SD"){

    #Estimate Measure all datapoints
    da <- dplyr::group_by(data, id)
    da <- dplyr::mutate(da, relsd.happ=relativeVariability::relativeSD(Happiness, MIN = min, MAX=max))
    da <- dplyr::slice_head(da)

    res <- c()

    for (i in 1:perm) {

      #Simulate variable
      da$sim <- faux::rnorm_pre(da$relsd.happ, mu = 0, sd = 1, r = r)

      #Resample individuals
      rand.id <- sample(da$id, size=individuals)
      data.res <- data[data$id %in% rand.id,]

      #Resample observations per individual
      data.res <- dplyr::group_by(data.res, id)
      data.res <- dplyr::slice_sample(data.res, n=samples)

      #Re-estimate measure
      data.res <- dplyr::mutate(data.res, relsd.happ2=relativeVariability::relativeSD(Happiness, MIN = min, MAX=max))
      data.res <- dplyr::slice_head(data.res)
      data.res <- data.res[,c("id", "relsd.happ2")]

      #Merge
      data.res <- merge(data.res, da, by="id")

      #Check if significant
      a <- stats::cor.test(data.res$relsd.happ2, data.res$sim)
      res[i]<-ifelse(a$estimate>0 & a$p.value<p.value,1,0)

    }}


  #RMSSD
  ###############################################################################

  if(metric=="RMSSD"){

    #Estimate Measure all datapoints
    da <- dplyr::group_by(data, id)
    da <- dplyr::arrange(da,time, by_group=T)
    da <- dplyr::mutate(da, diff = Happiness - dplyr::lag(Happiness))
    da$diff2<-da$diff*da$diff
    da <- dplyr::mutate(da, rmssd.happ = sqrt(mean(diff2, na.rm=T)))
    da <- dplyr::slice_head(da)

    res <- c()

    for (i in 1:perm) {

      #Simulate variable
      da$sim <- faux::rnorm_pre(da$rmssd.happ, mu = 0, sd = 1, r = r)

      #Resample individuals
      rand.id <- sample(da$id, size=individuals)
      data.res <- data[data$id %in% rand.id,]

      #Resample observations per individual
      data.res <- dplyr::group_by(data.res, id)
      data.res <- dplyr::slice_sample(data.res, n=samples)

      #Re-estimate measure

      data.res <- dplyr::arrange(data.res,time, by_group=T)
      data.res <- dplyr::mutate(data.res, diff = Happiness - dplyr::lag(Happiness))
      data.res$diff2<-data.res$diff*data.res$diff
      data.res <- dplyr::mutate(data.res, rmssd.happ2 = sqrt(mean(diff2, na.rm=T)))
      data.res <- dplyr::slice_head(data.res)
      data.res <- data.res[,c("id", "rmssd.happ2")]

      #Merge
      data.res <- merge(data.res, da, by="id")

      #Check if significant
      a <- stats::cor.test(data.res$rmssd.happ2, data.res$sim)
      res[i]<-ifelse(a$estimate>0 & a$p.value<p.value,1,0)

    }}



  #TKEO
  ###############################################################################

  if(metric=="TKEO"){

    #Estimate Measure all datapoints
    da <- dplyr::group_by(data, id)
    da <- dplyr::arrange(da,time, by_group=T)
    da <- dplyr::mutate(da, Happiness.lag = dplyr::lag(Happiness))
    da <- dplyr::mutate(da, Happiness.lead = dplyr::lead(Happiness))
    da$tkeo <- (da$Happiness*da$Happiness) - (da$Happiness.lag * da$Happiness.lead)
    da <- dplyr::mutate(da, tkeo.happ = sum(tkeo, na.rm = T)/dplyr::n())
    da <- dplyr::slice_head(da)

    res <- c()

    for (i in 1:perm) {

      #Simulate variable
      da$sim <- faux::rnorm_pre(da$tkeo.happ, mu = 0, sd = 1, r = r)

      #Resample individuals
      rand.id <- sample(da$id, size=individuals)
      data.res <- data[data$id %in% rand.id,]

      #Resample observations per individual
      data.res <- dplyr::group_by(data.res, id)
      data.res <- dplyr::slice_sample(data.res, n=samples)

      #Re-estimate measure

      data.res <- dplyr::arrange(data.res,time, by_group=T)

      data.res <- dplyr::mutate(data.res, Happiness.lag = dplyr::lag(Happiness))
      data.res <- dplyr::mutate(data.res, Happiness.lead = dplyr::lead(Happiness))
      data.res$tkeo <- (data.res$Happiness*data.res$Happiness) - (data.res$Happiness.lag * data.res$Happiness.lead)
      data.res <- dplyr::mutate(data.res, tkeo.happ2 = sum(tkeo, na.rm = T)/dplyr::n())
      data.res <- dplyr::slice_head(data.res)
      data.res <- data.res[,c("id", "tkeo.happ2")]

      #Merge
      data.res <- merge(data.res, da, by="id")

      #Check if significant
      a <- stats::cor.test(data.res$tkeo.happ2, data.res$sim)
      res[i]<-ifelse(a$estimate>0 & a$p.value<p.value,1,0)

    }}



  #PAC
  ###############################################################################

  if(metric=="PAC"){

    #Estimate Measure all datapoints
    da <- dplyr::group_by(data, id)
    da <- dplyr::arrange(da,time, by_group=T)
    da <- dplyr::mutate(da, diff = Happiness - dplyr::lag(Happiness))
    da$abs.diff <- abs(da$diff)
    cut<-as.numeric(quantile(da$abs.diff, probs=0.9, na.rm = TRUE))
    da <- dplyr::mutate(da, pac.happ = sum(abs.diff>cut, na.rm = T)/(dplyr::n() - 1))
    da <- dplyr::slice_head(da)

    res <- c()

    for (i in 1:perm) {

      #Simulate variable
      da$sim <- faux::rnorm_pre(da$pac.happ, mu = 0, sd = 1, r = r)

      #Resample individuals
      rand.id <- sample(da$id, size=individuals)
      data.res <- data[data$id %in% rand.id,]

      #Resample observations per individual
      data.res <- dplyr::group_by(data.res, id)
      data.res <- dplyr::slice_sample(data.res, n=samples)

      #Re-estimate measure

      data.res <- dplyr::arrange(data.res,time, by_group=T)

      data.res <- dplyr::mutate(data.res, diff = Happiness - dplyr::lag(Happiness))


      data.res$abs.diff <- abs(data.res$diff)
      cut<-as.numeric(quantile(data.res$abs.diff, probs=0.9, na.rm = TRUE))
      data.res <- dplyr::mutate(data.res, pac.happ2 = sum(abs.diff>cut, na.rm = T)/(dplyr::n()-1))
      data.res <- dplyr::slice_head(data.res)
      data.res <- data.res[,c("id", "pac.happ2")]

      #Merge
      data.res <- merge(data.res, da, by="id")

      #Check if significant
      a <- stats::cor.test(data.res$pac.happ2, data.res$sim)
      res[i]<-ifelse(a$estimate>0 & a$p.value<p.value,1,0)

    }}

  #Autocorrelation
  ###############################################################################

  if(metric=="Autocorrelation"){

    #Estimate Measure all datapoints
    da <- dplyr::group_by(data, id)
    da <- dplyr::arrange(da,time, by_group=T)
    da <- dplyr::mutate(da, Happiness.t0 = dplyr::lag(Happiness))
    da <- dplyr::mutate(da, auto.happ = stats::cor(Happiness, Happiness.t0, use = "complete.obs"))
    da <- dplyr::slice_head(da)

    res <- c()

    for (i in 1:perm) {

      #Simulate variable
      da$sim <- faux::rnorm_pre(da$auto.happ, mu = 0, sd = 1, r = r)

      #Resample individuals
      rand.id <- sample(da$id, size=individuals)
      data.res <- data[data$id %in% rand.id,]

      #Resample observations per individual
      data.res <- dplyr::group_by(data.res, id)
      data.res <- dplyr::slice_sample(data.res, n=samples)

      #Re-estimate measure

      data.res <- dplyr::arrange(data.res,time, by_group=T)

      data.res <- dplyr::mutate(data.res, Happiness.t0 = dplyr::lag(Happiness))
      data.res <- dplyr::mutate(data.res, auto.happ2 = stats::cor(Happiness, Happiness.t0, use = "complete.obs"))
      data.res <- dplyr::slice_head(data.res)
      data.res <- data.res[,c("id", "auto.happ2")]

      #Merge
      data.res <- merge(data.res, da, by="id")

      #Check if significant
      a <- stats::cor.test(data.res$auto.happ2, data.res$sim)
      res[i]<-ifelse(a$estimate>0 & a$p.value<p.value,1,0)

    }}

  res.power<-mean(res, na.rm=T)

  if(metric=="Average"){metric.me<-"trait affect"}
  if(metric=="Rel.SD"){metric.me<-"the Relative Standard Deviation (Rel. SD) in affect"}
  if(metric=="SD"){metric.me<-"the Standard Deviation (SD) in affect"}
  if(metric=="RMSSD"){metric.me<-"the Root Mean Square of Successive Differences (RMSSD) in affect"}
  if(metric=="TKEO"){metric.me<-"the Teager-Kaiser Energy Operator (TKEO) of affect"}
  if(metric=="PAC"){metric.me<-"the Probability of Acute Change (PAC) in affect"}
  if(metric=="Autocorrelation"){metric.me<-"the Autocorrelation coefficient of affect"}


  cat("\nPower to detect a Pearson correlation of size r=", r,"between", metric.me,"and a given variable using a two-tailed t-test and an alpha of", p.value,".",
      "\n\nPower is estimated using the resampling and simulation approaches detailed in Pirla et al. (2021):\n")
  cat("\n -----------------------\nEmpirical power simulation under the sampling approach specified (",individuals, "individuals and", samples, "observations per individual):\n")
  cat("\n")
  cat("Aprox. Power=", res.power*100, "%")
  cat("\n -----------------------\nHow to report:\n")
  cat("\n")
  cat("Power analysis for affect dynamic studies (Pirla et al., 2021) suggests that our sampling strategy achieved a statistical power of", res.power*100, "% to detect a Pearson correlation of size r =",r,"using a two-tailed t-test with an alpha of",p.value,".")
  cat("\n -----------------------\nReference:\n")
  cat("\n")
  cat("Pirla, S., Taquet, M., & Quoidbach, J. (2021). Measuring Affect Dynamics: An Empirical Framework. https://doi.org/10.31219/osf.io/x2ywa")

  invisible(res.power)
}

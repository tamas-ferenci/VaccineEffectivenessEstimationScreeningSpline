Vaccine effectiveness (VE) estimation using the screening method,
upgraded with splines
================
Tamás Ferenci

## Overview

As already hinted in the original publication of Farrington in 1993 (C P
Farrington. Estimation of vaccine effectiveness using the screening
method. Int J Epidemiol. 1993 Aug;22(4):742-6. doi:
10.1093/ije/22.4.742,
<https://academic.oup.com/ije/article-abstract/22/4/742/664122>),
vaccine effectiveness (VE) estimation with the screening method can be
cast in the framework of logistic regression.

This comes as no surprise as VE (more precisely 1-VE) itself is an odds
ratio, since *VE = 1 - AR<sub>V</sub>/AR<sub>U</sub>*, where *AR* stands
for attack rate (number of infected divided by the size of the
population); *V* subscript indicates vaccinated, *U* denotes
unvaccinated subjects. By rearranging the terms, we arrive to the
following expression:

![1-=1-=1-=1-=1-=1-=1-](https://latex.codecogs.com/svg.image?VE%20=%201-%5Cfrac%7BAR_V%7D%7BAR_U%7D=1-%5Cfrac%7BI_V/N_V%7D%7BI_U/N_U%7D=1-%5Cfrac%7BI_V%7D%7BI_U%7D%5Ccdot%5Cfrac%7BN_U%7D%7BN_V%7D=1-%5Cfrac%7BI_V%7D%7BI-I_V%7D%5Ccdot%5Cfrac%7BN-N_V%7D%7BN_V%7D=1-%5Cfrac%7BI_V/I%7D%7B1-I_V/I%7D%5Ccdot%5Cfrac%7B1-N_V/N%7D%7BN_V/N%7D=1-%5Cfrac%7BPCV%7D%7B1-PCV%7D%5Ccdot%5Cfrac%7B1-PPV%7D%7BPPV%7D=1-%5Cfrac%7B%5Cfrac%7BPCV%7D%7B1-PCV%7D%7D%7B%5Cfrac%7BPPV%7D%7B1-PPV%7D%7D,)

where *I* and *N* denote the number of infected and total number of the
class indicated by the subscript, respectively, with *PCV* being the
proportion of cases vaccinated and *PPV* is the proportion of the
population vaccinated. Hence it can be directly seen that it corresponds
to a logistic regression where the number of vaccinated and unvaccinated
cases is the outcome, with a vaccination status indicator being the only
covariate. Logit of PPV should be entered as an offset (this assumes
that the PPV’s value is fixed and not a random variable).

This approach offers two advantages (as already pointed out by
Farrington in 1993). The first is that we simply obtain confidence
interval, the second is that further covariates can be added easily to
the regression, which represents a way to break down the analysis to
smaller strata, which allows some control of confounding (which is
otherwise the major problem of the screening method).

We make use of both of these advantages. The second is particularly
relevant here: we include calendar week, vaccine brand and age group as
covariates. Interaction was allowed between all three covariates, which
essentially means that separate time trends are allowed for every
vaccine brand and age group combination. While the failure to account
for comorbidities is one of the major limitations, the inclusion of age
somewhat alleviates this, as age is correlated with many chronic
diseases.

The inclusion of calendar week is of crucial importance, as it
implements the real-time monitoring. Here, we make a slight improvement
to the original screening method: we expand the calendar week with
splines, which essentially allows the evolution of VE over time to take
almost any – smooth – functional form, following the data. Splines are
rather flexible, but still use a limited number of degrees of freedom,
thus representing a compromise between parametric and fully
non-parametric models. This is, to our best knowledge, a novel solution
in this application.

Using splines in a logistic regression model essentially means that we
apply Generalized Additive Models (GAM).

Here we present the code that implements the above approach and a
simulation validation of it.

## Prerequisites

We will use R version 4.1.2 (2021-11-01) with packages `mgcv` (version
1.8.38), `data.table` (version 1.14.2), `ggplot2` (version 3.3.5).

Seed was set to 1 in this document.

## Code

The above approach is realized by the following function (the format of
the necessary input data will be made clear by the following section):

``` r
VEestim <- function(RawData, PPVData) {
  res <- merge(PPVData[TargetGroup!="ALL", .(UnVaccPPV = Denominator-sum(CumDose1),
                                             CumDose, Vaccine), .(TargetGroup, Year,
                                                                  Week = Week + 2)],
               merge(merge(CJ(Year = unique(RawData$Year), Week = unique(RawData$Week),
                              Vaccine = unique(RawData$Vaccine),
                              TargetGroup = unique(RawData$TargetGroup)),
                           RawData[, .(Vacc = sum(After)), .(Year, Week, Vaccine, TargetGroup)],
                           all.x = TRUE),
                     RawData[, .(UnVacc = sum(!AfterFirst)), .(Year, Week, TargetGroup)],
                     by = c("Year", "Week", "TargetGroup"))[,.(Vacc = sum(Vacc, na.rm = TRUE),
                                                               UnVacc = sum(UnVacc)),
                                                            .(Year, Week, Vaccine, TargetGroup)])
  
  res$PPV <- res$CumDose/(res$CumDose+res$UnVaccPPV)
  res$PCV <- res$Vacc/(res$Vacc+res$UnVacc)
  res$VE <- 1-(res$PCV/(1-res$PCV))*((1-res$PPV)/res$PPV)
  
  res <- res[PPV>0]
  
  fit <- mgcv::gam(cbind(Vacc, UnVacc) ~ s(Week, by = interaction(Vaccine, TargetGroup)) +
                     Vaccine * TargetGroup, offset = qlogis(PPV), data = res,
                   family = quasibinomial(), method = "REML", control = mgcv::gam.control(
                     nthreads = parallel::detectCores(logical = FALSE)-1))
  preds <- predict(fit, type = "link", se.fit = TRUE)
  res$VEsmoothed <- 1-exp(preds$fit)
  res$VEsmoothedUpr <- 1-exp(preds$fit-1.96*preds$se.fit)
  res$VEsmoothedLwr <- 1-exp(preds$fit+1.96*preds$se.fit)
  
  res
}
```

Thin plate regression splines are used, and estimation is done with
restricted maximum likelihood. The response distribution is
quasi-binomial to allow for potential overdispersion, which might be an
issue on real-life dataset (the case for allowing extra-Poisson
variability is already made by Farrington in 1993).

## Validation

### Creating a synthetic data set

We first simulate an epidemic (which will be artificial, but bears
resemblence to the usual epidemic curves with two ‘’waves’’). In more
detail, we take a baseline risk curve, and use multipliers to create the
epidemic curves for each age group. (I.e., there is no interaction, the
shape will be the same in each age group.)

``` r
RawData <- CJ(Week = 1:52,
              TargetGroup = c("Age<18", "Age18_24", "Age25_49", "Age50_59", "Age60_69",
                              "Age70_79", "Age80+"),
              Vaccine = c("COM", "MOD", "SPU", "AZ", "BECNBG", "JANSS"))
RawData <- merge(RawData, data.table(TargetGroup = c("Age<18", "Age18_24", "Age25_49", "Age50_59",
                                                     "Age60_69", "Age70_79", "Age80+"),
                                     Multiplier = c(0.45, 1, 1.1, 0.8, 0.4, 0.5, 0.6)),
                 by = "TargetGroup")
RawData <- merge(RawData, data.table(TargetGroup = c("Age<18", "Age18_24", "Age25_49", "Age50_59",
                                                     "Age60_69", "Age70_79", "Age80+"),
                                     Denominator = c(589350, 741094, 3490565, 1235206, 1295058,
                                                     859765, 438790)),
                 by = "TargetGroup")

RawData$UnvaccRisk <- (dnorm(RawData$Week, 12, 3)*0.05 +
                         dnorm(RawData$Week, 47, 3)*0.06 + 0.0001) * RawData$Multiplier

ggplot(RawData, aes(x = Week, y = UnvaccRisk, group = TargetGroup, color = TargetGroup)) +
  geom_line()
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

Next, we create VEs: each vaccine and each age group will have a pattern
over time (functional form). I.e., we use vaccines and age groups to
experiment with different ground truth values for VE. We will use ‘’age
group’’ to encode overall functional form (e.g., linear) and ‘’vaccine’’
to encode smaller variations within (such as different slopes of the
linear).

We try different functional forms to ensure robust validation: constant,
linear change, sinusoidal change, abrupt jump.

``` r
RawData$VID <- as.numeric(as.factor(RawData$Vaccine))
RawData$TID <- as.numeric(as.factor(RawData$TargetGroup))

RawData[, VE := switch(TID,
                       VID/10*Week/52,
                       VID/10*Week/52 + 0.3,
                       VID/10,
                       VID/10 - 0.3,
                       cos(Week/52*2*pi*VID/2)*0.5,
                       cos(Week/52*2*pi*2)*VID/10,
                       ifelse(Week<26, 0, VID/10)), .(TID)]
```

We download and reshape the vaccine uptake data (in contrast to the
previous, artifical data, these will be the actual values from the
ECDC):

``` r
if(!file.exists("PPVData2021.rds")) {
  PPVData <- fread("https://opendata.ecdc.europa.eu/covid19/vaccine_tracker/csv/data.csv")
  PPVData <- PPVData[Region=="HU"&TargetGroup%in%c("ALL", "Age<18", "Age18_24", "Age25_49",
                                                   "Age50_59", "Age60_69", "Age70_79", "Age80+")]
  PPVData$Year <- as.numeric(substring(PPVData$YearWeekISO, 1, 4))
  PPVData$Week <- as.numeric(substring(PPVData$YearWeekISO, 7, 8))
  PPVData[TargetGroup=="Age<18"]$Denominator <- 589350
  PPVData[TargetGroup=="ALL"]$Denominator <- 8649828
  PPVData <- PPVData[Year==2021]
  
  PPVData[, FirstDose := ifelse(TargetGroup!="ALL", FirstDose, sum(FirstDose[TargetGroup!="ALL"])),
          .(Vaccine, YearWeekISO)]
  PPVData[, SecondDose := ifelse(TargetGroup!="ALL", SecondDose,
                                 sum(SecondDose[TargetGroup!="ALL"])), .(Vaccine, YearWeekISO)]
  PPVData[, ThirdDose := ifelse(TargetGroup!="ALL", DoseAdditional1,
                                sum(DoseAdditional1[TargetGroup!="ALL"])), .(Vaccine, YearWeekISO)]
  PPVData[, CumDose1 := cumsum(FirstDose), .(Vaccine, TargetGroup)]
  PPVData[, CumDose2 := cumsum(SecondDose), .(Vaccine, TargetGroup)]
  PPVData[, CumDose3 := cumsum(ThirdDose), .(Vaccine, TargetGroup)]
  PPVData$CumDose <- ifelse(PPVData$Vaccine=="JANSS", PPVData$CumDose1, PPVData$CumDose2)
  saveRDS(PPVData, "PPVData2021.rds")
} else PPVData <- readRDS("PPVData2021.rds")
```

Finally, we simulate the number of people in the different categories
(unvaccinated, partially vaccinated, fully vaccinated) based on the
ground truth VE, and then convert the database to a case-based one
(noting that the vaccine’s effect’s onset is assumed to take two weeks):

``` r
GenSimData <- function(RawData, PPVData) {
  temp <- merge(RawData, PPVData[Year==2021&TargetGroup!="ALL",
                                 .(Week, TargetGroup, Vaccine, CumDose)],
                by = c("Week", "TargetGroup", "Vaccine"))
  temp$InfVacc <- rbinom(nrow(temp), temp$CumDose, temp$UnvaccRisk*(1-temp$VE))
  
  temp <- merge(temp, PPVData[Year==2021&TargetGroup!="ALL"&Vaccine!="JANSS",
                              .(Week, TargetGroup, Vaccine, CumDose1)],
                by = c("Week", "TargetGroup", "Vaccine"), all = TRUE)
  temp$InfPartvacc <- rbinom(nrow(temp), temp$CumDose-temp$CumDose1,
                             temp$UnvaccRisk*(1-temp$VE/2))
  temp[is.na(CumDose1)]$CumDose1 <- temp[is.na(CumDose1)]$CumDose
  temp[is.na(InfPartvacc)]$InfPartvacc <- 0
  
  temp2 <- rbind(
    temp[, .(rbinom(1, Denominator[1]-sum(CumDose1), UnvaccRisk[1])), .(Week, TargetGroup)][
      rep(seq_len(.N), V1), .(Week, TargetGroup, Vaccine = "", After = FALSE,
                              AfterFirst = FALSE)],
    temp[rep(seq_len(.N), InfPartvacc), .(Week, TargetGroup, Vaccine, After = FALSE,
                                          AfterFirst = TRUE)],
    temp[rep(seq_len(.N), InfVacc), .(Week, TargetGroup, Vaccine, After = TRUE,
                                      AfterFirst = TRUE)])
  temp2$Year <- 2021
  temp2$Week <- temp2$Week + 2
  temp2[Week<=52]
}
```

This concludes the creation of the synthetic data set.

### Single simulation

We estimate VE on the simulated data set (i.e., for which we know the
ground truth):

``` r
res <- GenSimData(RawData, PPVData)
res <- VEestim(res, PPVData)

res <- merge(res, RawData[, .(Week, TargetGroup, Vaccine, VEtrue = VE)],
             by = c("Week", "TargetGroup", "Vaccine"))
```

We create a plot to check, only visually for the time being, how well
the estimated VEs match the ground truth (dotted line shpws the
estimated, unsmoothed value, solid line indicates the estimated smoothed
value, dashed is the ground truth):

``` r
pargrid <- unique(res[, .(TargetGroup, Vaccine)])[order(TargetGroup)]
for(i in 1:nrow(pargrid))
  print(ggplot(res[TargetGroup==pargrid$TargetGroup[i]&Vaccine==pargrid$Vaccine[i]],
               aes(x = Week, ymin = VEsmoothedLwr, ymax = VEsmoothedUpr)) +
          geom_line(aes(y = VEsmoothed)) + geom_line(aes(y = VEtrue), linetype = "dashed") +
          geom_line(aes(y = VE), linetype = "dotted") +
          geom_ribbon(alpha = 0.2, aes(linetype = NA)) + coord_cartesian(ylim = c(-0.5, 1)) +
          labs(y = "Vaccine effectiveness [%]"))
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-3.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-4.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-5.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-6.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-7.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-8.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-9.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-10.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-11.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-12.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-13.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-14.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-15.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-16.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-17.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-18.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-19.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-20.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-21.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-22.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-23.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-24.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-25.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-26.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-27.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-28.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-29.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-30.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-31.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-32.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-33.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-34.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-35.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-36.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-37.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-7-38.png)<!-- -->

The method seems to correctly pick up every case, except the abrupt
change (unsurprisingly, as splines are meant to model smooth functions).
We also see some lag in picking up the sinusoidal signals, especially at
lower sample sizes.

Overall, the smoothing works well, but depends on the sample size (lower
between simulated waves).

### Replicated simulations

The above simulation however fails to investigate the impact of random
fluctuations, i.e., unbiasedness. (In terms of case numbers – vaccine
coverage is assumed to be fixed.) This can be resolved by re-running the
simulations several times (thin lines indicate individual values, bold
line is their pointwise median; unsmoothed estimates are no longer
shown):

``` r
res <- rbindlist(lapply(1:10, function(i)
  merge(VEestim(GenSimData(RawData, PPVData), PPVData),
        RawData[, .(Week, TargetGroup, Vaccine, VEtrue = VE)],
        by = c("Week", "TargetGroup", "Vaccine"))), idcol = TRUE)

for(i in 1:nrow(pargrid))
  print(ggplot(res[TargetGroup==pargrid$TargetGroup[i]&Vaccine==pargrid$Vaccine[i]],
               aes(x = Week, ymin = VEsmoothedLwr, ymax = VEsmoothedUpr)) +
          geom_line(aes(y = VEsmoothed, group = .id), alpha = 0.2) +
          geom_line(aes(y = VEtrue), linetype = "dashed") + coord_cartesian(ylim = c(-0.5, 1)) +
          labs(y = "Vaccine effectiveness [%]"))
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-4.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-5.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-6.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-7.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-8.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-9.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-10.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-11.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-12.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-13.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-14.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-15.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-16.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-17.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-18.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-19.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-20.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-21.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-22.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-23.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-24.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-25.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-26.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-27.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-28.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-29.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-30.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-31.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-32.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-33.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-34.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-35.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-36.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-37.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-8-38.png)<!-- -->

The overall picture is the same.

## Possible future research directions

-   Investigating the impact of the wrong specification for the lead
    time in PPV calculation.
-   Changing vaccination uptake (i.e., ECDC data) as well.
-   Verifying that the coverage of the confidence intervals is indeed
    95%.
-   Experimenting with other age structures (e. g., interaction between
    age and the epidemic curve).
-   Validating the GAM (e.g., *k* value, or different basis function; or
    `gam.check()`)
-   Trying adaptive basis (i.e., `bs = "ad"`); seems to have huge
    computational burden.

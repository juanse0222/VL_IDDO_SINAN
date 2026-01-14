library(mgcv)
library(survival) ## for data
col1 <- colon[colon$etype==1,] ## concentrate on single event
col1$differ <- as.factor(col1$differ)
col1$sex <- as.factor(col1$sex)

b <- gam(time~s(age,by=sex)+sex+s(nodes)+perfor+rx+obstruct+adhere,
         family=cox.ph(),data=col1,weights=status)

summary(b) 

plot(b,pages=1,all.terms=TRUE) ## plot effects

plot(b$linear.predictors,residuals(b))

## plot survival function for patient j...

np <- 300;j <- 6
newd <- data.frame(time=seq(0,3000,length=np))
dname <- names(col1)
for (n in dname) newd[[n]] <- rep(col1[[n]][j],np)
newd$time <- seq(0,3000,length=np)
fv <- predict(b,newdata=newd,type="response",se=TRUE)
plot(newd$time,fv$fit,type="l",ylim=c(0,1),xlab="time",ylab="survival")
lines(newd$time,fv$fit+2*fv$se.fit,col=2)
lines(newd$time,fv$fit-2*fv$se.fit,col=2)

## crude plot of baseline survival...

plot(b$family$data$tr,exp(-b$family$data$h),type="l",ylim=c(0,1),
     xlab="time",ylab="survival")
lines(b$family$data$tr,exp(-b$family$data$h + 2*b$family$data$q^.5),col=2)
lines(b$family$data$tr,exp(-b$family$data$h - 2*b$family$data$q^.5),col=2)
lines(b$family$data$tr,exp(-b$family$data$km),lty=2) ## Kaplan Meier

## Checking the proportional hazards assumption via scaled score plots as
## in Klein and Moeschberger Section 11.4 p374-376... 

ph.resid <- function(b,stratum=1) {
  ## convenience function to plot scaled score residuals against time,
  ## by term. Reference lines at 5% exceedance prob for Brownian bridge
  ## (see KS test statistic distribution).
  rs <- residuals(b,"score");term <- attr(rs,"term")
  if (is.matrix(b$y)) {
    ii <- b$y[,2] == stratum;b$y <- b$y[ii,1];rs <- rs[ii,]
  }
  oy <- order(b$y)
  for (i in 1:length(term)) {
    ii <- term[[i]]; m <- length(ii)
    plot(b$y[oy],rs[oy,ii[1]],ylim=c(-3,3),type="l",ylab="score residuals",
         xlab="time",main=names(term)[i])
    if (m>1) for (k in 2:m) lines(b$y[oy],rs[oy,ii[k]],col=k);
    abline(-1.3581,0,lty=2);abline(1.3581,0,lty=2) 
  }  
}
par(mfrow=c(2,2))
ph.resid(b)

## stratification example, with 2 randomly allocated strata
## so that results should be similar to previous....
col1$strata <- sample(1:2,nrow(col1),replace=TRUE) 
bs <- gam(cbind(time,strata)~s(age,by=sex)+sex+s(nodes)+perfor+rx+obstruct
          +adhere,family=cox.ph(),data=col1,weights=status)
plot(bs,pages=1,all.terms=TRUE) ## plot effects

## baseline survival plots by strata...

for (i in 1:2) { ## loop over strata
  ## create index selecting elements of stored hazard info for stratum i...
  ind <- which(bs$family$data$tr.strat == i)
  if (i==1) plot(bs$family$data$tr[ind],exp(-bs$family$data$h[ind]),type="l",
                 ylim=c(0,1),xlab="time",ylab="survival",lwd=2,col=i) else
                   lines(bs$family$data$tr[ind],exp(-bs$family$data$h[ind]),lwd=2,col=i)
  lines(bs$family$data$tr[ind],exp(-bs$family$data$h[ind] +
                                     2*bs$family$data$q[ind]^.5),lty=2,col=i) ## upper ci
  lines(bs$family$data$tr[ind],exp(-bs$family$data$h[ind] -
                                     2*bs$family$data$q[ind]^.5),lty=2,col=i) ## lower ci
  lines(bs$family$data$tr[ind],exp(-bs$family$data$km[ind]),col=i) ## KM
}


## Simple simulated known truth example...
ph.weibull.sim <- function(eta,gamma=1,h0=.01,t1=100) { 
  lambda <- h0*exp(eta) 
  n <- length(eta)
  U <- runif(n)
  t <- (-log(U)/lambda)^(1/gamma)
  d <- as.numeric(t <= t1)
  t[!d] <- t1
  list(t=t,d=d)
}
n <- 500;set.seed(2)
x0 <- runif(n, 0, 1);x1 <- runif(n, 0, 1)
x2 <- runif(n, 0, 1);x3 <- runif(n, 0, 1)
f0 <- function(x) 2 * sin(pi * x)
f1 <- function(x) exp(2 * x)
f2 <- function(x) 0.2*x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10
f3 <- function(x) 0*x
f <- f0(x0) + f1(x1)  + f2(x2)
g <- (f-mean(f))/5
surv <- ph.weibull.sim(g)
surv$x0 <- x0;surv$x1 <- x1;surv$x2 <- x2;surv$x3 <- x3

b <- gam(t~s(x0)+s(x1)+s(x2,k=15)+s(x3),family=cox.ph,weights=d,data=surv)

plot(b,pages=1)

## Another one, including a violation of proportional hazards for
## effect of x2...

set.seed(2)
h <- exp((f0(x0)+f1(x1)+f2(x2)-10)/5)
t <- rexp(n,h);d <- as.numeric(t<20)

## first with no violation of PH in the simulation...
b <- gam(t~s(x0)+s(x1)+s(x2)+s(x3),family=cox.ph,weights=d)
plot(b,pages=1)
ph.resid(b) ## fine

## Now violate PH for x2 in the simulation...
ii <- t>1.5
h1 <- exp((f0(x0)+f1(x1)+3*f2(x2)-10)/5)
t[ii] <- 1.5 + rexp(sum(ii),h1[ii]);d <- as.numeric(t<20)

b <- gam(t~s(x0)+s(x1)+s(x2)+s(x3),family=cox.ph,weights=d)
plot(b,pages=1)
ph.resid(b) ## The checking plot picks up the problem in s(x2) 


## conditional logistic regression models are often estimated using the 
## cox proportional hazards partial likelihood with a strata for each
## case-control group. A dummy vector of times is created (all equal). 
## The following compares to 'clogit' for a simple case. Note that
## the gam log likelihood is not exact if there is more than one case
## per stratum, corresponding to clogit's approximate method.
library(survival);library(mgcv)
infert$dumt <- rep(1,nrow(infert))
mg <- gam(cbind(dumt,stratum) ~ spontaneous + induced, data=infert,
          family=cox.ph,weights=case)
ms <- clogit(case ~ spontaneous + induced + strata(stratum), data=infert,
             method="approximate")
summary(mg)$p.table[1:2,]; ms
################################################################################################################


# (re-fit) fit the model again so it's linked to the datadist metadata
fit_frailty <- cph(Surv(dif_treat_outcome, death_vl) ~ rcs(dif_day_tx_symp, 4) + rcs(age_years, 4) + 
                     hiv_3,
                   data = data_8)

# now plotting should work
plot(Predict(fit_frailty, dif_day_tx_symp, fun = exp), xlab = "Delay days", ylab = "Hazard ratio")

summary(fit_frailty)



fit_coxme <- coxme(Surv(dif_treat_outcome, death_vl) ~ rcs(dif_day_tx_symp, 4) + rcs(age_years,4) + hiv_3 + 
                     (1|federative_unit),
                   data = data_8)
summary(fit_coxme)

# Option A: marginal model with cluster robust SEs
fit_cox_cluster <- coxph(Surv(dif_treat_outcome, death_vl) ~ rcs(dif_day_tx_symp,4) + rcs(age_years,4) + hiv_3,
                         data = data_8, cluster = federative_unit)

summary(fit_cox_cluster)

ggforest(fit_cox_cluster, data = data_8)

# Build newdata grid for prediction (100 points)
delay_seq <- seq(min(data_8$dif_day_tx_symp, na.rm=TRUE),
                 max(data_8$dif_day_tx_symp, na.rm=TRUE),
                 length.out = 100)
ref_age <- median(data_8$age_years, na.rm = TRUE)
ref_hiv <- if (is.factor(data_8$hiv_3)) levels(data_8$hiv_3)[1] else mode_val(data_8$hiv_3)

newd <- data.frame(dif_day_tx_symp = delay_seq,
                   age_years = rep(ref_age, length(delay_seq)),
                   hiv_3 = if (is.factor(data_8$hiv_3)) 
                     factor(rep(ref_hiv, length(delay_seq)), levels = levels(data_8$hiv_3)) 
                   else rep(ref_hiv, length(delay_seq)))

# Predict linear predictor and se
lp <- predict(fit_cox_cluster, newdata = newd, type = "lp", se.fit = TRUE)

newd$HR   <- exp(lp$fit)
newd$low  <- exp(lp$fit - 1.96*lp$se.fit)
newd$high <- exp(lp$fit + 1.96*lp$se.fit)

# Optional: make HR relative to delay = e.g. 30 days
ref_day <- 30
idx_ref <- which.min(abs(delay_seq - ref_day))
newd$HR_rel  <- newd$HR / newd$HR[idx_ref]
newd$low_rel <- newd$low / newd$HR[idx_ref]
newd$high_rel<- newd$high / newd$HR[idx_ref]

# Plot (relative HR)
plot(delay_seq, newd$HR_rel, type = "l", lwd = 2, xlab = "Delay days",
     ylab = paste0("HR (ref = ", delay_seq[idx_ref], " days)"),
     ylim = c(min(newd$low_rel, na.rm=TRUE), max(newd$high_rel, na.rm=TRUE)))
lines(delay_seq, newd$low_rel, lty = 2)
lines(delay_seq, newd$high_rel, lty = 2)
abline(h = 1, col = "grey50")



# fixed effects (log-hazard coefficients), se, HR and 95% CI
fixef <- coef(fit_coxme)                    # log(HR) for fixed effects
vcov_fix <- vcov(fit_coxme)                 # variance-covariance matrix of fixed effects
se_fix  <- sqrt(diag(vcov_fix))
HR      <- exp(fixef)
ci_low  <- exp(fixef - 1.96 * se_fix)
ci_high <- exp(fixef + 1.96 * se_fix)

fixed_effects_table <- data.frame(
  term = names(fixef),
  coef = as.numeric(fixef),
  se = se_fix,
  HR = HR,
  CI95_low = ci_low,
  CI95_high = ci_high,
  row.names = NULL
)
print(fixed_effects_table)


data_8.1 <- data_8[data_8$death_vl == 1,]

data_8.1 <- data_8.1 %>% 
  mutate(dif_treat_start_symp_cat = case_when(dif_treat_start_symp_cat == "< 30 days" ~ 1,
                                              dif_treat_start_symp_cat == "30-89 days" ~ 2,
                                              dif_treat_start_symp_cat == "90-119 days" ~ 3,
                                              dif_treat_start_symp_cat == "≥120 days" ~ 4,
                                              TRUE ~ NA))


b <- gam(dif_treat_outcome ~ s(dif_day_tx_symp) + s(age_years,by=gender_2) + hiv_3,
         family=cox.ph(), data=data_8.1, weights= dif_treat_start_symp_cat)

summary(b) 

plot(b,pages=1,all.terms=TRUE) ## plot effects

plot(b$linear.predictors,residuals(b))




library(survival)
library(splines)

fit_cox_ns <- coxph(Surv(dif_treat_outcome, death_vl) ~ ns(dif_day_tx_symp, df=3) + ns(age_years, df=3) + gender_2 + hiv_3,
                    data = data_8)

# create a grid of dif_day_tx_symp to predict HR (set other covariates to typical values)
newd <- data.frame(dif_day_tx_symp = seq(min(data_8$dif_day_tx_symp, na.rm=TRUE),
                                         max(data_8$dif_day_tx_symp, na.rm=TRUE), length=100),
                   age_years = median(data_8$age_years, na.rm=TRUE),
                   gender_2 = levels(data_8$gender_2)[1],
                   hiv_3 = levels(data_8$hiv_3)[1])

# get linear predictor and CI
lp <- predict(fit_cox_ns, newdata = newd, se.fit = TRUE, type = "lp")
newd$HR  <- exp(lp$fit)
newd$low <- exp(lp$fit - 1.96*lp$se.fit)
newd$high<- exp(lp$fit + 1.96*lp$se.fit)

# plot with base or ggplot
plot(newd$dif_day_tx_symp, newd$HR, type="l", ylim=c(min(newd$low), max(newd$high)),
     xlab="Delay days", ylab="Hazard ratio")
lines(newd$dif_day_tx_symp, newd$low, lty=2)
lines(newd$dif_day_tx_symp, newd$high, lty=2)





# 1) create sequence for delay (100 points)
delay_seq <- seq(min(data_8$dif_day_tx_symp, na.rm = TRUE),
                 max(data_8$dif_day_tx_symp, na.rm = TRUE),
                 length.out = 100)

# 2) choose reference values for other covariates
# median age
ref_age <- median(data_8$age_years, na.rm = TRUE)

# for categorical vars use the most frequent (mode). This handles numeric or factor columns.
mode_val <- function(x) {
  ux <- unique(na.omit(x))
  ux[which.max(tabulate(match(x, ux)))]
}
ref_gender <- mode_val(data_8$gender_2)
ref_hiv    <- mode_val(data_8$hiv_3)

# 3) ensure reference values use same factor levels as in the original data
# If the original columns are factors, create factors with the same levels
if (is.factor(data_8$gender_2)) {
  new_gender <- factor(rep(ref_gender, length(delay_seq)), levels = levels(data_8$gender_2))
} else {
  new_gender <- rep(ref_gender, length(delay_seq))
}

if (is.factor(data_8$hiv_3)) {
  new_hiv <- factor(rep(ref_hiv, length(delay_seq)), levels = levels(data_8$hiv_3))
} else {
  new_hiv <- rep(ref_hiv, length(delay_seq))
}

# 4) build newd with exactly 100 rows
newd <- data.frame(
  dif_day_tx_symp = delay_seq,
  age_years = rep(ref_age, length(delay_seq)),
  gender_2 = new_gender,
  hiv_3 = new_hiv
)

# 5) predict linear predictor and standard error from the fitted coxph
lp <- predict(fit_cox_ns, newdata = newd, se.fit = TRUE, type = "lp")

# 6) compute HR and 95% CI (using reference at the first delay point or relative to median)
# Here we present HR relative to baseline equal to reference covariates at dif_day_tx_symp = delay_seq[1].
# To show HR relative to a particular reference day (e.g., 30 days), adjust below.
newd$HR  <- exp(lp$fit)
newd$low <- exp(lp$fit - 1.96 * lp$se.fit)
newd$high<- exp(lp$fit + 1.96 * lp$se.fit)

# 7) (Optional) if you want HR relative to a chosen reference day, say ref_day <- 30:
# ref_day <- 30
# idx_ref <- which.min(abs(delay_seq - ref_day))
# newd$HR_rel  <- newd$HR / newd$HR[idx_ref]
# newd$low_rel <- newd$low / newd$HR[idx_ref]
# newd$high_rel<- newd$high / newd$HR[idx_ref]

# 8) Plot HR (on multiplicative scale)
plot(newd$dif_day_tx_symp, newd$HR, type = "l", lwd = 2,
     xlab = "Delay days", ylab = "Hazard ratio (relative scale)",
     ylim = c(min(newd$low, na.rm = TRUE), max(newd$high, na.rm = TRUE)))
lines(newd$dif_day_tx_symp, newd$low, lty = 2)
lines(newd$dif_day_tx_symp, newd$high, lty = 2)
abline(h = 1, col = "grey60")

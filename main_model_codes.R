####################################################################################################
# Example of code written by Jia Wei (Nuffield Department of Medicine, University of Oxford)       #
# for analyses of antibody waning following the first vaccination.                                 #
# Accompanying paper: SARS-CoV-2 antibody trajectories after a single COVID-19 vaccination with    # 
# and without prior infection.                                                                     #
####################################################################################################

library(tidyverse)
library(brms)

priors=c(
  set_prior("normal(8,2)",class="Intercept"),
  set_prior("normal(0,0.1)",coef="time",class="b"),
  set_prior("normal(0,1)",coef="age",class="b"),
  set_prior("normal(0,1)",coef="Male1",class="b"),
  set_prior("normal(0,1)",coef="ethnicity1",class="b"),
  set_prior("normal(0,1)",coef="lthc1",class="b"),
  set_prior("normal(0,1)",coef="hcw1",class="b"),
  set_prior("normal(0,1)",coef="IMD",class="b"),
  set_prior("normal(0,2)",coef="prior1",class="b"),
  
  set_prior("normal(0,0.1)",coef="time:age",class="b"),
  set_prior("normal(0,0.1)",coef="time:Male1",class="b"),
  set_prior("normal(0,0.1)",coef="time:ethnicity1",class="b"),
  set_prior("normal(0,0.1)",coef="time:lthc1",class="b"),
  set_prior("normal(0,0.1)",coef="time:hcw1",class="b"),
  set_prior("normal(0,0.1)",coef="time:IMD",class="b"),
  set_prior("normal(0,0.1)",coef="time:prior1",class="b"),
  
  set_prior("cauchy(0,0.01)",coef="time",class="sd",group="ID"),
  set_prior("cauchy(0,1)",coef="Intercept",class="sd",group="ID"),
  set_prior("cauchy(0,0.5)",class="sigma"),
  set_prior("lkj_corr_cholesky(1)",class="L")
)

set_inits=function(seed=1){
  set.seed(seed)
  list(
    Intercept=rnorm(1,8,1),
    b=c(rnorm(1,-0.01,0.01),#time
        rnorm(1,0,0.1),#age
        rnorm(1,0,0.1),#sex
        rnorm(1,0,0.1),#ethnicity
        rnorm(1,0,0.1),#lthc
        rnorm(1,0,0.1),#hcw
        rnorm(1,0,0.1),#IMD
        rnorm(1,0,0.1),#prior
        
        rnorm(1,0,0.01),#age
        rnorm(1,0,0.01),#sex
        rnorm(1,0,0.01),#ethnicity
        rnorm(1,0,0.01),#lthc
        rnorm(1,0,0.01),#hcw
        rnorm(1,0,0.01),#IMD
        rnorm(1,0,0.01)#prior
    ),
    sigma=runif(1,0,1),
    sd_1=c(runif(1,0,1),runif(1,0,0.02)),
    z_1=matrix(rep(c(7.5,-0.01),35638),2,35638)
  )
}

inits_list=list(
  set_inits(1),
  set_inits(2),
  set_inits(3),
  set_inits(4)
)

fit <- brm(formula=log2(assay)|cens(cen1)~
                1+time*age+time*Male+time*ethnicity+time*lthc+time*hcw+time*IMD+
                time*prior+(1+time|ID),
              data=data,cores=4,family=gaussian(),
              prior = priors,
              chains = 4, iter=4000, warmup = 2000, seed=42,inits=inits_list,
              control = list(adapt_delta=0.95))

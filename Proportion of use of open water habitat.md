# Proportion of time spent in open water of lakes

## 1. Load data and transform variables

:books:`library(tidyverse)`

```
open_water_prop <- data.table(read_csv("open_water_prop.csv"))
```
```
open_water_prop$up_lake <- as.factor(open_water_prop$up_lake)
open_water_prop$tag_id <- as.factor(open_water_prop$tag_id)
open_water_prop<-horiz_move[with(open_water_prop, order(tag_id,up_lake, date)),]
```

## 2. Multicolinearity tests

:books:`library(corrplot)`

```
corr<-cor(open_water_prop[,c(3,5,6,7,11,12)], use="pairwise", method="spearman")
corrplot(corr,method="number")
```
![Hab_prop](/Plots/Hab_prop_zoib0.png "Hab_prop")

:books:`library(ppcor)`

```
pcor(open_water_prop[,c(3,5,6,7,11,12)],method = "pearson")$estimate
```
```
               hab.prop       area95    temp.day  day_length       tl_mm       k1size
hab.prop    1.000000000 -0.007701849 -0.01214729  0.00564384 -0.02615858  0.005441638
area95     -0.007701849  1.000000000 -0.05276799  0.02259822  0.44379544  0.155270414
temp.day   -0.012147286 -0.052767993  1.00000000  0.16407162  0.02087628 -0.065199009
day_length  0.005643840  0.022598221  0.16407162  1.00000000  0.08018597 -0.088232886
tl_mm      -0.026158585  0.443795442  0.02087628  0.08018597  1.00000000 -0.061440686
k1size      0.005441638  0.155270414 -0.06519901 -0.08823289 -0.06144069  1.000000000
```

## 3. Fit Autoregressive Moving Average ARMA(p, q) model to data of proportion in open water use

### 3.1. Fit a linear mixed-effects model using package _nlme_

:books:`library(nlme)`

#### 3.1.1. Full model

```
model.ow <- lme(asin(sqrt(hab.prop)) ~ temp.day*up_lake +day_length*up_lake+up_lake*tl_mm+temp.day*day_length, data = open_water_prop, random = ~1|tag_id, method = "REML", correlation = corARMA(value = c(0.4276,  -0.9426),~date.num|tag_id, p = 1, q=1))
```

Lets compare estimates of regression coefficients with and without scaling predictors

#### 3.1.2. Final model (unscaled predictors)

```
model.ow.f <- lme(asin(sqrt(hab.prop)) ~  day_length*up_lake+tl_mm, data = open_water_prop ,random = ~1|tag_id, method = "REML", correlation = corARMA(value = c(0.4276,  -0.9426),~date.num|tag_id, p = 1, q=1))
```
```
summary(model.ow.f)
```
```
Linear mixed-effects model fit by REML
 Data: open_water_prop
        AIC       BIC   logLik
  -698.1296 -647.9909 358.0648

Random effects:
 Formula: ~1 | tag_id
        (Intercept)  Residual
StdDev:    0.191679 0.2437836

Correlation Structure: ARMA(1,1)
 Formula: ~date.num | tag_id
 Parameter estimate(s):
      Phi1     Theta1
 0.8167204 -0.3525310
Fixed effects: asin(sqrt(hab.prop)) ~ day_length * up_lake + tl_mm
                            Value  Std.Error   DF    t-value p-value
(Intercept)            -0.1881128 0.29847195 1919 -0.6302528  0.5286
day_length             -0.0182465 0.01120119 1919 -1.6289818  0.1035
up_lakeMost            -0.3786813 0.24774028   22 -1.5285413  0.1406
tl_mm                   0.0007454 0.00032128   22  2.3201986  0.0300
day_length:up_lakeMost  0.0347729 0.01566051 1919  2.2204160  0.0265
 Correlation:
                       (Intr) dy_lng up_lkM tl_mm
day_length             -0.534
up_lakeMost            -0.283  0.669
tl_mm                  -0.810 -0.023 -0.161
day_length:up_lakeMost  0.382 -0.715 -0.934  0.016

Standardized Within-Group Residuals:
       Min         Q1        Med         Q3        Max
-2.9547532 -0.3840864 -0.1747805  0.2897519  6.1412196

Number of Observations: 1946
Number of Groups: 25
```
```
vif(model.ow.f)
```
```
        day_length            up_lake              tl_mm day_length:up_lake
          2.047923           9.431874           1.201934          10.251848
```

- There is collinearity from _up_lake_ so we should scale predictors

#### 3.1.3. Final model (scaled predictors)

:books:`library(dplyr)`

```
open_water_prop_scale_p<-open_water_prop %>% mutate_at(c(7,11), funs(c(scale(.))))
```
```
model.ow.f_scale_p <- lme(asin(sqrt(hab.prop)) ~  day_length*up_lake+tl_mm, data = open_water_prop_scale_p ,random = ~1|tag_id, method = "REML", correlation = corARMA(value = c(0.4276,  -0.9426),~date.num|tag_id, p = 1, q=1))
```
```
summary(model.ow.f_scale_p)
```
```
Linear mixed-effects model fit by REML
 Data: open_water_prop_scale_p
        AIC      BIC   logLik
  -709.6716 -659.533 363.8358

Random effects:
 Formula: ~1 | tag_id
        (Intercept)  Residual
StdDev:    0.191679 0.2437836

Correlation Structure: ARMA(1,1)
 Formula: ~date.num | tag_id
 Parameter estimate(s):
      Phi1     Theta1
 0.8167204 -0.3525310
Fixed effects: asin(sqrt(hab.prop)) ~ day_length * up_lake + tl_mm
                             Value  Std.Error   DF   t-value p-value
(Intercept)             0.16635624 0.06217504 1919  2.675611  0.0075
day_length             -0.02764855 0.01697290 1919 -1.628982  0.1035
up_lakeMost             0.12896493 0.08845954   22  1.457897  0.1590
tl_mm                   0.10417117 0.04489753   22  2.320199  0.0300
day_length:up_lakeMost  0.05269050 0.02373001 1919  2.220416  0.0265
 Correlation:
                       (Intr) dy_lng up_lkM tl_mm
day_length             -0.032
up_lakeMost            -0.760  0.026
tl_mm                   0.353 -0.023 -0.410
day_length:up_lakeMost  0.022 -0.715 -0.032  0.016

Standardized Within-Group Residuals:
       Min         Q1        Med         Q3        Max
-2.9547532 -0.3840864 -0.1747805  0.2897519  6.1412196

Number of Observations: 1946
Number of Groups: 25
```
```
vif(model.ow.f_scale_p)
```
```
        day_length            up_lake              tl_mm day_length:up_lake
          2.047923           1.202525           1.201934           2.048239
```

### 3.2. Fit a beta regression model (distribution with logit/beta processes) using the package betareg

:books:`library(betareg)`

First, we need to check if there are NAs for the variable _hab.prop_. We can do this with the following code:
```
open_water_prop[is.na(open_water_prop$hab.prop),]
```
- There is not any NA so we can proceed with the following step
- The variable must be bound between 0 and 1 but excluding the extremes as the _betareg()_ function does not allow these values. We use two approaches:

**Note:** If we dont perform this transformation we get an error due to 0's and 1's contained in the extremes since the betareg function cant handle those values

#### 3.2.1. Transform the response variable with values higher than 0 and lower than 1 but in close decimal place

```
open_water_prop <- open_water_prop %>% mutate(hab.prop = replace(hab.prop, hab.prop == 1, 0.999))
open_water_prop <- open_water_prop %>% mutate(hab.prop = replace(hab.prop, hab.prop == 0, 0.001))

open_water_prop_scale_p <- open_water_prop_scale_p %>% mutate(hab.prop = replace(hab.prop, hab.prop == 1, 0.999))
open_water_prop_scale_p <- open_water_prop_scale_p %>% mutate(hab.prop = replace(hab.prop, hab.prop == 0, 0.001))
```

**Final model (unscaled predictors)**

```
model.ow.f_beta <- betareg(hab.prop ~  day_length*up_lake+tl_mm, data = open_water_prop ,random = ~1|tag_id, method = "Nelder-Mead", correlation = corARMA(value = c(0.4276,  -0.9426),~date.num|tag_id, p = 1, q=1))
```
```
summary(model.ow.f_beta)
```
```
Call:
betareg(formula = hab.prop ~ day_length * up_lake + tl_mm, data = open_water_prop, random = ~1 | tag_id, method = "Nelder-Mead",
    correlation = corARMA(value = c(0.4276, -0.9426), ~date.num | tag_id, p = 1, q = 1))

Standardized weighted residuals 2:
    Min      1Q  Median      3Q     Max
-1.5319 -0.5091 -0.0629  0.5215  3.0120

Coefficients (mean model with logit link):
                         Estimate Std. Error z value Pr(>|z|)
(Intercept)            -2.7825714  0.3421519  -8.133  4.2e-16 ***
day_length             -0.0490611  0.0223471  -2.195   0.0281 *
up_lakeMost            -0.7548783  0.4861765  -1.553   0.1205
tl_mm                   0.0018215  0.0001892   9.625  < 2e-16 ***
day_length:up_lakeMost  0.0744568  0.0329860   2.257   0.0240 *

Phi coefficients (precision model with identity link):
      Estimate Std. Error z value Pr(>|z|)
(phi)  1.78439    0.06611   26.99   <2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Type of estimator: ML (maximum likelihood)
Log-likelihood:  4311 on 6 Df
Pseudo R-squared: 0.192
Number of iterations: 4483 (Nelder-Mead) + 5 (Fisher scoring)
```

**Final model (scaled predictors)**

```
model.ow.f_beta_scaled <- betareg(hab.prop ~  day_length*up_lake+tl_mm, data = open_water_prop_scale_p ,random = ~1|tag_id, method = "Nelder-Mead", correlation = corARMA(value = c(0.4276,  -0.9426),~date.num|tag_id, p = 1, q=1))
```
```
summary(model.ow.f_beta_scaled)
```
```
Call:
betareg(formula = hab.prop ~ day_length * up_lake + tl_mm, data = open_water_prop_scale_p, random = ~1 | tag_id, method = "Nelder-Mead",
    correlation = corARMA(value = c(0.4276, -0.9426), ~date.num | tag_id, p = 1, q = 1))

Standardized weighted residuals 2:
    Min      1Q  Median      3Q     Max
-1.5319 -0.5091 -0.0629  0.5215  3.0120

Coefficients (mean model with logit link):
                       Estimate Std. Error z value Pr(>|z|)
(Intercept)            -1.98176    0.04547 -43.584  < 2e-16 ***
day_length             -0.07434    0.03386  -2.195   0.0281 *
up_lakeMost             0.33211    0.05252   6.324 2.55e-10 ***
tl_mm                   0.25454    0.02645   9.625  < 2e-16 ***
day_length:up_lakeMost  0.11282    0.04998   2.257   0.0240 *

Phi coefficients (precision model with identity link):
      Estimate Std. Error z value Pr(>|z|)
(phi)  1.78439    0.06611   26.99   <2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Type of estimator: ML (maximum likelihood)
Log-likelihood:  4311 on 6 Df
Pseudo R-squared: 0.192
Number of iterations: 2363 (Nelder-Mead) + 5 (Fisher scoring)
```

#### 3.2.2. Transform the response variable following the recommendations by Smithson & Verkuilen (2006)

The R _betareg_ documentation mentions the distribution proposed by these authors, which takes the form:

**y’ = (y*(N - 1) + .5)/N**, where N is the sample size

Calculate sample size

```
nrow(open_water_prop)/(length(levels(open_water_prop$tag_id)))
```
```
[1] 77.84
```

open_water_prop$hab.prop <- (open_water_prop$hab.prop*(77.84-1) + 0.5)/77.84

open_water_prop$hab.prop_t <- with(open_water_prop, hab.prop*(77.84-1) + 0.5)/77.84)

hab.prop_df <- data.frame(open_water_prop$hab.prop)
open_water_prop$hab.prop_t <- apply(hab.prop_df, 1, function(x) prod((hab.prop_df$open_water_prop.hab.prop*(77.84-1) + 0.5)/77.84))


### 3.3. Fit a zero- /one- inflated beta regression model

Now we fit a series of zero-inflated models using different packages on the original dataset without further transformation of the (0,1) bound

```
open_water_prop <- data.table(read_csv("open_water_prop.csv"))
```
```
open_water_prop$up_lake <- as.factor(open_water_prop$up_lake)
open_water_prop$tag_id <- as.factor(open_water_prop$tag_id)
```

#### 3.3.1. Using package _gamlss_

:books:`library(gamlss)`

```
model.ow.f_beta_zi_gamlss <- gamlss(hab.prop ~  day_length*up_lake+tl_mm, data = open_water_prop ,random = ~1|tag_id, family = BEZI, correlation = corARMA(value = c(0.4276,  -0.9426),~date.num|tag_id, p = 1, q=1))
```
```
GAMLSS-RS iteration 1: Global Deviance = -8493.518
GAMLSS-RS iteration 2: Global Deviance = -8577.033
GAMLSS-RS iteration 3: Global Deviance = -8606.78
GAMLSS-RS iteration 4: Global Deviance = -8616.876
GAMLSS-RS iteration 5: Global Deviance = -8620.174
GAMLSS-RS iteration 6: Global Deviance = -8621.226
GAMLSS-RS iteration 7: Global Deviance = -8621.556
GAMLSS-RS iteration 8: Global Deviance = -8621.659
GAMLSS-RS iteration 9: Global Deviance = -8621.691
GAMLSS-RS iteration 10: Global Deviance = -8621.7
GAMLSS-RS iteration 11: Global Deviance = -8621.703
GAMLSS-RS iteration 12: Global Deviance = -8621.704
```
```
summary(model.ow.f_beta_zi_gamlss)
```
```
******************************************************************
Family:  c("BEZI", "Zero Inflated Beta")

Call:  gamlss(formula = hab.prop ~ day_length * up_lake + tl_mm, family = BEZI,
    data = open_water_prop, random = ~1 | tag_id, correlation = corARMA(value = c(0.4276,
        -0.9426), ~date.num | tag_id, p = 1, q = 1))

Fitting method: RS()

------------------------------------------------------------------
Mu link function:  logit
Mu Coefficients:
                         Estimate Std. Error t value Pr(>|t|)
(Intercept)            -2.7814817  0.3402564  -8.175 5.29e-16 ***
day_length             -0.0490482  0.0223462  -2.195   0.0283 *
up_lakeMost            -0.7547071  0.4861934  -1.552   0.1208
tl_mm                   0.0018209  0.0001886   9.656  < 2e-16 ***
day_length:up_lakeMost  0.0744367  0.0329835   2.257   0.0241 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

------------------------------------------------------------------
Sigma link function:  log
Sigma Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  0.57851    0.02477   23.36   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

------------------------------------------------------------------
Nu link function:  logit
Nu Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)   -24.87    2266.88  -0.011    0.991

------------------------------------------------------------------
No. of observations in the fit:  1946
Degrees of Freedom for the fit:  7
      Residual Deg. of Freedom:  1939
                      at cycle:  12

Global Deviance:     -8621.704
            AIC:     -8607.704
            SBC:     -8568.69
******************************************************************
```

#### 3.3.2. Using package _zoid_

:books:`library(zoib)`
:books:`library(ggplot2)`
:books:`library(plotly)`

The_zoib_ function fits a zero/one inflated beta regression model and obtains the Bayesian Inference for the model via the MCMC approach implemented in JAGS (Liu & Kong, 2015) (see [zoib R documentation](https://cran.r-project.org/web/packages/zoib/zoib.pdf)

```
model.ow.f_beta_zi_zoib1 <- zoib(hab.prop ~ day_length*up_lake+tl_mm | 1 | 1,
                                            data=open_water_prop,
                                            joint = FALSE, random=1, EUID=open_water_prop$tag_id,
                                            zero.inflation = FALSE, one.inflation = FALSE,
                                            n.iter=3200, n.thin=15, n.burn=200)
```

Posterior samples of regression coefficients from the model

```
coeff <- model.ow.f_beta_zi_zoib1$coef
summary(coeff)
```
```
Iterations = 1:200
Thinning interval = 1
Number of chains = 2
Sample size per chain = 200

1. Empirical mean and standard deviation for each variable,
   plus standard error of the mean:

                            Mean       SD  Naive SE Time-series SE
(Intercept)            -3.345490 1.035671 5.178e-02      0.1093600
day_length             -0.066437 0.021998 1.100e-03      0.0019857
up_lakeMost            -2.193605 0.571631 2.858e-02      0.1022484
tl_mm                   0.002631 0.001201 6.004e-05      0.0001378
day_length:up_lakeMost  0.173945 0.031804 1.590e-03      0.0052916
(Intercept)             0.926012 0.039789 1.989e-03      0.0019919
sigma                   0.497519 0.199321 9.966e-03      0.0123749

2. Quantiles for each variable:

                             2.5%       25%       50%       75%     97.5%
(Intercept)            -5.4947877 -3.956496 -3.265265 -2.645798 -1.518609
day_length             -0.1101579 -0.080203 -0.065606 -0.052808 -0.023366
up_lakeMost            -3.3934272 -2.583951 -2.113921 -1.792193 -1.143217
tl_mm                   0.0005704  0.001793  0.002601  0.003373  0.005241
day_length:up_lakeMost  0.1119151  0.151534  0.171972  0.196671  0.234240
(Intercept)             0.8545880  0.897990  0.924192  0.952310  1.003833
sigma                   0.2336854  0.377350  0.453220  0.576084  0.955440
```
```
pred <- pred.zoib(model.ow.f_beta_zi_zoib1, data.frame(temp = seq(100, 600, 0.01)))
```

**Check convergence on the regression coefficients**

```
layout(matrix(1:6, nrow = 2))
traceplot(coeff)
autocorr.plot(coeff)
```
![Hab_prop](/Plots/Hab_prop_zoib1.png "Hab_prop")

![Hab_prop](/Plots/Hab_prop_zoib12.png "Hab_prop")

![Hab_prop](/Plots/Hab_prop_zoib13.png "Hab_prop")

```
check.psrf(coeff)
```
```
                       Point est. Upper C.I.
(Intercept)             1.0397452   1.087617
day_length              0.9964912   0.996502
up_lakeMost             1.0034921   1.003577
tl_mm                   1.0591390   1.178759
day_length:up_lakeMost  1.0005255   1.021111
(Intercept)             0.9966835   1.002369
sigma                   1.0046960   1.025672
[1] 1.067262
$psrf.s
                       Point est. Upper C.I.
(Intercept)             1.0397452   1.087617
day_length              0.9964912   0.996502
up_lakeMost             1.0034921   1.003577
tl_mm                   1.0591390   1.178759
day_length:up_lakeMost  1.0005255   1.021111
(Intercept)             0.9966835   1.002369
sigma                   1.0046960   1.025672

$psrf.m
[1] 1.067262

$psrf.s.summ
        Point est. Upper C.I.
Min.     0.9964912   0.996502
1st Qu.  0.9986045   1.002973
Median   1.0034921   1.021111
Mean     1.0143961   1.045087
3rd Qu.  1.0222206   1.056645
Max.     1.0591390   1.178759
```
![Hab_prop](/Plots/Hab_prop_zoib14.png "Hab_prop")

**Plot posterior mean of _hab.prop_ vs. observed _hab.prop_ to check on goodness of fit**

_ypred_: posterior predictive samples of hab.prop

```
layout(matrix(1:1, nrow = 1))
hab.prop_pred <- rbind(model.ow.f_beta_zi_zoib1$ypred[[1]],model.ow.f_beta_zi_zoib1$ypred[[2]])
post.mean <- apply(hab.prop_pred,2,mean);
plot(open_water_prop$hab.prop, post.mean, col='blue',pch=2);
abline(0,1,col='red')
```
![Hab_prop](/Plots/Hab_prop_zoib15.png "Hab_prop")

**Plot credible intervals of predictions**

xnew <- data.frame(day_length = c(15, 20), tl_mm = c(7, 15), up_lake = factor(c(1, 2), levels = 1:2))
predictions <- pred.zoib(model.ow.f_beta_zi_zoib1, xnew)

predictions <- data.frame(temp = seq(100, 600, 0.01), pred$summary)

ggplotly(
         ggplot() +
         geom_point(data = open_water_prop,
                           aes(x = day_length, y = hab.prop, fill = up_lake),
                           size = 4, shape = 21) +
                           xlim(100, 600) +
         geom_line(data = predictions, aes(y = mean, x = day_length), col="red") +
         geom_ribbon(data = predictions, aes(ymin= X2.5., ymax = X97.5., x = day_length), alpha = 0.3) +
         theme_classic())


```
model.ow.f_beta_zi_zoib2 <- zoib(hab.prop ~ day_length*up_lake+tl_mm|1|day_length*up_lake+tl_mm|1,
                                            data = open_water_prop, random = 1, EUID= open_water_prop$tag_id,
                                            zero.inflation = TRUE, one.inflation = FALSE, joint = FALSE,
                                            n.iter=5000, n.thin=20, n.burn=1000)
```
```
model.ow.f_beta_zi_zoib3 <- zoib(hab.prop ~ day_length*up_lake+tl_mm| 1 | 1,
                                            data=open_water_prop,
                                            joint = FALSE, random=1, EUID=open_water_prop$tag_id,
                                            zero.inflation = FALSE, one.inflation = FALSE,
                                            n.iter=3200, n.thin=15, n.burn=200)
```


### 3.4. Fit a Bayesian zero-inflated beta regression model with MCMC process

#### 3.4.1. Using package _brms_

:books:`library(brms)`

```
model.ow.f_beta_zi_brm1 <- brm(hab.prop ~ day_length*up_lake+tl_mm +(1|tag_id),
                                          data = open_water_prop,
                                          family = zero_one_inflated_beta,
                                          autocor = cor_arma(~date.num|tag_id, cov = TRUE, 1, 1),
                                          warmup = 1000, iter = 2000, cores = 4, chains = 1)
```
```
model.ow.f_beta_zi_brm2 <- brm(hab.prop ~ day_length*up_lake+tl_mm +(1|tag_id),
                                          data = open_water_prop,
                                          family = zero_one_inflated_beta,
                                          arma(~date.num|tag_id, cov = TRUE, p = 1, q = 1),
                                          warmup = 1000, iter = 2000, cores = 4, chains = 1)
```
#### 3.4.2. Using package _rstanarm_

:books:`library(rstanarm)`

```
open_water_prop <- as.data.table(open_water_prop)
```
model.ow.f_beta_zi_stan <- stan_betareg(hab.prop ~ day_length*up_lake+tl_mm +(date.num|tag_id),
                                                   data = open_water_prop,
                                                   link = "logit",
                                                   seed = 12345)

model.ow.f_beta_zi_stan <- stan_betareg(hab.prop ~ day_length*up_lake+tl_mm +(date.num|tag_id),
                                                   data = open_water_prop,
                                                   link = "logit",
                                                   seed = 12345)

model.ow.f_beta_zi_stan <- stan_betareg(hab.prop ~ day_length*up_lake+tl_mm +(date.num|tag_id), data = open_water_prop, link = "logit", link.phi = "log",
                     cores = 4, seed = 12345)
```

## References

Liu, F. & Kong, Y. 2015. ZOIB: an R Package for Bayesian Inferences in Beta and Zero One Inflated Beta Regression Models, The R Journal, 7(2):34-51

Smithson, M. & Verkuilen, J. 2006. "A Better Lemon Squeezer? Maximum-Likelihood Regression with Beta-Distributed Dependent Variables." Psychological Methods, 11 (1), 54–71



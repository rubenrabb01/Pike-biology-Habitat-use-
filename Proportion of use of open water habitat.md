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

**Final model (scaled predictors)**

:books:`library(dplyr)`

```
open_water_prop_scale_p<-open_water_prop %>% mutate_at(c(7,11), funs(c(scale(.))))
```
```
model.ow.f_scale_p <- lme(asin(sqrt(hab.prop)) ~ day_length*up_lake+tl_mm,
                                                 data = open_water_prop_scale_p,
                                                 random = ~1|tag_id, method = "REML",
                                                 correlation = corARMA(value = c(0.4276,  -0.9426),~date.num|tag_id, p = 1, q=1))
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

### 3.2. Fit a beta regression model (distribution with logit/beta processes) using the package _betareg_

:books:`library(betareg)`

- For beta regression (zero/one -inflated) models we need to perform some transformations to the response variable
- First, we need to check if there are NAs for the variable _hab.prop_. We can do this with the following code:

```
open_water_prop[is.na(open_water_prop$hab.prop),]
```
- There is not any NA so we can proceed with the following step
- The variable must be bound between 0 and 1 but excluding the extremes as the _betareg()_ function does not allow these values. We use two approaches:

**Note:** If we dont perform this transformation we get an error due to 0's and 1's contained in the extremes since the betareg function cant handle those values

#### Transform the response variable with values higher than 0 and lower than 1 but in close decimal place

Create a new variable _hab.prop_t_ with 0,1 values changed to 0.001,09999
```
open_water_prop_scale_p <- open_water_prop_scale_p %>% mutate(hab.prop_t = replace(hab.prop, hab.prop == 1, 0.999))
open_water_prop_scale_p <- open_water_prop_scale_p %>% mutate(hab.prop_t = replace(hab.prop_t, hab.prop == 0, 0.001))
```
**Note:** It is important that values are close to 0 and 1 because the odds will be much different (tending to infinite) if we add or substract more decimals (e.g. 0.00001, 0.99999)

**Final model (scaled predictors)**

```
model.ow.f_beta_scaled <- betareg(hab.prop_t ~  day_length*up_lake+tl_mm,
                                  data = open_water_prop_scale_p,
                                  random = ~1|tag_id,
                                  method = "Nelder-Mead",
                                  correlation = corARMA(value = c(0.4276,  -0.9426),~date.num|tag_id, p = 1, q=1))
```
```
summary(model.ow.f_beta_scaled)
```
```

Call:
betareg(formula = hab.prop_t ~ day_length * up_lake + tl_mm, data = open_water_prop_scale_p, random = ~1 |
    tag_id, method = "Nelder-Mead", correlation = corARMA(value = c(0.4276, -0.9426), ~date.num | tag_id,
    p = 1, q = 1))

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

#### Transform the response variable according to Smithson & Verkuilen (2006)

The R documentation for the **_betareg_** package mentions the distribution proposed by these authors, which takes the form:

**y’ = (y*(N - 1) + .5)/N**, where N is the sample size

I create custom function to transform the values of proportion of habitat use (_hab.prop_) following the former formula
```
df <- open_water_prop_scale_p
id <- open_water_prop_scale_p$tag_id
x <- open_water_prop_scale_p$hab.prop

zero_one <- function(x)
{
   N = nrow(df)/(length(levels(id)))
   d = (x*(N - 1) + .5)/N
   return(d)
}
zero_one(x)
```

Create a new variable _hab.prop_tf_ with range of values as provided by the previous function
```
open_water_prop_scale_p <- open_water_prop_scale_p %>% mutate(hab.prop_tf = zero_one(x))
```

**Final model (scaled predictors)**

```
model.ow.f_beta_scaled <- betareg(hab.prop_tf ~ day_length*up_lake+tl_mm,
                                                data = open_water_prop_scale_p,
                                                random = ~1|tag_id, method = "Nelder-Mead",
                                                correlation = corARMA(value = c(0.4276, -0.9426),~date.num|tag_id, p = 1, q=1))
```
```
summary(model.ow.f_beta_scaled)
```
```
Call:
betareg(formula = hab.prop_tf ~ day_length * up_lake + tl_mm, data = open_water_prop_scale_p, random = ~1 | tag_id,
    method = "Nelder-Mead", correlation = corARMA(value = c(0.4276, -0.9426), ~date.num | tag_id, p = 1, q = 1))

Standardized weighted residuals 2:
    Min      1Q  Median      3Q     Max
-1.4543 -0.4904 -0.0946  0.4322  3.1200

Coefficients (mean model with logit link):
                       Estimate Std. Error z value Pr(>|z|)
(Intercept)            -1.92357    0.04198 -45.820  < 2e-16 ***
day_length             -0.07949    0.03308  -2.403   0.0163 *
up_lakeMost             0.33857    0.05077   6.669 2.58e-11 ***
tl_mm                   0.24825    0.02561   9.692  < 2e-16 ***
day_length:up_lakeMost  0.11838    0.04837   2.448   0.0144 *

Phi coefficients (precision model with identity link):
      Estimate Std. Error z value Pr(>|z|)
(phi)   2.4997     0.0882   28.34   <2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Type of estimator: ML (maximum likelihood)
Log-likelihood:  2827 on 6 Df
Pseudo R-squared: 0.1923
Number of iterations: 1739 (Nelder-Mead) + 2 (Fisher scoring)
```

### 3.3. Fit a zero- /one- inflated beta regression model

Now we fit a series of zero-inflated models using different packages on the original dataset without further transformation of the (0,1) bound

#### 3.3.1. Using package _gamlss_

:books:`library(gamlss)`

The GAMLSS R documentation refers to the fitting of GAM models for Location Scale and Shape (Rigby and Stasinopoulos, 2005). The models use a distributional regression approach where all the parameters of the conditional distribution of the response variable are modelled using explanatory variables (see [zoib R documentation](https://cran.r-project.org/web/packages/gamlss/gamlss.pdf))

##### 3.3.1.1. Fit final model

- We specify family _BEZI_ instead of _BEINF_ since if use the non-transformed (0,1) distribution there is an error of "response variable out of range"
- If instead we use the transformed variable _hab_prop_t_ the model fit well but it is not allright for a zero/one inflated beta distribution

The GAMLSS documentantion states that:
The function BEINF() defines the beta inflated distribution, a four parameter distribution, for a gamlss.family object to be used in GAMLSS fitting using the function gamlss(). The beta inflated is similar to the beta but allows zeros and ones as values for the response variable. The two
extra parameters model the probabilities at zero and one.

Thus, we use famiLY _BEZI_ for fitting a zero-inflated distribution and family _BEINF_ for a zero/one -inflated distribution (see [zoib R documentation, p. 20](https://cran.r-project.org/web/packages/gamlss.dist/gamlss.dist.pdf))

```
model.ow.f_beta_zi_gamlss <- gamlss(hab.prop ~ day_length*up_lake+tl_mm,
                                                  data = open_water_prop_scale_p,
                                                  random = ~1|tag_id,
                                                  family = BEINF,
                                                  correlation = corARMA(value = c(0.4276, -0.9426),~date.num|tag_id, p = 1, q=1))
```
```
GAMLSS-RS iteration 1: Global Deviance = 2069.06
GAMLSS-RS iteration 2: Global Deviance = 2038.689
GAMLSS-RS iteration 3: Global Deviance = 2038.665
GAMLSS-RS iteration 4: Global Deviance = 2038.664
```
```
summary(model.ow.f_beta_zi_gamlss)
```
```
******************************************************************
Family:  c("BEINF", "Beta Inflated")

Call:  gamlss(formula = hab.prop ~ day_length * up_lake +
    tl_mm, family = BEINF, data = open_water_prop_scale_p,
    random = ~1 | tag_id, correlation = corARMA(value = c(0.4276,
        -0.9426), ~date.num | tag_id, p = 1, q = 1))

Fitting method: RS()

------------------------------------------------------------------
Mu link function:  logit
Mu Coefficients:
                       Estimate Std. Error t value Pr(>|t|)
(Intercept)            -1.37949    0.06337 -21.769  < 2e-16 ***
day_length             -0.17988    0.05219  -3.447 0.000579 ***
up_lakeMost             0.42789    0.07431   5.758 9.85e-09 ***
tl_mm                   0.25839    0.03773   6.848 1.00e-11 ***
day_length:up_lakeMost  0.23869    0.06811   3.504 0.000468 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

------------------------------------------------------------------
Sigma link function:  logit
Sigma Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  0.06401    0.03426   1.868   0.0619 .
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

------------------------------------------------------------------
Nu link function:  log
Nu Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  0.22198    0.04573   4.854 1.31e-06 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

------------------------------------------------------------------
Tau link function:  log
Tau Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  -4.4555     0.3181  -14.01   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

------------------------------------------------------------------
No. of observations in the fit:  1946
Degrees of Freedom for the fit:  8
      Residual Deg. of Freedom:  1938
                      at cycle:  4

Global Deviance:     2038.664
            AIC:     2054.664
            SBC:     2099.252
******************************************************************
```

**Plot model**

```
plot(model.ow.f_beta_zi_gamlss)
```
```
******************************************************************
	 Summary of the Randomised Quantile Residuals
                           mean   =  -0.00845896
                       variance   =  0.9614488
               coef. of skewness  =  0.02864597
               coef. of kurtosis  =  3.196871
Filliben correlation coefficient  =  0.9992626
******************************************************************
```
![Hab_prop](/Plots/Hab_prop_gamlss01.png "Hab_prop")

**Plot the ACF and PACF of the residuals**

```
acfResid(model.ow.f_beta_zi_gamlss)
```
![Hab_prop](/Plots/Hab_prop_gamlss02.png "Hab_prop")

**Plots the centile curves**

```
centiles(model.ow.f_beta_zi_gamlss,xvar=open_water_prop$hab.prop)
```
```
% of cases below  0.4 centile is  55.24152
% of cases below  2 centile is  55.24152
% of cases below  10 centile is  55.24152
% of cases below  25 centile is  55.24152
% of cases below  50 centile is  55.24152
% of cases below  75 centile is  77.74923
% of cases below  90 centile is  90.69887
% of cases below  98 centile is  97.79034
% of cases below  99.6 centile is  100
```
![Hab_prop](/Plots/Hab_prop_gamlss03.png "Hab_prop")

##### 3.3.1.2. Fit final model with cublic splines for _day_length_

```
model.ow.f_beta_zi_gamlss_c.spline <- gamlss(hab.prop ~  scs(day_length, by="up_lake")+tl_mm,
                                      data = open_water_prop_scale_p,
                                      random = ~1|tag_id,
                                      family = BEINF,
                                      correlation = corARMA(value = c(0.4276,  -0.9426),~date.num|tag_id, p = 1, q=1))
```
```
GAMLSS-RS iteration 1: Global Deviance = 2089.491
GAMLSS-RS iteration 2: Global Deviance = 2059.017
GAMLSS-RS iteration 3: Global Deviance = 2059.049
GAMLSS-RS iteration 4: Global Deviance = 2059.079
GAMLSS-RS iteration 5: Global Deviance = 2059.084
GAMLSS-RS iteration 6: Global Deviance = 2059.085
GAMLSS-RS iteration 7: Global Deviance = 2059.085
```
```
summary(model.ow.f_beta_zi_gamlss_c.spline)
```
```
******************************************************************
Family:  c("BEINF", "Beta Inflated")

Call:  gamlss(formula = hab.prop ~ scs(day_length, by = "up_lake") +
    tl_mm, family = BEINF, data = open_water_prop_scale_p,
    random = ~1 | tag_id, correlation = corARMA(value = c(0.4276,
        -0.9426), ~date.num | tag_id, p = 1, q = 1))

Fitting method: RS()

------------------------------------------------------------------
Mu link function:  logit
Mu Coefficients:
                                Estimate Std. Error t value Pr(>|t|)
(Intercept)                     -1.10235    0.04173 -26.418   <2e-16 ***
scs(day_length, by = "up_lake") -0.04332    0.03386  -1.279    0.201
tl_mm                            0.27717    0.03638   7.618    4e-14 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

------------------------------------------------------------------
Sigma link function:  logit
Sigma Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  0.08522    0.03411   2.499   0.0125 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

------------------------------------------------------------------
Nu link function:  log
Nu Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  0.22198    0.04573   4.854 1.31e-06 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

------------------------------------------------------------------
Tau link function:  log
Tau Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  -4.4555     0.3181  -14.01   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

------------------------------------------------------------------
NOTE: Additive smoothing terms exist in the formulas:
 i) Std. Error for smoothers are for the linear effect only.
ii) Std. Error for the linear terms maybe are not accurate.
------------------------------------------------------------------
No. of observations in the fit:  1946
Degrees of Freedom for the fit:  21.94018
      Residual Deg. of Freedom:  1924.06
                      at cycle:  7

Global Deviance:     2059.085
            AIC:     2102.966
            SBC:     2225.25
******************************************************************
```

**Plot model**

```
plot(model.ow.f_beta_zi_gamlss_c.spline)
```
```
******************************************************************
	 Summary of the Randomised Quantile Residuals
                           mean   =  -0.006583857
                       variance   =  0.9630601
               coef. of skewness  =  0.04277548
               coef. of kurtosis  =  3.11282
Filliben correlation coefficient  =  0.9991167
******************************************************************
```
![Hab_prop](/Plots/Hab_prop_gamlss04.png "Hab_prop")

**Plot the ACF and PACF of the residuals**

```
acfResid(model.ow.f_beta_zi_gamlss_c.spline)
```
![Hab_prop](/Plots/Hab_prop_gamlss05.png "Hab_prop")

**Plots the centile curves**

```
centiles(model.ow.f_beta_zi_gamlss_c.spline,xvar=open_water_prop$hab.prop)
```
```
% of cases below  0.4 centile is  55.24152
% of cases below  2 centile is  55.24152
% of cases below  10 centile is  55.24152
% of cases below  25 centile is  55.24152
% of cases below  50 centile is  55.24152
% of cases below  75 centile is  76.92703
% of cases below  90 centile is  91.00719
% of cases below  98 centile is  97.73895
% of cases below  99.6 centile is  100
```
![Hab_prop](/Plots/Hab_prop_gamlss06.png "Hab_prop")

**Plot fitted model values (day length)**

```
data_ow.f_beta_zi_gamlss_c.spline <- visreg(model.ow.f_beta_zi_gamlss_c.spline, breaks = c ("Chabarovice", "Most"), gg = TRUE, overlay = TRUE, jitter = TRUE, lwd = 0.5, rug = FALSE, partial = FALSE, plot = FALSE)
layout(matrix(1:2, nrow = 1))
plot(data_ow.f_beta_zi_gamlss_c.spline, plot.type="rgl") + theme_bw()
```
![Hab_prop](/Plots/Hab_prop_gamlss07.png "Hab_prop")

#### 3.3.2. Using package _zoid_

:books:`library(zoib)`
:books:`library(ggplot2)`
:books:`library(plotly)`

The **_zoib_** function fits a zero-one-inflated regression model and obtains the Bayesian Inference for the model via the MCMC approach implemented in JAGS (Liu & Kong, 2015) (see [zoib R documentation](https://cran.r-project.org/web/packages/zoib/zoib.pdf))

```
model.ow.f_beta_zi_zoib1 <- zoib(hab.prop ~ day_length*up_lake+tl_mm | 1 | 1,
                                            data=open_water_prop_scale_p,
                                            joint = FALSE, random=1, EUID=open_water_prop_scale_p$tag_id,
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

                       Mean       SD       Naive SE Time-series SE
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


### 3.4. Fit a Bayesian zero-one-inflated beta regression model with MCMC process

#### 3.4.1. Using package _brms_

:books:`library(brms)`

We specify family **_zero_one_inflated_beta_**

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
```
model.ow.f_beta_zi_stan <- stan_betareg(hab.prop ~ day_length*up_lake+tl_mm +(date.num|tag_id),
                                                   data = open_water_prop,
                                                   link = "logit",
                                                   seed = 12345)

model.ow.f_beta_zi_stan <- stan_betareg(hab.prop ~ day_length*up_lake+tl_mm +(date.num|tag_id),
                                                   data = open_water_prop,
                                                   link = "logit",
                                                   seed = 12345)

model.ow.f_beta_zi_stan <- stan_betareg(hab.prop ~ day_length*up_lake+tl_mm +(date.num|tag_id), data = open_water_prop, link = "logit", link.phi = "log", cores = 4, seed = 12345)
```

## References

- _Liu, F. & Kong, Y_. 2015. ZOIB: an R Package for Bayesian Inferences in Beta and Zero One Inflated Beta Regression Models, The R Journal, 7(2):34-51

- _Rigby, R. A. & Stasinopoulos, D. M_. 2005. Generalized additive models for location, scale and shape,(with discussion), Appl. Statist., 54, part 3, pp 507-554

- _Smithson, M. & Verkuilen, J_. 2006. A Better Lemon Squeezer? Maximum-Likelihood Regression with Beta-Distributed Dependent Variables. Psychological Methods, 11 (1), 54–71


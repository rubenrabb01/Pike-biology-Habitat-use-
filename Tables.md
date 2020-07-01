---
output: pdf_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE)
library(stargazer)
```

Here is the latex table in a PDF document:

```{r mylatextable, results = "asis"}
stargazer(m_vert_r_7,m_vert_r_2,m_vert_r_4,m_vert_r_1,m_vert_r_3,
          dep.var.labels=  "Vertical range use: log(V-KS)",
          covariate.labels=c("(Incercept)","Daily temperature × Lake(LSCL)","Day Length × Lake(LSCL)","Day Length × Body Length","Daily temperature","Day Length","Lake(LSCL)","Body Length"),
          type="latex",ci = F,intercept.bottom = T, align=TRUE, order=c("Constant","temp.day:up_lake0","day_length:up_lake0","day_length:tl_mm","temp.day","day_length","up_lake0","tl_mm"),title="Vertical range")
```

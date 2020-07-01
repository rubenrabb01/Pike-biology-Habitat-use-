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
% Table created by stargazer v.5.2.2 by Marek Hlavac, Harvard University. E-mail: hlavac at fas.harvard.edu
% Date and time: Wed, Jul 01, 2020 - 09:12:06 AM
\begin{table}[!htbp] \centering 
  \caption{} 
  \label{} 
\begin{tabular}{@{\extracolsep{5pt}}lc} 
\\[-1.8ex]\hline 
\hline \\[-1.8ex] 
 & \multicolumn{1}{c}{\textit{Dependent variable:}} \\ 
\cline{2-2} 
\\[-1.8ex] & log(k1size) \\ 
\hline \\[-1.8ex] 
 day\_length & 0.099$^{***}$ \\ 
  & (0.031) \\ 
  & \\ 
 tl\_mm & $-$0.020 \\ 
  & (0.066) \\ 
  & \\ 
 temp.day & $-$0.060$^{**}$ \\ 
  & (0.030) \\ 
  & \\ 
 up\_lake0 & 0.165 \\ 
  & (0.129) \\ 
  & \\ 
 day\_length:tl\_mm & 0.104$^{***}$ \\ 
  & (0.031) \\ 
  & \\ 
 Constant & 5.188$^{***}$ \\ 
  & (0.090) \\ 
  & \\ 
\hline \\[-1.8ex] 
Observations & 1,946 \\ 
Log Likelihood & $-$2,309.173 \\ 
Akaike Inf. Crit. & 4,638.345 \\ 
Bayesian Inf. Crit. & 4,694.080 \\ 
\hline 
\hline \\[-1.8ex] 
\textit{Note:}  & \multicolumn{1}{r}{$^{*}$p$<$0.1; $^{**}$p$<$0.05; $^{***}$p$<$0.01} \\ 
\end{tabular} 
\end{table} 

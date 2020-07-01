```{r, results='asis'}

stargazer(m_vert_r_7,m_vert_r_2,m_vert_r_4,m_vert_r_1,m_vert_r_3,
          dep.var.labels=  "Vertical range use: log(V-KS)",
          covariate.labels=c("(Incercept)","Daily temperature × Lake(LSCL)","Day Length × Lake(LSCL)","Day Length × Body Length","Daily temperature","Day Length","Lake(LSCL)","Body Length"),
          ci = F,intercept.bottom = T, align=TRUE, order=c("Constant","temp.day:up_lake0","day_length:up_lake0","day_length:tl_mm","temp.day","day_length","up_lake0","tl_mm"),title="Vertical range")
```
% Table created by stargazer v.5.2.2 by Marek Hlavac, Harvard University. E-mail: hlavac at fas.harvard.edu
% Date and time: Wed, Jul 01, 2020 - 09:12:35 AM
% Requires LaTeX packages: dcolumn 
\begin{table}[!htbp] \centering 
  \caption{Vertical range} 
  \label{} 
\begin{tabular}{@{\extracolsep{5pt}}lD{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} } 
\\[-1.8ex]\hline 
\hline \\[-1.8ex] 
 & \multicolumn{5}{c}{\textit{Dependent variable:}} \\ 
\cline{2-6} 
\\[-1.8ex] & \multicolumn{5}{c}{Vertical range use: log(V-KS)} \\ 
\\[-1.8ex] & \multicolumn{1}{c}{(1)} & \multicolumn{1}{c}{(2)} & \multicolumn{1}{c}{(3)} & \multicolumn{1}{c}{(4)} & \multicolumn{1}{c}{(5)}\\ 
\hline \\[-1.8ex] 
 (Incercept) & 5.188^{***} & 5.186^{***} & 5.187^{***} & 5.185^{***} & 5.185^{***} \\ 
  & (0.090) & (0.090) & (0.091) & (0.090) & (0.091) \\ 
  & & & & & \\ 
 Daily temperature × Lake(LSCL) &  & -0.048 &  & -0.048 &  \\ 
  &  & (0.060) &  & (0.060) &  \\ 
  & & & & & \\ 
 Day Length × Lake(LSCL) &  &  & -0.012 & -0.012 & -0.015 \\ 
  &  &  & (0.066) & (0.066) & (0.066) \\ 
  & & & & & \\ 
 Day Length × Body Length & 0.104^{***} & 0.105^{***} & 0.107^{***} & 0.107^{***} & 0.109^{***} \\ 
  & (0.031) & (0.031) & (0.033) & (0.033) & (0.033) \\ 
  & & & & & \\ 
 Daily temperature & -0.060^{**} & -0.036 & -0.060^{**} & -0.035 &  \\ 
  & (0.030) & (0.043) & (0.030) & (0.043) &  \\ 
  & & & & & \\ 
 Day Length & 0.099^{***} & 0.099^{***} & 0.105^{**} & 0.105^{**} & 0.106^{**} \\ 
  & (0.031) & (0.031) & (0.045) & (0.045) & (0.045) \\ 
  & & & & & \\ 
 Lake(LSCL) & 0.165 & 0.165 & 0.167 & 0.166 & 0.176 \\ 
  & (0.129) & (0.128) & (0.129) & (0.129) & (0.130) \\ 
  & & & & & \\ 
 Body Length & -0.020 & -0.022 & -0.021 & -0.023 & -0.022 \\ 
  & (0.066) & (0.065) & (0.066) & (0.065) & (0.066) \\ 
  & & & & & \\ 
\hline \\[-1.8ex] 
Observations & \multicolumn{1}{c}{1,946} & \multicolumn{1}{c}{1,946} & \multicolumn{1}{c}{1,946} & \multicolumn{1}{c}{1,946} & \multicolumn{1}{c}{1,946} \\ 
Log Likelihood & \multicolumn{1}{c}{-2,309.173} & \multicolumn{1}{c}{-2,308.853} & \multicolumn{1}{c}{-2,309.157} & \multicolumn{1}{c}{-2,308.836} & \multicolumn{1}{c}{-2,311.165} \\ 
Akaike Inf. Crit. & \multicolumn{1}{c}{4,638.345} & \multicolumn{1}{c}{4,639.705} & \multicolumn{1}{c}{4,640.314} & \multicolumn{1}{c}{4,641.673} & \multicolumn{1}{c}{4,642.330} \\ 
Bayesian Inf. Crit. & \multicolumn{1}{c}{4,694.080} & \multicolumn{1}{c}{4,701.014} & \multicolumn{1}{c}{4,701.623} & \multicolumn{1}{c}{4,708.555} & \multicolumn{1}{c}{4,698.066} \\ 
\hline 
\hline \\[-1.8ex] 
\textit{Note:}  & \multicolumn{5}{r}{$^{*}$p$<$0.1; $^{**}$p$<$0.05; $^{***}$p$<$0.01} \\ 
\end{tabular} 
\end{table} 

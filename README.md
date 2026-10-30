# BRI_assets_model
R code and supporting data for evaluating bounded and unbounded growth models using the long-term asset trajectory of PT Bank Rakyat Indonesia (Persero) Tbk.

Models examined:
- Linear
- Quadratic
- Exponential
- Verhulst Logistic
- Gompertz
- Richards

The repository includes:
- Expanding-window out-of-sample validation\
- Full-sample model calibration
- Information-theoretic model comparison (AIC, BIC, ΔAIC, Akaike weights)
- Recursive parameter-stability analysis
- Residual diagnostics
- Figure and table generation

The empirical data consist of quarterly consolidated assets (IDR million) covering 2013Q1–2026Q1 and were compiled from publicly available financial reports released by PT Bank Rakyat Indonesia (Persero) Tbk.:
https://bri.co.id/web/guest/report-detail

## Requirements 
- R
- minpack.lm
- ggplot2
- stats
- patchwork

## Running the Code Execute the main R script to reproduce all analyses and figures. 

## Author
Dhidhi Pambudi,
Vernanda Kisty Faradilla,
Sachiko Nakano,
Aivin B. Balquedra

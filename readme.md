# ARMA System Identification of Realized Volatility Across Market Regimes

This repository contains the MATLAB implementation of a system identification and econometric analysis of monthly **Realized Volatility (RV)** across six major global equity indices from **2000 to 2026**. 

The project investigates volatility dynamics and persistence across distinct macroeconomic regimes, comparing regime-specific models against a long-run cross-period pooled benchmark.

---

## Dataset & Preprocessing

Daily closing prices are filtered for missing/invalid observations, sorted chronologically, and converted into percentage daily log-returns:
$$r_t = 100 \cdot \left(\log P_t - \log P_{t-1}\right)$$

Monthly realized volatility is aggregated across common trading months and log-transformed:
$$\text{RV}_m = \sqrt{\sum_{t \in m} r_t^2}, \qquad y_m = \log(\text{RV}_m + \varepsilon)$$

### Analyzed Equity Indices
* **TSX** (`^tsx_d.csv`) — S&P/TSX Composite Index (Canada)
* **DAX** (`^dax_d.csv`) — DAX Performance Index (Germany)
* **Nikkei 225** (`^nkx_d.csv`) — Nikkei 225 (Japan)
* **KOSPI** (`^kospi_d.csv`) — KOSPI Composite Index (South Korea)
* **Nasdaq 100** (`^ndx_d.csv`) — Nasdaq 100 (USA)
* **FTSE-MIB** (`^fmib_d.csv`) — FTSE MIB (Italy)

### Macroeconomic Regimes
The 2000–2026 historical timeline is partitioned into four sub-periods:
1. **P1 (2000–2006):** Post dot-com recovery
2. **P2 (2007–2010):** Global Financial Crisis
3. **P3 (2011–2017):** Post-crisis quantitative easing and bull market
4. **P4 (2018–2026):** COVID-19 shock and rate-hiking cycle

---

## Methodology

1. **System Identification & Order Selection:**
   * Estimation of $\text{ARMA}(n_a, n_c)$ models via Prediction Error Methods (PEM) using `armax`.
   * Order selection carried out over the grid $n_a \in \{1,\dots,6\}$ and $n_c \in \{0,\dots,4\}$.
   * Multi-criteria comparison between **BIC (Bayesian Information Criterion)** on the full regime segment and **Temporal Hold-out Cross-Validation** (70% training, 30% validation).

2. **Cross-Period (Pooled) vs. Regime-Specific Fitting:**
   * Compares models fitted within isolated regimes against a global cross-period pooled model.
   * Compares residual innovation noise variance across regimes and the pooled benchmark.

3. **Dynamical Analysis & Persistence:**
   * Analysis of system stability via pole-zero constellation maps in the complex plane.
   * Tracking of the dominant pole magnitude $\rho_{\max} = \max_i |p_i|$ and the primary autoregressive parameter $\phi_1 = -a_1$ across time regimes to assess structural changes in volatility memory.

---

## Requirements & Usage

* **MATLAB** with the **System Identification Toolbox** (`iddata`, `armax`, `compare`).

### Running the Analysis
1. Place the six index CSV files (`^tsx_d.csv`, `^dax_d.csv`, `^nkx_d.csv`, `^kospi_d.csv`, `^ndx_d.csv`, `^fmib_d.csv`) in the working directory.
2. Run the main script:
   ```matlab
   run('main_rv_identification.m')
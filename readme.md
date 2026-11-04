# ARMA System Identification of Realized Volatility Across Market Regimes

This repository contains the MATLAB implementation of a system identification and econometric analysis of monthly **Realized Volatility (RV)** across six major global equity indices from **2000 to 2026**. 

The project investigates volatility persistence across distinct macroeconomic regimes, testing whether high long-run autocorrelation in pooled models is structural or partially an artifact of unmodelled level shifts between market regimes (**spurious persistence**).

---

## Dataset & Market Regimes

Daily closing prices are processed to compute log-returns and aggregated into monthly realized volatility (normalized by the number of trading days per month to remove deterministic seasonal effects):

$$\text{RV}_m = \sqrt{\frac{1}{N_m} \sum_{t \in m} r_t^2}, \qquad y_m = \log(\text{RV}_m + \varepsilon)$$

### Analyzed Equity Indices
* **TSX** (`^tsx_d.csv`) — S&P/TSX Composite Index (Canada)
* **DAX** (`^dax_d.csv`) — DAX Performance Index (Germany)
* **Nikkei 225** (`^nkx_d.csv`) — Nikkei 225 (Japan)
* **KOSPI** (`^kospi_d.csv`) — KOSPI Composite Index (South Korea)
* **Nasdaq 100** (`^ndx_d.csv`) — Nasdaq 100 (USA)
* **FTSE-MIB** (`^fmib_d.csv`) — FTSE MIB (Italy)

### Macroeconomic Regimes
The 2000–2026 sample is partitioned into four historical regimes:
1. **P1 (2000–2006):** Dot-com bubble aftermath and recovery
2. **P2 (2007–2010):** Global Financial Crisis and European sovereign debt crisis
3. **P3 (2011–2017):** Post-crisis quantitative easing and bull market
4. **P4 (2018–2026):** COVID-19 shock, inflation surge, and interest rate hikes

---

## Methodology & Key Features

1. **Multi-Criteria Order Selection ($\text{ARMA}(n_a, n_c)$):**
   * **BIC (Bayesian Information Criterion)** over the grid $n_a \in \{1,\dots,6\}$, $n_c \in \{0,\dots,4\}$.
   * **GCV (Generalized Cross-Validation)** to approximate leave-one-out error without refitting $N$ models.
   * **Exact PRESS (Leave-One-Out CV)** for pure AR models via the **Woodbury matrix identity** leveraging the hat matrix diagonal $h_t = \phi_t^\top (\Phi^\top \Phi)^{-1} \phi_t$, avoiding sample starvation on short regime segments.
   * **Hold-out Cross-Validation** (70/30 split) as a benchmark.

2. **Model Diagnostics & Over-Parameterization Checks:**
   * **Residual Whiteness:** 99% confidence bounds ($\pm 2.58/\sqrt{N}$) on the residual ACF and **Ljung–Box portmanteau test** ($p$-value via Wilson–Hilferty $\chi^2$ approximation).
   * **Pole–Zero Cancellation Detection:** Automatic flagging of near-cancellations ($\min |p_i - z_j| < \text{tol}$) and cross-market pole dispersion analysis across regimes.
   * **Finite-Sample Bias Correction:** Analytical correction of the $\text{AR}(1)$ persistence parameter $\hat{a}_1$ via the Marriott–Pope / Kendall formula:
     $$\hat{a}_1^{\text{corr}} = \hat{a}_1 + \frac{1 + 3\hat{a}_1}{N}$$

3. **Spurious Persistence & Multi-Experiment Pooling:**
   * Compares a standard **Global Pooled ARMA** model (single global mean) against a **Demeaned Multi-Experiment Pooled ARMA** model (using MATLAB's `merge` across regime-demeaned `iddata` objects so regressors never cross regime boundaries).
   * Quantifies how much of the apparent persistence $\hat{a}_1$ and innovation variance in the full sample is driven by regime level shifts rather than true intra-regime dynamics.

4. **Time-Varying Parameter Tracking (RPEM / Forgetting-Factor RLS):**
   * Implements a **Recursive Prediction Error Method / Recursive Least Squares (RLS)** estimator with forgetting factor $\lambda = 0.97$ (effective memory $\approx 33$ months) on the augmented regressor $\varphi(t) = [y(t-1),\, 1]^\top$ to jointly track time-varying persistence $a_1(t)$ and local volatility level without prior detrending.

---

## Requirements & Usage

* **MATLAB** with the **System Identification Toolbox** (`iddata`, `armax`, `resid`, `merge`, `getpvec`).

### Running the Code
1. Clone this repository and ensure the six index CSV files (`^tsx_d.csv`, `^dax_d.csv`, `^nkx_d.csv`, `^kospi_d.csv`, `^ndx_d.csv`, `^fmib_d.csv`) are in the root directory.
2. Run the main script in MATLAB:
   ```matlab
   run('main_rv_identification.m')
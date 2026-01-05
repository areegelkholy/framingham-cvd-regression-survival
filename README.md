# Cardiovascular Disease Risk Modeling  
### Logistic Regression, Survival Analysis, and Simulation  
**Framingham Heart Study**

## Overview
This project investigates the relationship between cardiovascular disease (CVD) and major clinical and demographic risk factors using data from the Framingham Heart Study. The analysis combines exploratory data analysis, logistic regression, survival models, and simulation-based theoretical validation to study both the **occurrence** and **timing** of CVD events.

The work was completed as a graded exam group project for the **Regression (2025)** course at the **University of Copenhagen**.

---

## Data
- **Dataset**: Framingham Heart Study  
- **Source**: `riskCommunicator` R package  
- **Sample**: 4,434 individuals, up to 3 examination periods, 24-year follow-up  
- **Outcomes**:
  - `CVD`: Binary indicator of cardiovascular disease
  - `TIMECVD`: Time to first CVD event (right-censored)

---

## Methods
The project is structured into three main components:

### 1. Exploratory Data Analysis
- Assessment of missingness and structural missing data
- Baseline distributions of clinical and demographic predictors
- Correlation analysis and predictor screening
- Consistency and plausibility checks

### 2. Binary Regression Models
- Logistic regression models for 24-year CVD risk
- Variable selection using likelihood ratio tests
- Model diagnostics (residuals, calibration)
- Predictive performance evaluation (AUC ≈ 0.75)
- Exploration of interaction terms and non-linear effects

### 3. Survival Analysis & Theory
- Cox proportional hazards models for `TIMECVD`
- Assessment of proportional hazards assumptions
- Theoretical link between logistic regression at fixed time horizons and survival models
- Simulation study validating:
  - Time-invariant logistic slopes
  - Logarithmic time-dependent intercepts
  - Bias under covariate-dependent censoring

---

## Key Findings
- Age, sex, smoking, hypertension, diabetes, and prior cardiac conditions are strong predictors of CVD
- Logistic models showed good calibration and moderate predictive power
- Cox regression confirmed consistent effects in continuous time
- Simulation demonstrated that logistic regression fails under informative censoring, while survival models remain robust



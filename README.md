# Cardiovascular Disease Risk Modeling  
### Logistic Regression, Survival Analysis, and Simulation  
**Framingham Heart Study**

## Project Overview
This project analyzes cardiovascular disease (CVD) risk using the Framingham Heart Study, combining **statistical modeling, survival analysis, and simulation-based theory validation**.  
The objective is to understand both **whether** and **when** individuals develop CVD, and to examine the relationship between fixed-horizon logistic models and continuous-time survival models.

The work was completed as a graded **exam project** for the *Regression (2025)* course at the **University of Copenhagen**.

---

## Data
- **Source**: Framingham Heart Study (`riskCommunicator` R package)
- **Cohort**: 4,434 individuals, followed for up to 24 years
- **Structure**: Longitudinal examinations (1–3 observations per subject)
- **Outcomes**:
  - `CVD`: Binary indicator of cardiovascular disease
  - `TIMECVD`: Time to first CVD event (right-censored)

---

## Methodology

### 1. Exploratory Data Analysis
- Missingness assessment and identification of structural missing data
- Baseline-only cohort construction to avoid repeated-measure bias
- Distributional analysis of clinical and demographic predictors
- Correlation screening and biologically motivated variable selection
- Data consistency and plausibility checks

### 2. Binary Risk Modeling
- Logistic regression models for 24-year CVD risk
- Model refinement using likelihood ratio tests
- Diagnostics via residual analysis and calibration plots
- Predictive performance assessment (AUC ≈ 0.75)
- Evaluation of interaction effects and non-linear terms

### 3. Survival Analysis & Theory
- Cox proportional hazards models for time-to-CVD
- Assessment of proportional hazards assumptions
- Analytical derivation linking:
  - Fixed-time logistic regression  
  - Underlying survival functions and hazards
- Simulation study validating:
  - Time-invariant logistic slopes
  - Logarithmic time-dependent intercepts
  - Bias introduced by covariate-dependent censoring

---

## Key Findings
- Age, sex, smoking, hypertension, diabetes, and prior cardiac conditions are dominant CVD risk factors
- Logistic regression models showed good calibration and moderate discrimination
- Survival models confirmed consistent effects in continuous time
- Simulation demonstrated that logistic regression fails under informative censoring, while Cox models remain robust

---

## Technical Stack
- **Language**: R  
- **Models**: Logistic regression, Cox proportional hazards  
- **Methods**: Likelihood ratio testing, cross-validation, simulation via inverse transform sampling  
- **Visualization**: ggplot2  

---

## Individual Contribution
My contributions included:
- Exploratory data analysis and missing-data strategy
- Logistic regression modeling and diagnostics
- Survival analysis and theoretical derivations
- Design and implementation of the simulation study
- Interpretation and academic writing of analytical sections

---



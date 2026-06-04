# DiD-Toolbox for Matlab

## Overview

Code and Description for a Difference-in-Differences (DiD) Toolbox for Matlab.

**Requirements**

* Matlab from Version 2021+

* Statistics and Machine Learning Toolbox

## Purpose

The DiD Toolbox is a set of Matlab tools designed for applied statisticians and econometricians to conduct Difference-in-Differences (DiD) analyses, particularly focusing on designs involving staggered treatment timing. The primary goal of this toolbox is to address these methodological challenges by providing modern, robust estimators that yield valid causal estimates in complex multi-period, staggered adoption scenarios.

Implemented estimators are based on:

* **Goodman-Bacon (2021):** Provides the mathematical foundation for the decomposition of the TWFE estimator, explaining how it operates as a weighted average of DiD estimates and identifying the source of bias from "negative weights" due to treatment effect heterogeneity.

* **Wooldridge (2021):** Establishes the algebraic equivalence between the Two-Way Fixed Effects (TWFE) estimator and the Two-Way Mundlak (TWM) regression, enabling flexible implementation using pooled OLS.

* **Borusyak, Jaravel, and Spiess (2024):** Derives the efficient and robust imputation estimator (BJS) for staggered DiD, which estimates counterfactual outcomes using only untreated observations to calculate heterogeneous causal effects, providing efficiency and avoiding spurious identification.

* **de Chaisemartin and D'Haultfoeuille (2020):** Proposes the DID-M estimator, which estimates a robust Average Treatment Effect across switching cells and introduces robustness measures for assessing TWFE bias, particularly in designs where weights may be negative. A unique feature is that the estimator can handle an on/off treatment, while all other estimators assume that treatment is an absorptive state.

* **Callaway / Sant'Anna (2021):** CS focus on cohort/time based analyses, looking at cohort-wise treatment effects with a focus on taking covariates into account. They derive identification, estimation, and inference strategies assuming parametric nuisance models (i.e. linear outcome regression, logit for propensity scores) that admit standard regularity conditions. In an extension, the interaction weighted estimator by Sun/Abraham (2021) builds on this work, showing how to get unbiased event studies in this setting.

* **Rambachan/Roth (2023):** The study suggests a sensitivity analysis for the parallel trends assumption of DiD analyses by asking "how large could a violation be before our conclusion changes?" They propose bounding the treatment effect estimates under various degrees of pre-trend differences.

* **Arkhangelsky et al. (2021):** The study discusses synthetic difference-in-differences (SDID) methods, which generate parallel trends by reweighting units to match their pre-exposure trends.

## Install

**Recommended (toolbox installer)**: Download the latest MATLAB toolbox package:

[DIDToolbox.mltbx](https://github.com/relsas/DiDToolbox/releases/latest/download/DIDToolbox.mltbx)

Open the downloaded `.mltbx` file in MATLAB to install the toolbox.

**Alternative (source)**: Clone/download the repository and in MATLAB run:

```matlab
addpath(genpath('/path/to/DIDToolbox'))
```

The source tree is useful for inspecting the code, examples, tests, and documentation. For normal installation, the `.mltbx` file is the intended download.




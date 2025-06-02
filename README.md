# JMbdirect: Joint Model for Longitudinal and Multiple Time-to-Event Data
<img src="man/figures/JMbdirectlogo.png" align="right" alt="JMbdirect logo" width="180"/>

<!-- badges: start -->

[![CRAN Status](https://www.r-pkg.org/badges/version/JMbdirect)](https://CRAN.R-project.org/package=JMbdirect)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/JMbdirect)](https://cran.r-project.org/package=JMbdirect)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://cran.r-project.org/web/licenses/GPL-3)

<!-- badges: end -->

> **JMbdirect** provides tools for jointly modeling longitudinal outcomes and multiple time-to-event (survival) outcomes using a variety of statistical engines like **JMbayes2**, **FastJM**, **joineRML**, and **rstanarm**.

---

## ✨ Features

- 🔄 Fit **joint models** with **multiple survival outcomes**
- 🧠 Compatible with popular modeling engines:
  - `JMbayes2`
  - `FastJM`
  - `joineRML`
  - `rstanarm`
- 🔎 Includes functions for:
  - Model fitting (`jmbB()`, `jmcsB()`, `jmrmlB()`, `jmstB()`)
  - Survival prediction
  - Posterior plots
  - Bootstrapped confidence intervals
- 📦 Ready for **big data applications** using the `BIGdata = TRUE` flag

---

## 📦 Installation

### From GitHub (development version)

```r
# If not already installed
install.packages("remotes")
remotes::install_github("kumarbhrigu/JMbdirect")
```

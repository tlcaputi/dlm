# Installation

## Stata

### From GitHub (recommended)

```stata
net install dlm, from("https://raw.githubusercontent.com/tlcaputi/dlm-stata/main/")
```

This installs `dlm.ado`, `dlm_gen_data.ado`, and the help file.

### Prerequisites

The Stata package requires:

- **Stata 15** or later
- **reghdfe** — install with `ssc install reghdfe`

### Verify installation

```stata
which dlm
help dlm
```

### Uninstall

```stata
ado uninstall dlm
```

---

## R

### From GitHub

```r
# install.packages("devtools")  # if needed
devtools::install_github("tlcaputi/dlm")
```

### Prerequisites

The R package requires:

- **R 4.0+**
- **fixest** — installed automatically as a dependency
- **dplyr**, **ggplot2**, **glue**, **logger**, **scales** — installed automatically

### Verify installation

```r
library(dlm)
?distributed_lags_model
```

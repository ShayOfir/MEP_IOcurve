# MEP I/O Curve Analysis - Hierarchical Bayesian Modeling

Bayesian hierarchical sigmoidal model for comparing TMS coil efficacy in motor evoked potential (MEP) recruitment, with explicit handling of equipment-induced missing-not-at-random (MNAR) censoring.

## Overview

Motor Evoked Potentials (MEPs) measure cortical excitability by recording muscle responses to Transcranial Magnetic Stimulation (TMS). The recruitment curve—MEP amplitude as a function of stimulation intensity—provides crucial insights into corticospinal excitability.

There present Repository contains the code and data of a hierarchical Bayesian model for analyzing MEP recruitment curves across different TMS coil configurations, for a study aimed to comparatively assess the neurophysiological properties of motor evoked potentials (MEPs) elicited from the first dorsal interosseous (FDI) muscle following uni-directioma; TMS (via Figure-of-8 and H7 coils) vs. multidirectional rfTMS. In this study, 10 healthy adult subjects received TMS via the three coil configurations in a random order.

In short, each subject underwent resting motor threshold (rMT) determination followed by MEP recruitment curve acquisition at multiple stimulation intensities (from 90% to 150% of rMT for that side, coil and subject). The MEP amplitudes were recorded from the contralateral FDI muscle using surface electromyography (EMG). Approximately 5 trials were collected for each intensity. Condition order was randomize to counter-balance potential order effects.

More detail will be provided in the upcoming publication:

*Wonderman Bar Sela et al., 2026: "MEPs elicited by multidirectional rotational-field TMS show marked differences compared to unidirectional Figure-of-8 and H7 coils",* PLoS ONE\* (in review).\*

This repository will focus on the methodological aspects of modeling MEP recruitment curves with Bayesian hierarchical models.

### Why baysian heirarchial modelling?

Modeling of MEP amplitudes in general and in the current study faces several challenges: 1. There is a considerable inter-subject variability 2. The relationship between stimulation intensity and MEP amplitude is sigmoidal, requiring non-linear modeling 3. There is an inherent censoring bias, when sampling the logisitic curve at high intensities, due to safety limitations of TMS stimulation (one cannot exceed the maximal stimulator output).

Therefore, this data is both structurally-complex, suffer from structural missingness, and requires non-linear modeling. Bayesian hierarchical modeling provides a principled framework to address these challenges by allowing partial pooling of information across subjects and conditions, flexible non-linear model specification, and explicit modeling of the censoring mechanism. Moreover, bayesian inference yields full posterior distributions for parameters of interest, allowing quantification of uncertainty in estimates and predictions.

The analysis follows modern Bayesian workflow principles [Gelman et al., 2020](https://arxiv.org/abs/2011.01808), including comprehensive prior/posterior predictive checks and simulation-based parameter recovery validation. All code and pre-fitted models are provided for full reproducibility.

## Features

-   **Hierarchical Bayesian sigmoidal model** for MEP recruitment curves with non-centered parameterization
-   **Explicit MNAR censoring mechanism** accounting for coil-specific equipment limits
-   **Weakly informative priors** calibrated through iterative prior predictive checks
-   **Comprehensive model validation**: convergence diagnostics (R̂, ESS, divergences), posterior predictive checks
-   **Simulation-based sensitivity analysis**: Type I error control and parameter recovery
-   **Reproducible workflow**: pre-compiled Stan models with SHA-256 integrity verification, pre-fitted results
-   **Full transparency**: all data, code, and intermediate results included

## Repository Structure

```         
├── Main.ipynb                      # Main analysis notebook
├── IOcurve_model_functions.R       # Model fitting and data processing functions
├── make_hash.R                     # File integrity verification
├── CHECKSUMS.md                    # SHA-256 hashes for compiled models
├── data/
│   ├── cleaned_MEP_data.csv        # MEP amplitude data
│   └── motor_thresholds.csv        # rMT measurements
├── stan_models/
│   └── IOcurve9.stan               # Stan model specification
├── compiled_models/
│   └── IOcurve9.exe                # Pre-compiled Stan model (Windows)
├── fitted_models/
│   ├── prefitted_actual.RDS        # Fitted model (actual data)
│   ├── sim0_fit.RDS                # Null model simulation
│   └── sim1_fit.RDS                # Full effects simulation
```

## Requirements

### R Packages

``` r
required_packages <- c(
  "cmdstanr", "posterior", "dplyr", "ggplot2", 
  "tidyr", "rethinking", "purrr", "tibble", 
  "vctrs", "bayesplot", "dagitty", "digest", "knitr"
)
```

### System Requirements

-   R \>= 4.5.1
-   **For model compilation (optional)**:
    -   Windows: [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
    -   Mac/Linux: C++ toolchain

**Note**: Pre-compiled models are provided to avoid compilation requirements.

## Installation

### 1. Clone Repository with Git LFS

``` sh
# Clone the repository
git clone https://github.com/yourusername/MEP_IOcurve.git
cd MEP_IOcurve

# Install Git LFS (if not already installed)
git lfs install

# Pull large files (fitted models, results)
git lfs pull
```

### 2. Install R Dependencies

Open R/RStudio/Positron and run:

``` r
# The first cell in Main.ipynb will automatically install all required packages
# Or manually install:

required_packages <- c(
  "cmdstanr", "posterior", "dplyr", "ggplot2", 
  "tidyr", "rethinking", "purrr", "tibble", 
  "vctrs", "bayesplot", "dagitty", "digest", "knitr"
)

# Install from CRAN
install.packages(setdiff(required_packages, "cmdstanr"))

# Install cmdstanr from GitHub
if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("stan-dev/cmdstanr")

# Install CmdStan (for model compilation, optional if using pre-compiled models)
cmdstanr::install_cmdstan()
```

Note: Installing CmdStan may take some time as it compiles the CmdStan binaries. Compiling stan models requires a working C++ toolchain, with correcr configuration, which is not always trivial. Please refer to [CmdStanR installation guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html#installation) for troubleshooting.

### 3. Verify Installation

``` r
source("IOcurve_model_functions.R")

# Verify compiled model integrity
verify_file_hash("compiled_models/IOcurve9.exe")
# Should return: TRUE
```

## Usage

### Quick Start

Open `Main.ipynb` in Positron/RStudio and run cells sequentially.

### Key Analysis Steps

1.  **Load data and functions**

    ``` r
    source('IOcurve_model_functions.R')
    p2p <- load_and_preprocess_MEP_data()
    ```

2.  **Prior predictive checks**

    ``` r
    pp <- prior_predictive(mod, data_actual, priors)
    ```

3.  **Fit the model**

    ``` r
    # Use pre-fitted results
    fit <- readRDS("fitted_models/prefitted_actual.RDS")

    # Or refit (takes ~30 min)
    FIT_MODEL <- TRUE  # in Main.ipynb
    ```

4.  **View results**

    ``` r
    res <- show_results(fit, data_actual, ci = 0.89)
    ```

## Model Specification

The model uses a **hierarchical sigmoidal function** to describe MEP recruitment:

**Equation 1: Sigmoidal recruitment curve**

```         
μ(i,j,k) = A(i,j,k) / (1 + exp(-s(i,j,k) · (I - θ(i,j,k))))
```

Where: - `A`: Plateau amplitude (maximal recruitment) - `s`: Slope (recruitment curve steepness) - `θ`: Threshold (intensity for 50% maximal response) - `i`: Subject, `j`: Coil type, `k`: Hemisphere

**Equation 2: Log-normal likelihood**

```         
Y ~ LogNormal(ln(μ), σ_ε)
```

**Equations 3-4: Non-centered parameterization**

```         
φ(i,j,k) = β₀ + δ_coil(j) + δ_side(k) + u(i)
u(i) ~ Normal(0, τ)
```

Where `φ ∈ {ln(A), θ, s}` with sum-to-zero constraints on δ effects.

**Equations 5-6: MNAR censoring**

```         
Only include observations where: I ≤ MSO_limit(j)
Y | I ≤ MSO_limit ~ LogNormal(ln(μ), σ_ε)
```

This explicitly models the truncation mechanism, avoiding selection bias from equipment limits.

See `Main.ipynb` sections 1-3 for complete mathematical specification and derivations.

# ...existing code...

## Analysis Pipeline

The complete Bayesian workflow includes:

1.  **Data Preparation** (Section 4 in Main.ipynb)
    -   Load and preprocess MEP amplitude data
    -   Prepare rMT (resting motor threshold) measurements
    -   Format data for Stan model
2.  **Model Compilation & Verification** (Section 4)
    -   Verify integrity of pre-compiled model (SHA-256 checksums)
    -   Or compile from source using RTools/C++ toolchain
3.  **Prior Predictive Checks** (Section 4)
    -   Sample from prior distribution only
    -   Verify that the majority of prior mass overlaps observed data
    -   Confirm weakly informative priors allow data to dominate
4.  **Model Fitting** (Section 6)
    -   HMC sampling: 4 chains × 2000 iterations (1000 warmup)
    -   Non-centered parameterization for efficiency
    -   Adaptive delta = 0.999, max tree depth = 20
5.  **Convergence Diagnostics** (Section 6)
    -   R̂ \< 1.01, ESS \> 400 for all parameters
    -   Zero divergences, no max tree depth hits
    -   Traceplots
6.  **Posterior Predictive Checks** (Section 6)
    -   Generate replicated data from posterior
    -   Verify coverage of observed data
7.  **Inference** (Section 7)
    -   Extract posterior distributions for θ, s, log(A)
    -   Compute contrasts between coil types and hemispheres
    -   89% HPDIs for all parameters
8.  **Sensitivity Analysis** (Section 8)
    -   Simulation-based parameter recovery
    -   Type I error control (null model)
    -   Statistical power assessment (full effects model)

## Security & Reproducibility

### Compiled Model Verification

Pre-compiled models include SHA-256 checksums for integrity verification:

``` r
verify_file_hash("compiled_models/IOcurve9.exe")
```

If verification fails, recompile from source:

``` r
mod <- cmdstan_model("stan_models/IOcurve9.stan")
```

### Reproducibility

-   All fitted models saved as `.RDS` files
-   Package versions recorded in session info

## Troubleshooting

### Common Issues

**1. Git LFS files not downloading**

``` sh
# Ensure Git LFS is installed
git lfs install

# Manually pull large files
git lfs pull
```

**2. Model verification fails**

``` r
# If SHA-256 hash doesn't match, recompile from source
mod <- cmdstan_model("stan_models/IOcurve9.stan")

# Generate new checksum
source("make_hash.R")
compute_and_save_hash("compiled_models/IOcurve9.exe")
```

**3. CmdStan compilation errors** - **Windows**: Ensure [Rtools](https://cran.r-project.org/bin/windows/Rtools/) is installed and on PATH - **Mac/Linux**: Install C++ toolchain (gcc/clang) - See [CmdStanR troubleshooting guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html#troubleshooting)

**4. Memory issues during model fitting**

``` r
# Reduce parallel chains
sampler_params <- make_sampler_params_list(parallel_chains = 2)

# Or fit chains sequentially
sampler_params <- make_sampler_params_list(parallel_chains = 1)
```

**5. Package installation failures**

``` r
# For cmdstanr, ensure remotes is installed
install.packages("remotes")
remotes::install_github("stan-dev/cmdstanr")

# For rethinking, ensure dependencies are met
install.packages(c("coda", "mvtnorm", "loo"))
```

## Development

### Repository Maintenance

**Adding new compiled models:**

``` r
# Compile model
new_mod <- cmdstan_model("stan_models/new_model.stan")

# Generate checksum
source("make_hash.R")
compute_and_save_hash("compiled_models/new_model.exe")
```

**Updating fitted models:**

``` r
# Fit model
fit <- fit_TMS_model(mod, data_list, priors_list, sampler_params,
                     out_file = "fitted_models/my_fit.RDS")

# Track with Git LFS before committing
```

**Before committing large files:**

``` sh
# Ensure Git LFS is tracking the file type
git lfs track "*.RDS"
git lfs track "*.exe"

# Verify tracking
git lfs ls-files

# Then commit normally
git add .gitattributes
git commit -m "Configure LFS tracking"
git add fitted_models/*.RDS
git commit -m "Add fitted models via LFS"
```

### File Organization

-   **`data/`**: Raw and preprocessed data (CSV format)
-   **`stan_models/`**: Stan model specifications (.stan files)
-   **`compiled_models/`**: Pre-compiled Stan executables (.exe)
-   **`fitted_models/`**: Saved model fits (.RDS files, tracked by Git LFS)
-   **`results/`**: Output figures and tables
-   **`IOcurve_model_functions.R`**: Core modeling functions
-   **`make_hash.R`**: Security/integrity verification utilities

## Citation

If you use this code in your research, please cite:

*(to be added upon publication)*

## References

-   Gelman et al. (2020). Bayesian Workflow. arXiv:2011.01808
-   Koponen et al. (2024). Sigmoidal modeling of MEP recruitment curves. Front. Hum. Neurosci.
-   Stan Development Team (2024). Stan Modeling Language Users Guide

## License

BDSD 3-Clause License. See `LICENSE` file for details.

## Contact

-   **Author**: Shay Ofir-Geva
-   **Email**: shinofir\@gmail.com; ofirsh3\@clalit.org.il
-   **Issues**: [GitHub Issues](https://github.com/yourusername/MEP_IOcurve/issues)

## Acknowledgments

This repository was built with assistence of Claude Sonnet 4.5 and ChatGPT 5.1 via GitHub Copilot. All AI-generated content has been reviewed and approved by the author.
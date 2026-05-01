# ============================================================================
# WEEK 1 - R Package Structure
# Run this script to generate complete project structure
# ============================================================================


########################################
# 0. Install packages and load libraries
########################################

#library(pacman)
#pacman::p_load("fs", "usethis", "devtools")

install.packages(c("cli", "usethis", "devtools", "fs"))
library(usethis)
library(devtools)
library(cli)


#########################################
## 1. Create a package
#########################################

# Create a project, set up folder structure, open and switch to new RStudio session inside package
create_package("~/Documents/R/r-packages-by-me/musicianEMG")


#########################################
# 2. Setup infrastructure
#########################################

# Git
use_git()

# Testing
use_testthat()
use_mit_license()  # Or other choice

# Documentation
use_roxygen_md()
use_readme_rmd()
use_news_md()

# Package dependencies
use_package("dplyr")
use_package("tidyr")
use_package("purrr")
use_package("ggplot2")
use_package("signal")
use_package("R.matlab")
use_pipe() # Enables %>%


#########################################
# 3. Create function files
#########################################

# File: R/emg_processing.R
use_r("emg_processing")
# Paste this content:

cat('
#\' Bandpass filter for EMG signals
#\'
#\' Applies 4th order Butterworth bandpass filter following ISEK standards
#\'
#\' @param signal Numeric vector of raw EMG amplitude
#\' @param fs Sampling frequency in Hz
#\' @param low Low cutoff frequency (Hz), default 20
#\' @param high High cutoff frequency (Hz), default 450
#\' @param order Filter order, default 4
#\'
#\' @return Filtered signal (same length as input)
#\' @export
#\'
#\' @examples
#\' raw_emg <- rnorm(1000)
#\' filtered <- bandpass_emg(raw_emg, fs = 1000)
bandpass_emg <- function(signal, fs, low = 20, high = 450, order = 4) {

  # Input validation
  stopifnot(
    "Signal must be numeric vector" = is.numeric(signal) && is.vector(signal),
    "Signal cannot be empty" = length(signal) > 0,
    "Signal must not contain NA/Inf" = all(is.finite(signal)),
    "Sampling rate must be positive" = fs > 0,
    "Cutoff frequencies invalid" = low > 0 && high < fs/2 && low < high
  )

  # Design Butterworth filter
  bf <- signal::butter(order, c(low, high) / (fs/2), type = "pass")

  # Zero-phase filtering
  filtered <- signal::filtfilt(bf, signal)

  return(as.numeric(filtered))
}

#\' Calculate RMS envelope
#\'
#\' @param signal Numeric vector (typically rectified EMG)
#\' @param window Window size in samples
#\'
#\' @return RMS envelope
#\' @export
calculate_rms <- function(signal, window = 100) {

  stopifnot(
    is.numeric(signal),
    length(signal) > 0,
    window > 0,
    window <= length(signal)
  )

  # Moving RMS calculation
  n <- length(signal)
  rms <- numeric(n)

  for(i in 1:n) {
    start_idx <- max(1, i - window + 1)
    end_idx <- i
    window_vals <- signal[start_idx:end_idx]
    rms[i] <- sqrt(mean(window_vals^2))
  }

  return(rms)
}

#\' Complete EMG processing pipeline
#\'
#\' Bandpass filter -> Rectify -> RMS envelope
#\'
#\' @param signal Raw EMG signal
#\' @param fs Sampling frequency
#\' @param low Low cutoff (Hz)
#\' @param high High cutoff (Hz)
#\' @param rms_window RMS window size (samples)
#\'
#\' @return Processed RMS signal
#\' @export
#\'
#\' @examples
#\' emg_processed <- process_emg(raw_signal, fs = 1000)
process_emg <- function(signal, fs = 1000,
                        low = 20, high = 450,
                        rms_window = 100) {

  signal %>%
    bandpass_emg(fs = fs, low = low, high = high) %>%
    abs() %>%
    calculate_rms(window = rms_window)
}
', file = "R/emg_processing.R")

# File: R/emg_normalisation.R
use_r("emg_#########################################")

cat('
#\' Normalize EMG to MVC
#\'
#\' @param signal Processed EMG signal
#\' @param mvc_value Maximum voluntary contraction value
#\'
#\' @return Normalized EMG (% MVC)
#\' @export
normalize_mvc <- function(signal, mvc_value) {
  stopifnot(
    is.numeric(signal),
    is.numeric(mvc_value),
    mvc_value > 0
  )

  (signal / mvc_value) * 100
}

#\' Normalize to peak activation
#\'
#\' @param signal EMG signal
#\'
#\' @return Normalized to 0-100%
#\' @export
normalize_peak <- function(signal) {
  (signal / max(signal, na.rm = TRUE)) * 100
}
', file = "R/emg_normalisation.R")

# File: R/emg_analysis.R
use_r("emg_analysis")

cat('
#\' Calculate median frequency (fatigue indicator)
#\'
#\' @param signal Raw EMG signal
#\' @param fs Sampling frequency
#\'
#\' @return Median frequency (Hz)
#\' @export
median_frequency <- function(signal, fs) {

  # Power spectral density
  psd <- signal::pwelch(signal, fs = fs)

  # Find median
  cumsum_power <- cumsum(psd$spec)
  total_power <- sum(psd$spec)
  median_idx <- which(cumsum_power >= total_power/2)[1]

  psd$freq[median_idx]
}

#\' Muscle coactivation index
#\'
#\' @param agonist Agonist muscle EMG
#\' @param antagonist Antagonist muscle EMG
#\'
#\' @return Coactivation index (%)
#\' @export
coactivation_index <- function(agonist, antagonist) {

  # Normalize both to 0-1
  ag_norm <- (agonist - min(agonist)) / (max(agonist) - min(agonist))
  ant_norm <- (antagonist - min(antagonist)) / (max(antagonist) - min(antagonist))

  # Overlap area
  overlap <- pmin(ag_norm, ant_norm)
  (sum(overlap) / length(overlap)) * 100
}
', file = "R/emg_analysis.R")


#########################################
# 4. Create tests
#########################################

use_test("emg_processing")

cat('
test_that("bandpass_emg works with valid input", {
  signal <- sin(2*pi*50*seq(0, 1, length.out = 1000))
  result <- bandpass_emg(signal, fs = 1000)

  expect_type(result, "double")
  expect_equal(length(result), 1000)
  expect_true(all(is.finite(result)))
})

test_that("bandpass_emg rejects invalid input", {
  expect_error(bandpass_emg(NULL, fs = 1000))
  expect_error(bandpass_emg(c(1, NA, 3), fs = 1000))
  expect_error(bandpass_emg(c(1, 2, 3), fs = -1000))
})

test_that("process_emg pipeline works", {
  signal <- rnorm(1000)
  result <- process_emg(signal, fs = 1000)

  expect_type(result, "double")
  expect_equal(length(result), 1000)
  expect_true(all(result >= 0))  # RMS always positive
})
', file = "tests/testthat/test-emg_processing.R")


#########################################
# 5. Document
#########################################
document()


#########################################
# 6. Check package
#########################################
check()


#########################################
# 7. Build read-me file
#########################################
cat('
# musicianEMG

EMG analysis tools for musician health research.

## Installation

```r
# install.packages("devtools")
devtools::install_github("yourusername/musicianEMG")
```

## Quick Start

```r
library(musicianEMG)

# Load EMG data
raw_emg <- read_emg("path/to/file.mat")

# Process
processed <- process_emg(raw_emg, fs = 1000)

# Normalize
normalized <- normalize_mvc(processed, mvc_value = 150)

# Analyze fatigue
med_freq <- median_frequency(raw_emg, fs = 1000)
```

## Features

- Bandpass filtering (ISEK standards)
- RMS envelope calculation
- MVC normalisation
- Fatigue indices
- Coactivation analysis

', file = "README.md")


#########################################
# 8. Initial git commit
#########################################

system("git add .")
system("git commit -m \"Initial package structure\"")

print("✓ R package created: ~/musicianEMG")
print("✓ Run devtools::check() to verify")
print("✓ Start adding your functions to R/ directory")

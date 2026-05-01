#' Bandpass filter for EMG signals
#'
#' Applies a Butterworth bandpass filter following typical EMG standards.
#'
#' @param signal Numeric vector of raw EMG amplitude
#' @param fs Sampling frequency in Hz
#' @param low Low cutoff frequency (Hz), default 20
#' @param high High cutoff frequency (Hz), default 450
#' @param order Filter order, default 4
#'
#' @return Filtered signal (same length as input)
#' @export
#'
#' @examples
#' raw_emg <- rnorm(1000)
#' filtered <- bandpass_emg(raw_emg, fs = 1000)
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
  bf <- signal::butter(order, c(low, high) / (fs / 2), type = "pass")

  # Zero-phase filtering
  filtered <- signal::filtfilt(bf, signal)

  return(as.numeric(filtered))
}


#' Calculate RMS envelope
#'
#' Computes moving RMS of a signal.
#'
#' @param signal Numeric vector (typically rectified EMG)
#' @param window Window size in samples
#'
#' @return RMS envelope (same length as input)
#' @export
calculate_rms <- function(signal, window = 100) {

  stopifnot(
    is.numeric(signal),
    length(signal) > 0,
    window > 0,
    window <= length(signal)
  )

  n <- length(signal)
  rms <- numeric(n)

  for (i in seq_len(n)) {
    start_idx <- max(1, i - window + 1)
    window_vals <- signal[start_idx:i]
    rms[i] <- sqrt(mean(window_vals^2))
  }

  return(rms)
}


#' Complete EMG processing pipeline
#'
#' Bandpass filter -> Rectify -> RMS envelope
#'
#' @param signal Raw EMG signal
#' @param fs Sampling frequency
#' @param low Low cutoff (Hz)
#' @param high High cutoff (Hz)
#' @param rms_window RMS window size (samples)
#'
#' @return Processed RMS signal
#' @export
#'
#' @examples
#' raw <- rnorm(1000)
#' processed <- process_emg(raw, fs = 1000)
process_emg <- function(signal, fs = 1000,
                        low = 20, high = 450,
                        rms_window = 100) {

  filtered <- bandpass_emg(signal, fs = fs, low = low, high = high)
  rectified <- abs(filtered)
  rms <- calculate_rms(rectified, window = rms_window)

  return(rms)
}

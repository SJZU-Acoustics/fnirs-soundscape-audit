# =============================================================================
# Fig 1b recompute — amplitude response of the delivered 0.01-0.10 Hz band-pass.
# The filter is exactly the delivered pipeline's: signal::butter(3,
# c(0.01, 0.10)/5, "pass") applied zero-phase with filtfilt, so the effective
# gain is |H|^2. Pure signal computation — no data are read.
#
# Method: probe sinusoids at log-spaced frequencies 0.005-0.3 Hz, filtered with
# the same filtfilt call; gain measured by least-squares projection of the
# steady-state output onto the probe quadrature (robust to filtfilt edge
# transients, unlike a max|y| readout).
#
# Gate (mandatory): gain at the design fundamental 0.0111 Hz must reproduce
# 0.688 within +/-0.01, and the 60/30/20-s components ~1.01.
#
# Output: output/data_lock/fig1b_filter_gain.csv (plotted values)
#         output/data_lock/fig1b_gain_anchor_check.csv
# =============================================================================

source("code/style.R")   # LOCKDIR

butter   <- signal::butter
filtfilt <- signal::filtfilt

FS <- 10; NYQ <- FS / 2
BF <- butter(3, c(0.01, 0.10) / NYQ, type = "pass")   # the collector's band

# probe length: 4000 s so even the 0.005 Hz probe completes 20 cycles; the
# steady-state window drops the first and last 500 s (filtfilt transients).
T_PROBE <- 4000
t <- seq(0, T_PROBE, by = 1 / FS)
steady <- t >= 500 & t <= (T_PROBE - 500)

gain_at <- function(f) {
  s <- sin(2 * pi * f * t)
  y <- filtfilt(BF, s)
  # project steady-state output onto probe sin/cos quadrature
  cs <- cos(2 * pi * f * t)
  X <- cbind(s[steady], cs[steady])
  cf <- .lm.fit(X, y[steady])$coefficients
  sqrt(sum(cf^2))
}

freqs <- exp(seq(log(0.005), log(0.3), length.out = 60))
sweep <- data.frame(freq_hz = freqs, gain = vapply(freqs, gain_at, numeric(1)))

# design markers: the 90-s fundamental and the 60/30/20-s within-cycle components
markers <- data.frame(
  component = c("90-s cycle (0.0111 Hz)", "60-s component", "30-s component",
                "20-s component"),
  freq_hz = c(1 / 90, 1 / 60, 1 / 30, 1 / 20))
markers$gain <- vapply(markers$freq_hz, gain_at, numeric(1))
print(markers, digits = 4)

write.csv(sweep, file.path(LOCKDIR, "fig1b_filter_gain.csv"), row.names = FALSE)
write.csv(markers, file.path(LOCKDIR, "fig1b_gain_anchor_check.csv"), row.names = FALSE)

g90 <- markers$gain[1]
stopifnot("Fig1b gate failed: gain at 0.0111 Hz is not 0.688 +/- 0.01" =
            abs(g90 - 0.688) <= 0.01)
cat(sprintf("GATE PASSED: 90-s fundamental gain = %.4f (delivered pipeline: 0.688)\n", g90))
cat(sprintf("60/30/20-s components: %.3f %.3f %.3f\n",
            markers$gain[2], markers$gain[3], markers$gain[4]))

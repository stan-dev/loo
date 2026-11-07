library(ggplot2)
set.seed(1)

scale <- match.arg(commandArgs(TRUE)[1], c("relative", "absolute"))

mean_ci <- function(x, level = 0.8) {
  boot <- replicate(2000, mean(sample(x, replace = TRUE)))
  tail <- (1 - level) / 2
  c(mean = mean(x), lo = quantile(boot, tail, names = FALSE), hi = quantile(boot, 1 - tail, names = FALSE))
}

measurement_files <- function(dir) {
  files <- list.files(dir, pattern = "^measurements.csv$", recursive = TRUE, full.names = TRUE)
  setNames(files, basename(dirname(files)))
}

summarise_runs <- function(seconds, benchmark, impl, workers, reference) {
  ci <- mean_ci(seconds / reference)
  data.frame(benchmark = benchmark, impl = impl, workers = workers, mean = ci[["mean"]], lo = ci[["lo"]], hi = ci[["hi"]])
}

scaling_points <- function(platform, workers, reference_file) {
  rows <- list()
  for (w in workers) {
    for (file in measurement_files(file.path(platform, paste0("w", w)))) {
      name <- basename(dirname(file))
      reference_path <- reference_file(name)
      if (!file.exists(reference_path)) next
      reference <- if (scale == "relative") mean(read.csv(reference_path)$baseline_seconds) else 1
      m <- read.csv(file)
      for (impl in c("baseline", "candidate")) {
        rows[[length(rows) + 1]] <- summarise_runs(m[[paste0(impl, "_seconds")]], name, impl, w, reference)
      }
    }
  }
  do.call(rbind, rows)
}

serial_points <- function(benchmarks) {
  do.call(rbind, lapply(benchmarks, function(name) {
    m <- read.csv(file.path("linux/serial1/w1", name, "measurements.csv"))
    reference <- if (scale == "relative") mean(m$baseline_seconds) else 1
    rbind(
      summarise_runs(m$baseline_seconds, name, "baseline", 1, reference),
      summarise_runs(m$candidate_seconds, name, "candidate", 1, reference)
    )
  }))
}

draw <- function(d, file, reference_label) {
  d$impl <- factor(d$impl, c("baseline", "candidate"), c("Before (cores)", "This PR (parallel)"))
  ylab <- if (scale == "relative") paste("Time relative to", reference_label, "(below 1 = faster)") else "Seconds (mean, 80% bootstrap interval)"
  p <- ggplot(d, aes(workers, mean, colour = impl, fill = impl)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, colour = NA) +
    geom_line() +
    geom_point(size = 0.8) +
    facet_wrap(~benchmark, scales = "free_y", ncol = 6) +
    scale_x_continuous(trans = "log2", breaks = unique(d$workers)) +
    labs(x = "Workers", y = ylab, colour = NULL, fill = NULL) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom", panel.spacing = unit(0.8, "lines"))
  if (scale == "relative") p <- p + geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40")
  ggsave(file, p, width = 16, height = 11, dpi = 130)
}

order_by_serial_time <- function(d, serial_seconds) {
  d$benchmark <- factor(d$benchmark, names(sort(serial_seconds, decreasing = TRUE)))
  d
}

mean_baseline <- function(path) mean(read.csv(path)$baseline_seconds)
serial_path <- function(name) file.path("linux/serial1/w1", name, "measurements.csv")
first_windows_path <- function(name) file.path("windows/w2", name, "measurements.csv")

linux <- scaling_points("linux", c(2, 4, 8, 24, 32), serial_path)
benchmarks <- unique(linux$benchmark)
linux <- rbind(linux, serial_points(benchmarks))
serial_seconds <- sapply(benchmarks, function(name) mean_baseline(serial_path(name)))
draw(order_by_serial_time(linux, serial_seconds), paste0("scaling_linux_", scale, ".png"), "sequential `lapply`")

windows <- scaling_points("windows", c(2, 4, 8), first_windows_path)
benchmarks <- unique(windows$benchmark)
windows_seconds <- sapply(benchmarks, function(name) mean_baseline(first_windows_path(name)))
draw(order_by_serial_time(windows, windows_seconds), paste0("scaling_windows_", scale, ".png"), "the old code with 2 workers")

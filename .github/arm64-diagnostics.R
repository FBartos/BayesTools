# Temporary Apple silicon diagnostics (removed before the pull request is ready).
library(BayesTools)
cat("long double digits:", .Machine$longdouble.digits, " capabilities:", capabilities("long.double"), "\n")
hex <- function(x) paste(sprintf("%a", x), collapse = " ")

cat("\n== sd() and mean() paths of the centered predictor ==\n")
set.seed(5)
x <- rnorm(36, mean = 10, sd = 2)
xc <- x - mean(x)
cat("mean:", hex(mean(xc)), " mean na.rm:", hex(mean(xc, na.rm = TRUE)), "\n")
cat("sd:", hex(stats::sd(xc)), " sd na.rm:", hex(stats::sd(xc, na.rm = TRUE)), "\n")

cat("\n== truncated exponential draws (prior-exp-2) ==\n")
set.seed(1)
draws <- rng(prior("exp", list(2), list(1, 3)), 100000)
cat("sum:", hex(sum(draws)), " first:", hex(draws[1:3]), "\n")
cat("hist counts:", paste(hist(draws, breaks = 50, plot = FALSE)$counts, collapse = " "), "\n")

cat("\n== ordered singular-allocation prior curve ==\n")
common <- new.env()
withr::with_dir("tests/testthat", sys.source("common-functions.R", envir = common))
fixture <- common$ordered_plot_test_fixture(prior("normal", list(0, 1)))
cat("draws digest:", hex(sum(fixture$draws)), "\n")
p <- plot_posterior(samples = fixture$samples, parameter = "mu_f", prior = TRUE,
  plot_type = "ggplot", n_points = 64L, n_samples = 128L)
built <- ggplot2::ggplot_build(p)
for (i in seq_along(built$data)) {
  d <- built$data[[i]]
  if (!all(c("x", "y") %in% names(d))) next
  cat("layer", i, "n", nrow(d), "x sum", hex(sum(d$x)), "y sum", hex(sum(d$y)), "\n")
}
d <- built$data[[2]]
print(data.frame(x = sprintf("%.17g", d$x[28:38]), y = sprintf("%.17g", d$y[28:38])))

cat("\n== JAGS reading of numeric constants ==\n")
texts <- c("1.9599639845400538", "123456789.12345679", "1.3229858936207395", "1e-300")
syntax <- paste0("model{\n", paste0("  x[", seq_along(texts), "] <- ", texts, "\n", collapse = ""), "  y ~ dnorm(0, 1)\n}")
model <- rjags::jags.model(textConnection(syntax), n.chains = 1, quiet = TRUE)
jags <- rjags::coda.samples(model, "x", n.iter = 1, progress.bar = "none")[[1L]][1L, paste0("x[", seq_along(texts), "]")]
print(data.frame(text = texts, jags = sprintf("%a", jags),
  strtod = sprintf("%a", BayesTools:::.JAGS_read_number(texts)), r = sprintf("%a", as.numeric(texts))))
cat("JAGS equals strtod reader:", identical(unname(jags), BayesTools:::.JAGS_read_number(texts)), "\n")

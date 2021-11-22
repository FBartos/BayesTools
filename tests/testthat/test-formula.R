

df <- data.frame(
  y      = rnorm(60),
  x_cont = rnorm(60),
  x_bin  = rbinom(60, 1, .5),
  x_fac3 = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C")),
  x_fac4 = factor(rep(c("A", "B", "C", "D"), 15), levels = c("A", "B", "C", "D"))
)
str(df)


summary(lm(y ~ x_cont + x_bin + x_fac3 + x_fac4, data = df))


str(model.frame(y ~ x_cont + x_bin + x_fac3 + x_fac4, data = df))
str(model.frame(y ~ x_cont - x_bin + x_fac3 + x_fac4, data = df))

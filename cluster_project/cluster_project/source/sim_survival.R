# n: sample size
# lambda: mean survival time
sim_survival = function(n, lambda){
  # true underlying survival time
  t = rexp(n, 1/lambda)

  # Compute censoring time threshold to achieve desired censoring proportion
  c = runif(n, min = min(t),
                         max = 500)

  y = pmin(t, c)
  delta = as.numeric(t <= c)

  tibble(t = t, c = c, y = y, delta = delta)
}

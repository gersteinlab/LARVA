# Copyright (c) 2008 Rick Wash <rwash@umich.edu>
# 
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
# 
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#
# powerlaw.R
#
# Contains code for the continuous and discrete powerlaw distributions
# Also contains a model fit for powerlaw data
# Mostly based on a paper from Mark Newman
#  M.EE.EJ. Newman. Power laws, pareto distributions and zipf's law. Contemporary Physics, 46:323Ð351, 2005.
#
# dpowerlaw(), ppowerlaw(), qpowerlaw(), and rpowerlaw():
#  PDF, CDF, Quantile, and RNG functions for the continuous powerlaw distribution
#
# ddpowerlaw(), pdpowerlaw(), qdpowerlaw(), and rdpowerlaw():
#  PDF, CDF, Quantile, and RNG functions for the discrete powerlaw distribution
#
# These functions all have the same parameters as the corresponding functions for other
# distributions in R, except for the RNG:
#
# rpowerlaw(n, alpha, xmin=1)
# - Generate n random numbers from a continuous powerlaw distribution.
# 
# rdpowerlaw(n, alpha, xmin=1)
# - Generate n random numbers from a discrete powerlaw distribution
#
# Model fit data to a powerlaw distribution:
#
# m <- plm(data)
# - Fit a powerlaw distribution to data.   data can be either continuous
# or discrete.
# 
# m <- plm(data, xmin=1)
# - Fit a powerlaw distribution to data, with xmin pre-specified instead
# of being estimated.   Much faster.
# 
# print(m)
# m
# - Print the estimated powerlaw distribution information.
# 
# plot(m, ...)
# - Plot a graph of the original data and the fitted powerlaw
# distribution.   ... are passed as options to xyplot.
# 
# summary(m)
# - Run a Kolomogorov-Smirnov goodness-of-fit test for the powerlaw
# distribution.   Low p values indicate that you can reject the powerlaw
# distribution as a possible distribution.  Somewhat slow, as the
# p-value requires monte-carlo simulations to estimate.
# 
# summary(m, epsilon=0.025)
# - Run a KS goodness-of-fit test.   Estimate the p-value to
# approximately epsilon accuracy.   Note that running time of
# calculating the p-value is O(epsilon^2) (meaning, don't estimate it
# all that accurately.   Default is epsilon=0.1


library(gsl)
library(numDeriv)

# PDF, CDF, Quantile, and RNG functions for continuous power law distribution
dpowerlaw <- function(x, alpha=2, xmin=1, log=F) {
  if (log)
    log(alpha-1) - log(xmin) - alpha * log(x / xmin)
  else
    ((alpha - 1) / xmin) * ((x / xmin) ^ (-alpha))
}

ppowerlaw <- function(q, alpha=2, xmin=1, lower.tail=T, log.p = F) {
  p <- (q / xmin) ^ (- alpha + 1)
  if (lower.tail)
    p <- 1-p
  if (log.p)
    p <- log(p)
  p
}

qpowerlaw <- function(p, alpha=2, xmin=1, lower.tail=T, log.p = F) {
  if (!lower.tail)
    p <- 1-p
  if (log.p)
    p <- exp(p)
  xmin * ((1 - p) ^ (-1 / (alpha - 1)))
}  

rpowerlaw <- function(n, alpha=2, xmin=1) {
  qpowerlaw(runif(n, 0, 1), alpha, xmin)
}


ddpowerlaw <- function(x, alpha = 2, xmin=1, log=F) {
  if (log)
    -log(hzeta(alpha, xmin)) - alpha * log(x)
  else
    x ^ (-alpha) / hzeta(alpha, xmin)
}


# PDF, CDF, Quantile, and RNG functions for discrete power law distributions
pdpowerlaw <- function(q, alpha=2, xmin=1, lower.tail=T, log.p=F) {
  p <- hzeta(alpha, q) / hzeta(alpha, xmin)
  if (lower.tail)
    p <- 1-p
  if (log.p)
    p <- log(p)
  p
}

  # Helper function
individual_qdpowerlaw <- function(r, alpha, xmin) {
  x2 <- xmin
  # Doubling Up
  repeat {
    x1 <- x2
    x2 <- 2*x1
    if (pdpowerlaw(x2, alpha, xmin, lower.tail=F) < 1-r) break
  }
  # Binary Search
  repeat {
    temp <- (x2 + x1) / 2
    if (pdpowerlaw(temp, alpha, xmin, lower.tail=F) < 1-r) {
      x2 <- temp
    } else {
      x1 <- temp
    }
    if (ceiling(x2) - floor(x1) < 2) break
  }
  floor(x1)
}

qdpowerlaw <- function(p, alpha=2, xmin=1, approx=F, lower.tail=T, log.p=F) {
  if (approx) {
    return(floor(qpowerlaw(p, alpha, xmin - 0.5, lower.tail, log.p) + 0.5))
  }
  if (log.p)
    p <- exp(p)
  if (!lower.tail)
    p <- 1-p
  sapply(p, function(x) { individual_qdpowerlaw(x, alpha, xmin) })
}

rdpowerlaw <- function(n, alpha=2, xmin=1, approx=F) {
  #if (approx) {
    #return(floor(rpowerlaw(n, alpha, xmin - 0.5) + 0.5))
  #}
  r <- runif(n, 0, 1)
  qdpowerlaw(r, alpha, xmin, approx, lower.tail=FALSE, log.p=FALSE)
}


# Numerically fit a discrete powerlaw distribution via maximum likelihood
plm.fit.discrete.numeric <- function(data, xmin=1, alpha.starting = 2) {
  # nlm the negative log likelihood: n log zeta(a, xmin) + a sum log x_i
  N <- length(data)

  #dpowerlaw.nloglik <- function(alpha) { N * log(hzeta(alpha, xmin)) + alpha * sum(log(data)) }
  dpowerlaw.nloglik <- function(alpha) { -sum(ddpowerlaw(data, alpha, xmin, log=T)) }

  alpha.info <- nlm(dpowerlaw.nloglik, alpha.starting)
  alpha <- alpha.info$estimate

  x2 <- grad(hzeta, alpha, q=xmin) / hzeta(alpha, xmin)
  x1 <- hessian(hzeta, alpha, q=xmin) / hzeta(alpha, xmin)
  sigma <- 1 / (N * (x1 - x2^2))
  list(alpha = alpha, sigma = sigma, N = N, xmin=xmin, xmin.estimated=F)
}

# Fit a discrete power law distribution using a continuous approximation
plm.fit.discrete.approx <- function(data, xmin = 1, ...) {
  N <- length(data)
  x <- sum(log(data / (xmin - 0.5)))
  alpha <- 1 + N / x

  sigma <- (alpha - 1) / sqrt(N)

  list(alpha = alpha, sigma = sigma, N = N, xmin=xmin, xmin.estimated=F)
}

plm.fit.discrete <- function(data, approx=F, xmin, ...) {
  data <- data[data >= xmin]
  if (approx)
    plm.fit.discrete.approx(data, xmin, ...)
  else
    plm.fit.discrete.numeric(data, xmin, ...)
}

# Fit a continuous power law distribution via maximum likelihood (direct calculation)
plm.fit.continuous <- function(data, xmin = 1) {
  data <- data[data >= xmin]
  N <- length(data)
  x <- sum(log(data / xmin))
  alpha <- 1 + N / x
  sigma <- (alpha - 1) / sqrt(N)

  list(alpha = alpha, sigma = sigma, N = N, xmin=xmin, xmin.estimated=F)
}

# Helper function to calculate KS statistic
ks.plm <- function(data.cdf, x, alpha, xmin, discrete=T) {
  data.cdf <- data.cdf[x >= xmin]
  x <- x[x >= xmin]
  fit.cdf <- if (discrete)
               pdpowerlaw(x, alpha, xmin, lower.tail=F)
             else
               ppowerlaw(x, alpha, xmin, lower.tail=F)
  KS <- max(abs(data.cdf - fit.cdf))
}

# Use search to simultaneously estimate xmin and alpha (discrete)
plm.fit.xmin.discrete <- function(data, approx = F, xmin.starting = 1, alpha.starting = 2) {
  N <- length(data)
  x <- sort(unique(data))    
  cdf <- function(n) { sum(data >= n) / N }
  cdf.data <- sapply(x, cdf)

  try_xmin <- function(xmin) {
    if (xmin >= 1) {
      fit <- plm.fit.discrete(data, approx, xmin=xmin, alpha.starting=alpha.starting)
      KS <- ks.plm(cdf.data, x, fit$alpha, xmin, discrete=T)
    } else {
      fit <- NULL
      KS <- 2 # A high number so it never gets chosen
    }
    list(xmin = xmin, fit=fit, KS=KS)
  }
  
  current <- try_xmin(xmin.starting)
  left <- try_xmin(current$xmin - 1)
  right <- try_xmin(current$xmin + 1)

  repeat {
    # Test if we found the optimal xmin
    if ((current$KS <= left$KS) & (current$KS <= right$KS))
      break
    # Test if we have gone too far
    if (right$xmin > max(data))
      break
    
    # Move left or right appropriately
    if ((current$KS - left$KS) > (current$KS - right$KS)) {
      # Move left one
      right <- current
      current <- left
      left <- try_xmin(current$xmin - 1)
    } else {
      # Move right one
      left <- current
      current <- right
      right <- try_xmin(current$xmin + 1)
    }
  }
  
  list(alpha=current$fit$alpha, sigma=current$fit$sigma, N=current$fit$N, xmin=current$xmin, xmin.estimated=T, KS=current$KS)
}

# Use search to simultaneously estimate xmin and alpha (continuous)
plm.fit.xmin.continuous <- function(data, xmin.starting) {
  N <- length(data)
  x <- sort(unique(data))    
  cdf <- function(n) { sum(data >= n) / N }
  cdf.data <- sapply(x, cdf)

  try_xmin <- function(xmin) {
    fit <- plm.fit.continuous(data, xmin)
    KS <- ks.plm(cdf.data, x, fit$alpha, xmin, discrete=F)
    KS
  }

  out <- nlm(try_xmin, xmin.starting)
  fit <- plm.fit.continuous(data, out$estimate)
  
  list(alpha = fit$alpha, sigma=fit$sigma, N=fit$N, xmin=fit$xmin, xmin.estimated=T, KS=out$minimum)
}

# Generic function for a power law model -- Fit data to a power-law distribution
plm <- function(data, ...) {
  UseMethod("plm")
}

plm.default <- function(data, xmin=NULL, discrete = NULL, approx=FALSE, alpha.starting = 2, xmin.starting=1) {
  if (is.null(discrete)) {
    discrete = all(data == as.integer(data))
  }
  if (is.null(xmin)) {
    # Estimate xmin as well as alpha
    info <- if (discrete)
              plm.fit.xmin.discrete(data, approx, xmin.starting=xmin.starting, alpha.starting=alpha.starting)
            else
              plm.fit.xmin.continuous(data, xmin.starting=xmin.starting)
  } else {
    info <- if (discrete)
              plm.fit.discrete(data, approx, xmin, alpha.starting=alpha.starting)
            else 
              plm.fit.continuous(data, xmin)
  }
  info$data <- data
  info$discrete <- discrete
  info$approx <- approx
  class(info) <- "plm"
  info
}

plm.factor <- function(data, ...) {
  data <- tapply(data, data, length)
  NextMethod("plm")
#  plm.default(counts, discrete=T, ...)
}

# Helper functions for the "plm" class
print.plm <- function(x) {
  if (x$discrete) {
    cat("Discrete Power Law Distribution Fit\n")
  } else {
    cat("Power Law Distribution Fit\n")   
  }
  cat(paste("  Exponent:", x$alpha, "+-", x$sigma, "\n"))
  cat(paste("  Minimum X:", x$xmin, "\n"))
  if (x$xmin.estimated)
    cat(paste("    (estimated: KS =", x$KS, ")\n"))
  cat(paste("  N:", x$N, "\n"))
}

# Not sure if this works
# Usage:
#   lrtest.plm(m, dexp, 0.125)        To test against exponential
#   lrtest.plm(m, dlnorm, c(0.5, 2))  To test against log-normal
lrtest.plm <- function(x, dist_func, start = list()) {
  if (x$discrete)
    loglik.powerlaw <- function(alpha) { ddpowerlaw(x$data, alpha, x$xmin, log=T) }
  else
    loglik.powerlaw <- function(alpha) { dpowerlaw(x$data, alpha, x$xmin, log=T) }
  loglik.other <- function(y) { do.call(dist_func, c(list(x=x$data, log=T), y)) }

  nloglik.other <- function(y) { -sum(loglik.other(y)) }

  other.estimate <- nlm(nloglik.other, start)

  #R <- (sum(loglik.powerlaw(x$alpha)) - sum(loglik.other(other.estimate$estimate))) / length(x$data)
  R <- sum(loglik.powerlaw(x$alpha)) - sum(loglik.other(other.estimate$estimate))
  pl.mean <- mean(loglik.powerlaw(x$alpha))
  other.mean <- mean(loglik.other(other.estimate$estimate))
  sigsq <- sum(((loglik.powerlaw(x$alpha) - loglik.other(other.estimate$estimate)) - (pl.mean - other.mean))^2) / length(x$data)

  # Copied from help(Normal)
  erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
  
  p <- erfc(R / sqrt(2 * length(x$data) * sigsq))
  l <- sum(loglik.powerlaw(x$alpha))
  list(R = R, sigsq = sigsq, p=p, fit=other.estimate, loglik=l)
}

summary.plm <- function(x, num_datasets = NULL, epsilon=0.1) {
  # Calculate goodness of fit statistic and p-value
  N <- length(x$data)
  if (x$xmin.estimated)
    KS.real <- x$KS
  else {
    y <- sort(unique(x$data))    
    cdf <- function(n) { sum(x$data >= n) / N }
    cdf.data <- sapply(y, cdf)
    KS.real <- ks.plm(cdf.data, y, x$alpha, x$xmin, x$discrete)
  }

  # Set number of datasets based on intended accuracy of p-value
  if (is.null(num_datasets))
    num_datasets <- 0.25 * epsilon ^ (-2)

  # Probability of being below xmin
  data.low <- x$data[x$data < x$xmin]
  p <- length(data.low) / length(x$data)
  run_dataset <- function(N) {
    # Generate data
    if (length(data.low) > 0) {
      choices <- rbinom(N, 1, p)
      N.info <- tapply(rep(1,N), choices, length)
      N.low <- N.info[2]
      N.high <- N.info[1]
      rd.low <- data.low[floor(runif(N.low, 1, length(data.low)+1))]
    } else {
      N.high <- N
      rd.low <- numeric(0)
    }
    rd.high <- if (x$discrete)
                 rdpowerlaw(N.high, x$alpha, x$xmin)
               else
                 rpowerlaw(N.high, x$alpha, x$xmin)
    # Estimate model and KS statistic
    rd <- c(rd.low, rd.high)
    if (x$xmin.estimated) {
      fit <- plm(rd, xmin=NULL)
      KS <- fit$KS
    } else {
      fit <- plm(rd, xmin=x$xmin, discrete=x$discrete)
      y <- sort(unique(rd))    
      cdf <- function(n) { sum(rd >= n) / N }
      cdf.data <- sapply(y, cdf)
      KS <- ks.plm(cdf.data, y, x$alpha, x$xmin, x$discrete)
    }
    KS  
  }

  # Get num_datasets KS statistics
  KS.sim <- sapply(rep(N, num_datasets), run_dataset)

  p <- sum(KS.sim > KS.real) / length(KS.sim)
  
  out <- list(fit = x, KS = KS.real, p = p, sim.results = KS.sim)
  class(out) <- c("summary.plm", class(out))
  out
}

print.summary.plm <- function(x) {
  if (x$fit$discrete) {
    cat("Discrete Power Law Fit")
  } else {
    cat("Power Law Fit")   
  }
  cat(paste(" with alpha =", x$fit$alpha, "+-", x$fit$sigma, "\n"))
  cat(paste("\nKolmogorov-Smirnov Goodness-of-Fit test\n  Null Hypothesis is that the data fits a power-law\n"))
  cat(paste("  KS =", x$KS, "\n    (p-value:", x$p, ")\n"));
}

# Plots the cumulative distribution function for both the data and the fit
plot.plm <- function(x, hypo.col = "red", xlab="Value", ylab="Probability X > x", scales=list(x=list(log=T), y=list(log=T)), ...) {
  N <- length(x$data)
  plot_x <- sort(unique(x$data))    
  cdf <- function(n) { sum(x$data >= n) / N }
  real_y <- sapply(plot_x, cdf)
  if (x$discrete) {
    hypo_y <- pdpowerlaw(plot_x, x$alpha, x$xmin, lower.tail=F)
  } else {
    hypo_y <- ppowerlaw(plot_x, x$alpha, x$xmin, lower.tail=F)
  }
  df <- data.frame(X = plot_x, Data = real_y, Fit = hypo_y)
  df <- melt(df, measure.var = c("Data", "Fit"))
  xyplot(value ~ X, data=df, groups=variable, scales=scales, ylab = ylab, xlab=xlab, ...)
}





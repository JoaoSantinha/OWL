## Based on VLAD NICULAE's python implementation: https://github.com/vene/pyowl
library(Matrix)
library(Iso)
library(phonTools)

fista <- function(sfunc, nsfunc, x0, max.iter=500, max.linesearch=20, eta=2.0, tol=1e-3, verbose=0) {
  y <- x0
  x <- y
  L <- 1.0
  t <- 1.0

  for (iter in seq.int(1,max.iter)) {
    output <- sfunc(y, TRUE)
    f_old <- output[1][[1]][[1]]
    grad <- output[2][[1]]

    for (ls in seq.int(max.linesearch)) {
      y_proj <- nsfunc(y - grad / L, L)

      diff   <- as.vector(t(y_proj - y)) # python is collumn major and r is row major so to ravel we need to transpose
      sqdist <- diff %*% diff
      dist   <- sqrt(sqdist)

      F <- sfunc(y_proj)
      Q <- f_old + diff %*% as.vector(t(grad)) + 0.5 * L * sqdist

      if (F <= Q) {
        break
      }

      L <- L * eta
    }
    if (ls == max.linesearch - 1 & verbose) {
      print("Line search did not converge.")
    }

    if (verbose) {
      sprint("%d. %f", (it + 1), dist)
    }

    if (dist <= tol) {
      if (verbose) {
        print("Converged.")
      }
      break
    }

    x_next <- y_proj
    t_next <- (1 + sqrt(1 + 4 * t ** 2)) / 2.
    y <- x_next + (t-1) / t_next * (x_next - x)
    t <- t_next
    x <- x_next
  }
  return (y_proj)
}

fit_owl_fista <- function(X, y, w, loss, max.iter=500, max.linesearch=20, eta=2.0, tol=1e-3, verbose=FALSE) {

  #least squares loss
  sfunc <- function(coef, grad=FALSE) {

    y_scores = as.matrix(X) %*% coef #safe_sparse_dot(X, coef)

    if (grad) {

      output_loss <- loss(y, y_scores, TRUE)
      obj <- output_loss[1]
      lp <- output_loss[2][[1]]

      grad <- t(X) %*% lp

      return (list(obj, grad))
    } else {
      return (loss(y, y_scores, FALSE))
    }
  }

  nsfunc <- function(coef, L) {
    return (prox_owl(coef, w / L))
  }

  coef <- rep(0,dim(X)[2])

  return(fista(sfunc, nsfunc, coef, max.iter, max.linesearch, eta, tol, verbose))
}


prox_owl <- function(v, w) {

  v_abs <- abs(v)

  sorting <- sort(v_abs, decreasing = TRUE, index.return= TRUE)
  ix <- sorting$ix

  v_abs <- v_abs[ix]

  v_abs <- pava(v_abs - w, decreasing = TRUE)
  v_abs[v_abs < 0] <- 0

  # undo the sorting
  inv.ix <- phonTools::zeros(ix)
  inv.ix[ix] <- seq(length(v))
  v_abs <- v_abs[inv.ix]

  return (sign(v) * v_abs)
}

oscar_weights <- function(lambda1, lambda2, size) {
  w <- seq(size - 1, 0, -1)
  w <- w * lambda2
  w <- w + lambda1
  return (w)
}

squared_loss <-function(y_true=y_true, y_pred=y_pred, return_derivative=FALSE) {
  diff_val <- y_pred - y_true
  obj <- 0.5 * (t(diff_val) %*% diff_val)
  if (return_derivative) {
    return (list(obj, diff_val))
  } else {
    return (obj)
  }

}

squared_hinge_loss <-function(y_true=y_true, y_scores=y_scores, return_derivative=FALSE) {
  #labels in (-1, 1)
  z <- pmax(0, 1- y_true * y_scores)
  obj <- sum(z^2)

  if (return_derivative) {
    return (list(obj, -2 * y_true * z))
  } else {
    return (obj)
  }
}

OWL <- function(x,
                y,
                family= c("gaussian", "binomial", "multinomial"#, "coefficients"
                          ),
                intercept = TRUE,
                max.iter=500,
                max.linesearch=20,
                eta=2.0,
                tol=1e-3,
                lambda1=1,
                lambda2=100) {

  ocall <- match.call()

  family <- match.arg(family)

  if (!is.null(levels(y))) {
    class_names <- levels(y)
  } else {
    class_names <- NA_character_
  }

  if (is.factor(y)) {
    y <- as.numeric(y) -1
  }

  n <- NROW(x)
  p <- NCOL(x)

  if (family == "gaussian") {
    loss <- squared_loss
  } else {
    loss <- squared_hinge_loss
  }

  weights_vec <- oscar_weights(lambda1, lambda2, p)

  coefficients <- fit_owl_fista(x, y, weights_vec, loss = loss, max.iter = max.iter,
                                max.linesearch = max.linesearch, eta = eta, tol = tol,
                                verbose = FALSE)

  nonzeros <- apply(coefficients, c(1), function(x) abs(x) > 0)

  owl_class <- switch(family,
                      gaussian = "GaussianOWL",
                      binomial = "BinomialOWL",
                      multinomial = "MultinomialOWL")

  structure(list(coefficients = coefficients,
                 nonzeros = nonzeros,
                 lambda1 = lambda1,
                 lambda2 = lambda2,
                 class_names = class_names,
                 family = family,
                 call = ocall),
            class = c(owl_class, "OWL"))
}

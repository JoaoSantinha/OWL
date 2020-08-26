caretOWL <- function() {
  list(
    label = "OWL",
    library = c("Matrix"),
    type = c("Regression", "Classification"),

    parameters = data.frame(parameter = c("lambda1", "lambda2"),
                            class = c("numeric", "numeric"),
                            label = c("lambda 1", "lambda 2")),

    grid = function(x, y, len = NULL, search = "grid") {

      lambda1 <- 2^stats::runif(len, min = -6, max = 3)
      lambda2 <- 2^stats::runif(len, min = -6, max = 3)

      if (search == "grid") {
        out <- expand.grid(lambda1 = lambda1,
                           lambda2 = lambda2)
      } else {


        out <- data.frame(lambda1 = lambda1,
                          lambda2 = lambda2)
      }

      out
    },

    loop = function(grid) {

      lmbd1 <- unique(grid$lambda1)
      loop <- data.frame(lambda1 = lmbd1)
      loop$lambda2 <- NA

      submodels <- vector("list", length = length(lmbd1))

      for (i in seq_along(lmbd1)) {
        np <- grid[grid$lambda1 == lmbd1[i], "lambda2"]
        loop$lambda2[loop$lambda1 == lmbd1[i]] <- np[which.max(np)]
        submodels[[i]] <- data.frame(lambda2 = np[-which.max(np)])
      }

      list(loop = loop, submodels = submodels)
    },

    fit = function(x, y, wts=NULL, param, lev, last, classProbs, ...) {

      dots <- list(...)

      numLev <- if (is.character(y) | is.factor(y)) {
        if (is.character(y)) {
          y <- as.factor(y)
        }
        length(levels(y))
      } else
        NA

      if (all(names(dots) != "family")) {
        if (!is.na(numLev)) {
          fam <- ifelse(numLev > 2, "multinomial", "binomial")
        } else {
          fam <- "gaussian"
        }

        dots$family <- fam
      }

      if(!is.null(wts)) dots$weights <- wts

      if(!(class(x)[1] %in% c("matrix", "sparseMatrix")))
        x <- Matrix::as.matrix(x)

      modelArgs <- c(list(x = x,
                          y = y,
                          lambda1 = param$lambda1),
                     dots)

      out <- do.call(OWL, modelArgs)
      if(!is.na(param$lambda2[1])) out$lambda2Opt <- param$lambda2[1]
      out
    },

    predict = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
      if (!is.matrix(newdata))
        newdata <- Matrix::as.matrix(newdata)
      obsLevels <- if("class_names" %in% names(modelFit))
        modelFit$class_names
      print(obsLevels)
      print(modelFit$class_names)
      if (length(obsLevels) < 2) {
        out <- predict.GaussianOWL(modelFit,
                                   newdata,
                                   lambda1 = modelFit$lambda1Opt,
                                   lambda2 = modelFit$lambda2Opt)
      } else {
        if (length(obsLevels) == 2) {
          out <- predict.BinomialOWL(modelFit,
                                     newdata,
                                     lambda1 = modelFit$lambda1Opt,
                                     lambda2 = modelFit$lambda2Opt,
                                     type = "class")
        } else {
          out <- predict.MultinomialOWL(modelFit,
                                        newdata,
                                        lambda1 = modelFit$lambda1Opt,
                                        lambda2 = modelFit$lambda2Opt,
                                        type = "class")
        }
      }

      if (is.matrix(out))
        out <- out[, 1]

      if (!is.null(submodels)) {
        if (length(modelFit$obsLevels) < 2) {
          tmp <- predict(modelFit,
                         newdata,
                         lambda1 = submodels$lambda1,
                         lambda2 = submodels$lambda2)
          tmp <- as.list(as.data.frame(tmp))
        } else {
          tmp <- predict(modelFit,
                         newdata,
                         lambda1 = submodels$lambda1,
                         lambda2 = submodels$lambda2,
                         type = "class")
          tmp <- if (is.matrix(tmp))
            as.data.frame(tmp, stringsAsFactors = FALSE)
          else
            as.character(tmp)
          tmp <- as.list(tmp)
        }
        out <- c(list(out), tmp)
      }
      out
    },

    prob = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
      obsLevels <- if("class_names" %in% names(modelFit))
        modelFit$class_names
      else
        NULL

      probs <- predict(modelFit,
                       Matrix::as.matrix(newdata),
                       lambda1 = modelFit$lambda1Opt,
                       lambda2 = modelFit$lambda2Opt,
                       type = "response")

      if (length(obsLevels) == 2) {
        probs <- as.vector(probs)
        probs <- as.data.frame(cbind(1-probs, probs), stringsAsFactors = FALSE)
        colnames(probs) <- modelFit$obsLevels
      } else {
        probs <- as.data.frame(probs[,,1,drop = FALSE], stringsAsFactors = FALSE)
        names(probs) <- modelFit$obsLevels
      }

      if (!is.null(submodels)) {
        tmp <- predict(modelFit,
                       Matrix::as.matrix(newdata),
                       lambda1 = submodels$lambda1Opt,
                       lambda2 = submodels$lambda2Opt,
                       type = "response")

        if(length(obsLevels) == 2) {
          tmp <- as.list(as.data.frame(tmp, stringsAsFactors = TRUE))
          tmp <- lapply(tmp,
                        function(x, lev) {
                          x <- as.vector(x)
                          tmp <- data.frame(1-x, x)
                          names(tmp) <- lev
                          tmp
                        },
                        lev = obsLevels)
        } else tmp <- apply(tmp, 3, function(x) data.frame(x))
        probs <- if (is.list(tmp)) c(list(probs), tmp) else list(probs, tmp)
      }
      probs
    },

    predictors = function(x, lambda1 = NULL, lambda2 = NULL, ...) {
      if (is.null(lambda1) || is.null(lambda2)) {
        if (length(lambda1) > 1 || length(lambda2) > 1)
          stop("Only one value of lambda1 and lambda2 are allowed")

        if (!is.null(x$lambda1Opt)) {
          lambda1 <- x$lambda1Opt
        } else  {
          stop("must supply a value of lambda1")
        }
        if (!is.null(x$lambda2Opt)) {
          lambda2 <- x$lambd2Opt
        } else  {
          stop("must supply a value of lambda2")
        }
      }

      allVar <- rownames(stats::coef(x, simplify = FALSE))

      out <- apply(abs(allVar) > 0, 3, sum)

      out <- unique(out)
      if (length(out) > 0) {
        out <- out[!is.na(out)]
        out <- allVar[out]
      }

      out
    },

    varImp = function(object, lambda1 = NULL, lambda2 = NULL, ...) {
      if (is.null(lambda1) || is.null(lambda2)) {
        if (length(lambda1) > 1 || length(lambda2) > 1)
          stop("Only one value of lambda1 and lambda2 are allowed")

        if (!is.null(object$lambda2Opt) && !is.null(object$lambda2Opt)) {
          lambda1 <- object$lambda1Opt
          lambda2 <- object$lambda2Opt
        } else {
          stop("must supply a value of lambda1  and lambda2")
        }
      }

      beta <- stats::coef(object, lambda1 = lambda1, lambda2 = lambda2, simplify = TRUE)
      beta <- as.data.frame(beta)
      out <- as.data.frame(Overall = beta[, 1])
      out <- abs(out[rownames(out) != "(Intercept)", , drop = FALSE])
      out
    },

    levels = function(x) {
      if (any(names(x) == "obsLevels"))
        x$obsLevels
      else
        NULL
    },

    sort = function(x) {
      x[order(x$lambda1, x$lambda2),]
    },

    trim = function(x) {
      x$call <- NULL
      x
    },

    tags = c("Generalized Linear Model",
             "Implicit Feature Selection",
             "L1 Regularization",
             "Linear Classifier",
             "Linear Regression")
  )
}

#' @rdname predict.OWL
#' @export
predict.OWL <- function(object,
                        x,
                        lambda1 = NULL,
                        lambda2 = NULL,
                        type=c("link","response","coefficients","class"),
                        simplify = TRUE,
                        ...) {

  type=match.arg(type)
  # This method (the base method) only generates linear predictors
  if (inherits(x, "sparseMatrix"))
    x <- methods::as(x, "dgCMatrix")

  if (inherits(x, "data.frame"))
    x <- as.matrix(x)

  beta <- object$coefficients

  if (type == "coefficients") {
    beta <- stats::coef(object, lambda1 = lambda1, lambda2 = lambda2, simplify = FALSE)
    return(beta)
  } else {

    n <- NROW(x)
    p <- NROW(beta)
    m <- NCOL(beta)
    n_penalties <- dim(beta)[3]

    print(NROW(beta))
    print(NCOL(x))
    stopifnot(p == NCOL(x))

    if (is.na(dim(beta)[3])) {
      lin_pred <- array(dim = c(n, m),
                        dimnames = list(rownames(x),
                                        dimnames(beta)[[2]]))

      lin_pred <- as.matrix(x %*% beta)
    } else {
      lin_pred <- array(dim = c(n, m, n_penalties),
                        dimnames = list(rownames(x),
                                        dimnames(beta)[[2]],
                                        dimnames(beta)[[3]]))

      for (i in seq_len(n_penalties))
        lin_pred[, , i] <- as.matrix(x %*% beta[, , i])
    }
    return(lin_pred)
  }
}

#' @rdname predict.OWL
#' @export
predict.GaussianOWL <- function(object,
                                x,
                                lambda1 = NULL,
                                lambda2 = NULL,
                                type = c("link", "response", "coefficients"),
                                simplify = TRUE,
                                ...) {
  type <- match.arg(type)

  out <- NextMethod(object, x, lambda1, lambda2, type = type) # always linear predictors

  if(type != "coefficients") {
    if (simplify)
      out <- drop(out)
  }
  out
}

#' @rdname predict.OWL
#' @export
predict.BinomialOWL <- function(object,
                                x,
                                lambda1 = NULL,
                                lambda2 = NULL,
                                type = c("link", "response", "class", "coefficients"),
                                simplify = TRUE,
                                ...) {

  type <- match.arg(type)

  nbeta <- object$coefficients
  lin_pred=as.matrix(x)%*%nbeta

  if (type != "coefficients") {
    out <- switch(
      type,
      link = lin_pred,
      response = 1 / (1 + exp(-lin_pred)),
      class = {
        cnum <- ifelse(lin_pred > 0, 2, 1)
        clet <- object$class_names[cnum]

        if (is.matrix(cnum))
          clet <- array(clet, dim(cnum), dimnames(cnum))

        clet
      }
    )

    if (simplify)
      out <- drop(out)
  } else {
    out <- lin_pred
  }
  out
}

#' @export
#' @rdname predict.OWL
predict.MultinomialOWL <- function(object,
                                   x,
                                   lambda1 = NULL,
                                   lambda2 = NULL,
                                   type = c("link", "response", "class", "coefficients"),
                                   exact = FALSE,
                                   simplify = TRUE,
                                   ...) {
  type <- match.arg(type)

  lin_pred <- NextMethod(object, type = type, simplify = FALSE)
  m <- NCOL(lin_pred)

  if (type != "coefficients") {
    out <- switch(
      type,

      response = {
        n <- nrow(lin_pred)
        m <- ncol(lin_pred)
        path_length <- dim(lin_pred)[3]

        tmp <- array(0, c(n, m + 1, path_length))
        tmp[, 1:m, ] <- lin_pred

        aperm(apply(tmp, c(1, 3), function(x) exp(x)/sum(exp(x))), c(2, 1, 3))
      },

      link = lin_pred,

      class = {
        response <- stats::predict(object, x, type = "response", simplify = FALSE)
        tmp <- apply(response, c(1, 3), which.max)
        class_names <- object$class_names

        predicted_classes <-
          apply(tmp,
                2,
                function(a) factor(a, levels = 1:(m+1), labels = class_names))
        colnames(predicted_classes) <- dimnames(lin_pred)[[3]]
        predicted_classes
      }
    )
    if (simplify)
      out <- drop(out)
  } else {
    out <- lin_pred
  }

  out
}

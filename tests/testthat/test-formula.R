

df <- data.frame(
  y      = rnorm(60),
  x_cont = rnorm(60),
  x_bin  = rbinom(60, 1, .5),
  x_fac3 = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C")),
  x_fac4 = factor(rep(c("A", "B", "C", "D"), 15), levels = c("A", "B", "C", "D"))
)


library(stanova)
fit_stanova <- stanova::stanova_lm(y ~ x_cont + x_fac3, data = df)
fit_stanova <- stanova::stanova_lm(y ~ x_fac3 * x_fac4, data = df)

debugonce(stanova::stanova)
debugonce(rstanarm::stan_glm)


summary(fit_stanova)
fit_stanova$coefficients

JAGS_formula <- function(formula, parameter, data, prior_list){

  if(!is.language(formula))
    stop("'formula' must be a formula")
  check_char("parameter", parameter)
  if(!is.data.frame(data))
    stop("'data' must be a data.frame")
  check_list(prior_list, "prior_list")
  if(any(!sapply(prior_list, is.prior)))
    stop("'prior_list' must be a list of priors.")


  # remove the specified response (would crash the model.frame if not included)
  formula <- .remove_response(formula)

  # obtain predictors characteristics factors
  formula_terms    <- stats::terms(formula)
  has_intercept    <- attr(formula_terms, "intercept") == 1
  model_terms      <- c(if(has_intercept) "intercept", attr(formula_terms, "term.labels"))
  model_terms_type <- sapply(model_terms, function(model_terms){
    if(model_terms == "intercept")
      return("intercept")
    if(grepl(":", model_terms)){
      return("interaction")
    }else{
      return("simple")
    }
  })
  predictors       <- as.character(attr(formula_terms, "variables"))[-1]
  predictors_type  <- sapply(predictors, function(predictor){
    if(is.null(data[,predictor]))
      stop("The predictor variable is missing in the data set.")
    if(is.factor(data[,predictor]) | is.character(data[,predictor])){
      return("factor")
    }else{
      return("continuous")
    }
  })

  # check that all predictors have a prior distribution
  check_list(prior_list, "prior_list", name = model_terms, allow_other = FALSE, all_objects = TRUE)

  # assign factor contrasts to the data based on prior distributions
  if(any(predictors_type == "factor")){
    for(factor in names(predictors_type[predictors_type == "factor"])){
      if(is.prior.dummy(prior_list[[factor]])){
        stats::contrasts(data[,factor]) <- "contr.treatment"
      }else if(is.prior.orthonormal(prior_list[[factor]])){
        stats::contrasts(data[,factor]) <- "contr.orthonormal"
      }else{
        stop("Unsupported prior distribution defined for a factor variable. See '?prior_factor' for details.")
      }
    }
  }

  # get the default design matrix
  model_frame  <- stats::model.frame(formula, data = data)
  model_matrix <- stats::model.matrix(model_frame, formula = formula, data = data)

  # prepare syntax & data based on the formula
  formula_syntax <- NULL
  prior_syntax   <- NULL
  JAGS_data      <- list()
  JAGS_data[[paste0("N_", parameter)]] <- nrow(data)

  # add intercept and prepare the indexing vector
  if(has_intercept){
    terms_indexes    <- attr(model_matrix, "assign") + 1
    terms_indexes[1] <- NA

    formula_syntax   <- c(formula_syntax, paste0(parameter, "_intercept"))
    prior_syntax     <- c(prior_syntax, .JAGS_prior.simple(prior_list[["intercept"]], paste0(parameter, "_intercept")))
  }else{
    terms_indexes    <- attr(model_matrix, "assign")
  }

  # add remaining terms (omitting the intercept indexed as NA)
  for(i in unique(na.omit(terms_indexes))){

    if(all(predictors_type[predictors == model_terms[i]] == "continuous")){
      # continuous variables or interactions of continuous variables are simple predictors

      JAGS_data[[paste0(parameter, "_data_", model_terms[i])]] <- model_matrix[,terms_indexes == i]
      formula_syntax <- c(formula_syntax, paste0(parameter, "_", model_terms[i]), " * ", paste0(parameter, "_data_", model_terms[i], "[i]"))
      prior_syntax   <- c(prior_syntax, .JAGS_prior.simple(prior_list[["intercept"]], paste0(parameter, "_", model_terms[i])))

    }else if(any(predictors_type[predictors == model_terms[i]] == "factor")){
      # factor variables or interactions with a factor requires factor style prior

      JAGS_data[[paste0(parameter, "_data_", model_terms[i])]] <- model_matrix[,terms_indexes == i, drop = FALSE]
      formula_syntax <- c(formula_syntax, "inprod(", paste0(parameter, "_", model_terms[i]), ", ", paste0(parameter, "_data_", model_terms[i], "[i,]"), ")")
      prior_syntax   <- c(prior_syntax, .JAGS_prior.factor(prior_list[["intercept"]], paste0(parameter, "_", model_terms[i]), levels = sum(terms_indexes == i) + 1))

    }else{
      stop("Unrecognized model term.")
    }

  }


  # finish the syntax
  formula_syntax <- paste0(
    "for(i in 1:N_", parameter, "){\n",
    parameter, "[i] = ", paste0(formula_syntax, collapse = " + "), "\n",
    "}")
  prior_syntax   <- paste0(prior_syntax, collapse = " ")


  return(list(
    formula_syntax = formula_syntax,
    priors_syntax  = "",
    data           = JAGS_data,
    inits          = "",
    parameters     = "",
  ))
}

.remove_response  <- function(formula){
  if(attr(stats::terms(formula), "response")  == 1){
    formula[2] <- NULL
  }
  return(formula)
}
contr.orthonormal <- function(n, contrasts = TRUE){
  # based on: stanova::contr.bayes
  if(length(n) <= 1L){
    if(is.numeric(n) && length(n) == 1L && n > 1L){
      return(TRUE)
    }else{
      stop("Not enough degrees of freedom to define contrasts.")
    }
  }else{
    n <- length(n)
  }

  cont <- diag(n)
  if(contrasts){
    a       <- n
    I_a     <- diag(a)
    J_a     <- matrix(1, nrow = a, ncol = a)
    Sigma_a <- I_a - J_a/a
    cont    <- eigen(Sigma_a)$vectors[, seq_len(a - 1), drop = FALSE]
  }

  return(cont)
}

# add syntax for the formula
# add priors for the formula
# add data for the formula
# add inits for the formula
# variables to monitor
JAGS_formula(y ~ x_cont + x_fac3 * x_fac4, data = df, NULL)



formula <- formula(~ x_fac3 * x_fac4)
formula <- formula(~ x_fac3 * x_fac4 -1)
formula <- formula(z ~ x_cont + x_fac3 * x_fac4)
data <- df
parameter <- "mu"
prior_list <- list(
  "intercept"     = prior("normal", list(0, 1)),
  "x_cont"        = prior("normal", list(0, .1)),
  "x_fac3"        = prior_factor("normal", "orthonormal", list(1)),
  "x_fac3"        = prior_factor("normal", "orthonormal", list(1)),
  "x_fac4"        = prior_factor("normal", "orthonormal", list(2)),
  "x_fac3:x_fac4" = prior_factor("normal", "orthonormal", list(3))
)
length(formula)

attr(terms(formula), "response")
debug(lhs)
lhs(formula)
summary(lm(y ~ x_cont + x_bin + x_fac3 + x_fac4, data = df))


str(model.frame(y ~ x_cont + x_bin + x_fac3 + x_fac4, data = df))
str(model.frame(y ~ x_cont - x_bin + x_fac3 + x_fac4, data = df))





library(runjags)

counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
d.AD <- data.frame(treatment, outcome, counts)
glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())
temp <- template.jags(counts ~ outcome * treatment, data = d.AD, family = "poisson")
fit  <- run.jags("JAGSmodel.txt")

fit$model

bayestestR::contr.orthonorm




contrast_orthonorm <- function(n){
  # based on: bayestestR::contr.orthonorm
  # contr   <- stats::contr.treatment(n, contrasts = FALSE, base = 1, sparse = FALSE)
  # n       <- ncol(contr)
  # I_a     <- diag(n)
  # J_a     <- matrix(1, nrow = n, ncol = n)
  # Sigma_a <- I_a - J_a/n
  # contr   <- eigen(Sigma_a)$vectors[, seq_len(n - 1), drop = FALSE]

  cont    <- diag(n)
  a       <- n
  I_a     <- diag(a)
  J_a     <- matrix(1, nrow = a, ncol = a)
  Sigma_a <- I_a - J_a/a
  cont    <- eigen(Sigma_a)$vectors[,seq_len(a-1), drop = FALSE]


  return(cont)
}

Q2 <- contrast_orthonorm(2)
Q3 <- contrast_orthonorm(3)
Q5 <- contrast_orthonorm(5)

t(Q5[,c(2,1,3,4)])

t(bayestestR::contr.orthonorm(5))
t(stanova::contr.bayes(5))
t(stanova::contr.bayes(5))

diag(5) %*%  Q5

Q2 %*% diag(1) %*% t(Q2)
Q5 %*% diag(4) %*% t(Q5)

Q2 %*% t(Q2)
Q3 %*% t(Q3)
Q5 %*% t(Q5)

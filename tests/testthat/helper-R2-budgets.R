# A test shim that builds an R-squared budget from the fixed choices as
# separate arguments (law, convention, concentration, sigma, contrast) and an
# allocation given as concentrations or a Dirichlet prior, through the public
# prior_R2(model = R2_model(...), allocation = R2_dirichlet(...)) interface.
.R2_prior <- function(terms = NULL, total = prior("beta", list(1, 1)), allocation = NULL,
                      law = "exact", convention = "common", concentration = "inherited",
                      sigma = "sigma", contrast = "orthonormal", inclusion = NULL){

  if(is.null(contrast)){
    contrast <- "orthonormal"
  }
  if(is.prior(allocation)){
    allocation <- R2_dirichlet(allocation$parameters[["alpha"]])
  }else if(is.null(allocation) || is.numeric(allocation)){
    allocation <- R2_dirichlet(allocation)
  }
  prior_R2(terms = terms, total = total, allocation = allocation, inclusion = inclusion,
           model = R2_model(law = law, convention = convention, concentration = concentration,
                            sigma = sigma, contrast = contrast))
}

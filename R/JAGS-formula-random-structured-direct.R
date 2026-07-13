# Direct linear-time compilers for structured random-effect transforms.

.bt_JAGS_structured_dense_transform <- function(parameter, structure, K,
                                                n_groups, rho_name,
                                                sd_name,
                                                row_indexed_external_sd = FALSE,
                                                car_time_values = NULL){

  z_name      <- paste0(parameter, "_xRE_Zx")
  unit_name   <- paste0(parameter, "_xRE_UNIT_COEFx")
  coef_name   <- paste0(parameter, "_xRE_COEFx")
  syntax      <- character()

  if(K == 1L){
    syntax <- paste0(
      " for(g in 1:", n_groups, "){\n",
      "   ", unit_name, "[g,1] <- ", z_name, "[g,1]\n",
      if(isTRUE(row_indexed_external_sd)) "" else paste0(
        "   ", coef_name, "[g,1] <- ", sd_name, "[1] * ",
        unit_name, "[g,1]\n"
      ),
      " }\n"
    )
    return(syntax)
  }

  if(structure %in% c("cs", "hcs")){
    prefix_name <- paste0(parameter, "_xRE_CS_PREFIXx")
    diagonal_name <- paste0(parameter, "_xRE_CS_DIAGx")
    update_name <- paste0(parameter, "_xRE_CS_UPDATEx")
    syntax <- c(syntax, paste0(
      " for(i in 2:", K, "){\n",
      "   ", diagonal_name, "[i] <- sqrt((1 - ", rho_name,
      ") * (1 + (i - 1) * ", rho_name,
      ") / (1 + (i - 2) * ", rho_name, "))\n",
      "   ", update_name, "[i] <- ", rho_name,
      " * sqrt((1 - ", rho_name,
      ") / ((1 + (i - 2) * ", rho_name,
      ") * (1 + (i - 1) * ", rho_name, ")))\n",
      " }\n",
      " for(g in 1:", n_groups, "){\n",
      "   ", unit_name, "[g,1] <- ", z_name, "[g,1]\n",
      "   ", prefix_name, "[g,1] <- ", rho_name, " * ", z_name, "[g,1]\n",
      "   for(i in 2:", K, "){\n",
      "     ", unit_name, "[g,i] <- ", prefix_name,
      "[g,i - 1] + ", diagonal_name, "[i] * ", z_name, "[g,i]\n",
      "     ", prefix_name, "[g,i] <- ", prefix_name, "[g,i - 1] + ",
      update_name, "[i] * ", z_name, "[g,i]\n",
      "   }\n",
      if(isTRUE(row_indexed_external_sd)) "" else paste0(
        "   for(i in 1:", K, "){\n",
        "     ", coef_name, "[g,i] <- ", sd_name, "[i] * ", unit_name,
        "[g,i]\n",
        "   }\n"
      ),
      " }\n"
    ))
  }else{
    gap <- if(identical(structure, "car")){
      diff(car_time_values)
    }else{
      rep(1, K - 1L)
    }
    phi_name <- paste0(parameter, "_xRE_AR_PHIX")
    innovation_name <- paste0(parameter, "_xRE_AR_INNOVx")
    for(i in 2:K){
      syntax <- c(syntax, paste0(
        phi_name, "[", i, "] <- pow(", rho_name, ", ",
        .bt_JAGS_numeric_literal(gap[i - 1L]), ")\n",
        innovation_name, "[", i, "] <- sqrt(1 - pow(",
        phi_name, "[", i, "], 2))\n"
      ))
    }
    syntax <- c(syntax, paste0(
      " for(g in 1:", n_groups, "){\n",
      "   ", unit_name, "[g,1] <- ", z_name, "[g,1]\n",
      "   for(i in 2:", K, "){\n",
      "     ", unit_name, "[g,i] <- ", phi_name, "[i] * ",
      unit_name, "[g,i - 1] + ", innovation_name, "[i] * ", z_name,
      "[g,i]\n",
      "   }\n",
      if(isTRUE(row_indexed_external_sd)) "" else paste0(
        "   for(i in 1:", K, "){\n",
        "     ", coef_name, "[g,i] <- ", sd_name, "[i] * ", unit_name,
        "[g,i]\n",
        "   }\n"
      ),
      " }\n"
    ))
  }

  paste0(paste(syntax, collapse = ""), "\n")
}

.bt_JAGS_structured_local_transform <- function(parameter, layout, rho_name,
                                                sd_name,
                                                row_indexed_external_sd = FALSE){

  if(!inherits(layout, "BayesTools_random_effect_structured_local_layout")){
    stop("'layout' must be a structured local layout.", call. = FALSE)
  }
  check_bool(
    row_indexed_external_sd,
    "row_indexed_external_sd",
    allow_NA = FALSE
  )
  z_name      <- paste0(parameter, "_xRE_Zx")
  unit_name   <- paste0(parameter, "_xRE_UNIT_COEFx")
  coef_name   <- paste0(parameter, "_xRE_COEFx")
  prefix_name <- paste0(parameter, "_xRE_LOCAL_PREFIXx")
  syntax      <- character()
  prefix_flat <- 0L

  for(group in seq_along(layout$group_columns)){
    columns     <- layout$group_columns[[group]]
    if(length(columns) == 0L){
      next
    }
    coordinates <- layout$group_coordinates[[group]]
    first       <- columns[1L]
    syntax <- c(
      syntax,
      paste0(z_name, "[", group, ",", first, "] ~ dnorm(0, 1)"),
      paste0(unit_name, "[", group, ",", first, "] <- ",
             z_name, "[", group, ",", first, "]")
    )
    if(!isTRUE(row_indexed_external_sd)){
      syntax <- c(syntax, paste0(
        coef_name, "[", group, ",", first, "] <- ",
        sd_name, "[", first, "] * ", unit_name, "[", group, ",", first, "]"
      ))
    }
    if(layout$structure %in% c("cs", "hcs") && length(columns) > 1L){
      prefix_flat <- prefix_flat + 1L
      syntax <- c(syntax, paste0(
        prefix_name, "[", prefix_flat, "] <- ", rho_name, " * ",
        z_name, "[", group, ",", first, "]"
      ))
    }
    if(length(columns) == 1L){
      next
    }

    for(local in 2:length(columns)){
      column      <- columns[local]
      previous    <- columns[local - 1L]
      syntax <- c(syntax, paste0(
        z_name, "[", group, ",", column, "] ~ dnorm(0, 1)"
      ))
      if(layout$structure %in% c("cs", "hcs")){
        diagonal <- paste0(
          "sqrt((1 - ", rho_name, ") * (1 + ", local - 1L, " * ", rho_name,
          ") / (1 + ", local - 2L, " * ", rho_name, "))"
        )
        update <- paste0(
          rho_name, " * sqrt((1 - ", rho_name, ") / ((1 + ",
          local - 2L, " * ", rho_name, ") * (1 + ", local - 1L,
          " * ", rho_name, ")))"
        )
        syntax <- c(
          syntax,
          paste0(unit_name, "[", group, ",", column, "] <- ",
                 prefix_name, "[", prefix_flat, "] + ", diagonal, " * ",
                 z_name, "[", group, ",", column, "]")
        )
        if(local < length(columns)){
          previous_prefix <- prefix_flat
          prefix_flat     <- prefix_flat + 1L
          syntax <- c(syntax, paste0(
            prefix_name, "[", prefix_flat, "] <- ", prefix_name, "[",
            previous_prefix, "] + ", update, " * ", z_name, "[", group,
            ",", column, "]"
          ))
        }
      }else{
        gap <- coordinates[local] - coordinates[local - 1L]
        phi <- paste0(
          "pow(", rho_name, ", ", .bt_JAGS_numeric_literal(gap), ")"
        )
        syntax <- c(syntax, paste0(
          unit_name, "[", group, ",", column, "] <- ", phi, " * ",
          unit_name, "[", group, ",", previous, "] + sqrt(1 - pow(",
          phi, ", 2)) * ", z_name, "[", group, ",", column, "]"
        ))
      }
      if(!isTRUE(row_indexed_external_sd)){
        syntax <- c(syntax, paste0(
          coef_name, "[", group, ",", column, "] <- ", sd_name, "[", column,
          "] * ", unit_name, "[", group, ",", column, "]"
        ))
      }
    }
  }

  paste0(paste(syntax, collapse = "\n"), "\n")
}

.bt_JAGS_centered_independent_random <- function(parameter, K, n_groups,
                                                 sd_name){

  paste0(
    " for(g in 1:", n_groups, "){\n",
    "   for(i in 1:", K, "){\n",
    "     ", parameter, "_xRE_COEFx[g,i] ~ dnorm(0, pow(",
    sd_name, "[i], -2))\n",
    "     ", parameter, "_xRE_Zx[g,i] <- ",
    parameter, "_xRE_COEFx[g,i] / ", sd_name, "[i]\n",
    "   }\n",
    " }\n"
  )
}

.bt_JAGS_centered_correlated_random <- function(parameter, K, n_groups,
                                                sd_name, correlation_expression,
                                                cholesky_name){

  covariance_name <- paste0(parameter, "_xRE_COVx")
  syntax <- character()
  for(row in seq_len(K)){
    for(column in seq_len(K)){
      syntax <- c(syntax, paste0(
        covariance_name, "[", row, ",", column, "] <- ",
        sd_name, "[", row, "] * (", correlation_expression(row, column), ") * ",
        sd_name, "[", column, "]"
      ))
    }
  }
  syntax <- c(syntax, paste0(
    "for(g in 1:", n_groups, "){\n",
    "  ", parameter, "_xRE_COEFx[g,1:", K,
    "] ~ dmnorm.vcov(rep(0, ", K, "), ", covariance_name,
    "[1:", K, ",1:", K, "])\n",
    "  ", parameter, "_xRE_Zx[g,1] <- ",
    parameter, "_xRE_COEFx[g,1] / ", sd_name, "[1]\n"
  ))
  if(K > 1L){
    for(i in 2:K){
      previous <- paste0(
        cholesky_name, "[", i, ",", seq_len(i - 1L), "] * ",
        parameter, "_xRE_Zx[g,", seq_len(i - 1L), "]",
        collapse = " + "
      )
      syntax <- c(syntax, paste0(
        "  ", parameter, "_xRE_Zx[g,", i, "] <- ((",
        parameter, "_xRE_COEFx[g,", i, "] / ", sd_name, "[", i,
        "]) - (", previous, ")) / ", cholesky_name, "[", i, ",", i, "]"
      ))
    }
  }
  syntax <- c(syntax, "}")

  paste0(paste(syntax, collapse = "\n"), "\n")
}

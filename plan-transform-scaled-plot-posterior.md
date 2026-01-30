# Plan: Add `transform_scaled` to `plot_posterior`

## Overview

Add the ability to visualize prior and posterior distributions on both the **standardized** (scaled) and **original** (unscaled) scales when using formula-based models with auto-scaling.

## Problem Statement

When models use auto-scaling (standardizing predictors), the posterior samples are on the standardized scale. Users often want to visualize and interpret results on the original scale.

### Why Simple Linear Transformation Fails

For **coefficients** (main effects and interactions), a simple linear transformation `a + b*x` works because:
- Main effect: `β_orig = β_z / sd(x)` → `b = 1/sd(x)`, `a = 0`
- Interaction: `β_orig = β_z / (sd(x1) * sd(x2))` → `b = 1/(sd(x1)*sd(x2))`, `a = 0`

For the **intercept**, simple transformation FAILS because:
```
β₀_orig = β₀* - (mean(x₁)/sd(x₁)) * β₁* - (mean(x₂)/sd(x₂)) * β₂* - ...
```

The intercept on the original scale is a **weighted sum of multiple priors with different variances**:
- `Var(β₀_orig) = Var(β₀*) + (m₁/s₁)² * Var(β₁*) + (m₂/s₂)² * Var(β₂*) + ...`

If priors have different scales (e.g., intercept prior N(0,5) vs coefficient prior N(0,1)), the marginal prior for the intercept on the original scale is NOT a simple transformation of the intercept prior alone.

### Correct Approach

Generate samples from **ALL** priors simultaneously and apply the **same matrix transformation** used for posterior samples:
1. Sample from each prior (intercept, coefficients, interactions)
2. Form a matrix of prior samples (matching posterior structure)
3. Apply `transform_scale_samples()` to get prior samples on original scale
4. Extract marginal samples for the parameter of interest

This naturally handles all cases correctly, including the intercept.

---

## Proposed Solution

### New Function: `transform_prior_samples()`

Create a new function that generates prior samples and transforms them:

```r
#' Transform prior samples to original scale
#'
#' @param fit Fitted model object with prior_list and formula_scale attributes
#' @param n_samples Number of samples to generate (default: 10000)
#' @param seed Random seed for reproducibility
#' @return Matrix of prior samples on the original scale
transform_prior_samples <- function(fit, n_samples = 10000, seed = NULL)
```

**Implementation steps:**
1. Extract `prior_list` from fit
2. Extract `formula_scale` from fit
3. For each prior in `prior_list`, generate `n_samples` using `rng()` method
4. Arrange samples into matrix with same column structure as posterior
5. Apply `transform_scale_samples()` to transform to original scale
6. Return transformed prior samples matrix

### Modify `plot_posterior()`

Add `transform_scaled` argument to `plot_posterior()`:

```r
plot_posterior <- function(samples, parameter, ..., 
                           transform_scaled = FALSE,
                           formula_scale = NULL,
                           prior_samples = NULL)
```

When `transform_scaled = TRUE`:
1. Transform posterior samples using `transform_scale_samples()`
2. If `prior = TRUE`, use pre-generated `prior_samples` or generate via `transform_prior_samples()`
3. Plot both on the original scale

---

## Integration with `mix_posteriors` / `as_mixed_posteriors`

### Current State

- `mix_posteriors()` creates mixed posterior samples from multiple models
- Returns object with class `mixed_posteriors` containing:
  - `$samples` - mixed posterior samples
  - `$prior_list` - list of priors used
  - `$inference` - model weights and inference results

### Required Modifications

#### Option A: Extend `mix_posteriors()`

Add optional arguments:
```r
mix_posteriors(..., 
               transform_scaled = FALSE,
               formula_scale = NULL,
               generate_prior_samples = FALSE,
               n_prior_samples = 10000)
```

When `generate_prior_samples = TRUE`:
- Generate prior samples for each model
- Mix prior samples according to model weights (same as posterior mixing)
- Store in `$prior_samples` attribute

#### Option B: Create `as_mixed_posteriors()` enhancement

The function `as_mixed_posteriors()` can be enhanced to:
1. Accept single model fits and wrap them in mixed_posteriors format
2. Generate prior samples on demand
3. Handle formula_scale attribute propagation

**Note:** `as_mixed_posteriors()` may already support forcing prior samples for plotting - verify this capability and extend if needed.

---

## Implementation Plan

### Phase 1: Core Infrastructure

1. **Create `transform_prior_samples()`**
   - Location: `R/JAGS-formula.R` (near `transform_scale_samples`)
   - Inputs: fit object, n_samples, seed
   - Outputs: Matrix of transformed prior samples

2. **Create helper: `.generate_prior_sample_matrix()`**
   - Generate samples from all priors in prior_list
   - Match column structure to posterior samples
   - Handle different prior types (simple, factor, spike, etc.)

### Phase 2: Plot Integration

3. **Modify `plot_posterior()`** in `R/model-averaging-plots.R`
   - Add `transform_scaled` argument
   - Add `formula_scale` argument (optional, extracted from samples if available)
   - When `transform_scaled = TRUE`:
     - Transform posterior samples
     - Transform prior samples (if prior = TRUE)
     - Adjust axis labels to indicate "original scale"

### Phase 3: Mixed Posteriors Integration

5. **Extend `mix_posteriors()`**
   - Add `generate_prior_samples` argument
   - Store `prior_samples` in returned object
   - Propagate `formula_scale` attribute

6. **Extend `as_mixed_posteriors()`**
   - Handle single model case
   - Support prior sample generation
   - Ensure formula_scale is preserved

### Phase 4: Testing

7. **Unit tests** in `tests/testthat/test-JAGS-formula-scale.R`
   - Test `transform_prior_samples()` generates correct samples
   - Verify prior samples have expected distribution after transformation
   - Test with various model structures (main effects, interactions, intercept)

8. **Visual regression tests** in `tests/testthat/test-model-averaging-plots.R`
   - Side-by-side: scaled vs unscaled plots
   - Verify prior overlays match posterior distributions

---

## Technical Considerations

### Prior Sample Generation

Different prior types need different handling:
- **Simple priors** (normal, cauchy, etc.): Use `rng(prior, n_samples)`
- **Factor priors**: Generate for each level
- **Spike priors**: All samples = spike location
- **Mixture priors**: Sample component, then sample from component

### Column Matching

Prior sample matrix must match posterior column names exactly:
```r
# Posterior columns might be: c("mu_intercept", "mu_x1", "mu_x2", "mu_x1__xXx__x2", "sigma")
# Prior sample matrix must have same columns
```

### Missing Priors

Some posterior parameters may not have corresponding priors in `prior_list`:
- Handle gracefully (skip transformation or use identity)
- Document behavior clearly

### Performance

Generating many prior samples can be slow:
- Default to reasonable n_samples (10000)
- Allow caching/reuse of prior samples
- Consider lazy evaluation

---

## API Examples

### Basic Usage

```r
# Fit model with auto-scaling
fit <- JAGS_fit(data, formula = y ~ x1 * x2, ...)

# Extract posterior
posterior <- as.matrix(coda::as.mcmc.list(fit))

# Plot on standardized scale (current behavior)
plot_posterior(posterior, "mu_x1", prior = TRUE)

# Plot on original scale (new behavior)
plot_posterior(posterior, "mu_x1", prior = TRUE, 
               transform_scaled = TRUE, 
               fit = fit)  # or formula_scale = attr(fit, "formula_scale")
```

### With Mixed Posteriors

```r
# Mix posteriors from multiple models
mixed <- mix_posteriors(model_list, parameters = c("mu_x1", "mu_intercept"), ...)

# Plot on original scale with prior
plot_posterior(mixed, "mu_x1", 
               prior = TRUE,
               transform_scaled = TRUE)
```

---

## Open Questions

1. **Should `transform_scaled` default to `TRUE` or `FALSE`?**
   - FALSE maintains backward compatibility
   - TRUE might be more intuitive for users

2. **How to handle parameters not affected by scaling?**
   - Skip transformation (identity)
   - Or always apply (no-op for unscaled parameters)

3. **Axis labeling convention?**
   - Add "(original scale)" to axis labels?
   - Use parameter name transformation (e.g., "x1" vs "scale(x1)")?

4. **Integration with existing `transformation` argument?**
   - `transform_scaled` is separate from user-specified transformations (exp, log, etc.)
   - Should they be combinable? (first unscale, then apply user transformation)

---

## Files to Modify

| File | Changes |
|------|---------|
| `R/JAGS-formula.R` | Add `transform_prior_samples()`, `.generate_prior_sample_matrix()` |
| `R/model-averaging-plots.R` | Modify `plot_posterior()`, `plot_models()` |
| `R/model-averaging.R` | Extend `mix_posteriors()`, `as_mixed_posteriors()` |
| `tests/testthat/test-JAGS-formula-scale.R` | Add tests for prior sample transformation |
| `tests/testthat/test-model-averaging-plots.R` | Add visual tests for transformed plots |
| `man/plot_posterior.Rd` | Update documentation |
| `man/transform_prior_samples.Rd` | New documentation |

---

## Success Criteria

1. Prior and posterior distributions align correctly on both scales
2. Intercept transformation is correct (uses matrix transformation)
3. Coefficient transformations remain correct
4. Backward compatible (existing code works unchanged)
5. Clear documentation and examples
6. Comprehensive test coverage

remove incorrect: get_scale_transformation and the corresponding tests
---

## References

- Current implementation: `get_scale_transformation()` in `R/JAGS-formula.R`
- Matrix transformation: `.build_unscale_matrix()` in `R/JAGS-formula.R`
- Posterior transformation: `transform_scale_samples()` in `R/JAGS-formula.R`
- Plot functions: `R/model-averaging-plots.R`

---

## Detailed Implementation Analysis

### Existing Infrastructure (Code Review)

#### 1. `mix_posteriors()` (R/model-averaging.R:177)

**Current behavior:**
- Takes `model_list`, `parameters`, `is_null_list`
- Extracts `fits`, `priors` from each model
- Dispatches to type-specific helpers: `.mix_posteriors.simple()`, `.mix_posteriors.factor()`, etc.
- Each helper:
  - Samples from posterior based on `post_probs` (model weights)
  - For spike priors: uses `priors[[i]]$parameters[["location"]]` directly (no sampling)
  - Attaches `prior_list` as attribute to output

**Key insight:** The mixing logic samples from posteriors proportionally to model weights. The same logic can sample from priors.

**Relevant code pattern (`.mix_posteriors.simple`, line 370-378):**
```r
if(is.prior.point(priors[[i]])){
  # not sampling the priors as the samples would be already transformed
  samples <- c(samples, rep(priors[[i]]$parameters[["location"]], length(temp_ind)))
}else{
  samples <- c(samples, model_samples[temp_ind, parameter])
}
```

#### 2. `as_mixed_posteriors()` (R/model-averaging.R:707)

**Current behavior:**
- Takes single `BayesTools_fit` model
- Extracts `priors` from `attr(model, "prior_list")`
- Extracts `model_samples` from `coda::as.mcmc(model)`
- Applies conditioning if specified
- Dispatches to type-specific helpers: `.as_mixed_posteriors.simple()`, etc.
- Returns object with class `"as_mixed_posteriors"`, `"mixed_posteriors"`

**Key insight:** This wraps a single model's posterior as mixed_posteriors format. Ideal place to also generate prior samples.

**Note:** `force_plots` argument exists but appears minimally used (line 694-696):
> "temporal argument allowing to generate conditional posterior samples suitable for prior and posterior plots"

#### 3. `plot_posterior()` (R/model-averaging-plots.R:995)

**Current behavior:**
- Takes `samples` (mixed_posteriors object) and `parameter`
- If `prior = TRUE`:
  - Extracts `prior_list` from `attr(samples[[parameter]], "prior_list")`
  - Calls `.plot_data_prior_list.simple()` which uses `density()` on priors
  - The `density.prior()` function uses `rng()` when `force_samples = TRUE`

**Key insight:** Prior visualization already uses `rng()` via the `density()` method. The infrastructure for sampling exists.

#### 4. `rng.prior()` (R/priors.R:1129)

**Current behavior:**
- Generates random samples from any prior type
- Handles: simple priors, spike_and_slab, mixture, factor priors
- For mixtures: samples component first, then samples from that component

**Key insight:** This is the building block for generating prior samples. Already handles all prior types.

#### 5. `transform_scale_samples()` (R/JAGS-formula.R:1432)

**Current behavior:**
- Takes fit or matrix, extracts/uses `formula_scale`
- Calls `.apply_unscale_transform()` which uses `.build_unscale_matrix()`
- Applies full matrix transformation: `posterior_transformed = posterior %*% M^T`

**Key insight:** This can transform ANY matrix with the right column structure - posterior OR prior samples.

---

## Recommended Implementation Strategy

### Strategy: Minimal Changes with Maximum Reuse

The most maintainable approach leverages existing infrastructure:

1. **Reuse `rng.prior()`** for prior sample generation (already handles all types)
2. **Reuse `transform_scale_samples()`** for matrix transformation
3. **Extend `as_mixed_posteriors()`** to optionally generate prior samples
4. **Modify `plot_posterior()`** minimally to use transformed samples

### Implementation Details

#### Step 1: Create `.generate_prior_sample_matrix()` helper

Location: `R/JAGS-formula.R` (or `R/model-averaging.R` if preferred)

```r
.generate_prior_sample_matrix <- function(prior_list, n_samples, column_names = NULL, seed = NULL) {
  # Generate samples from all priors matching posterior column structure
  
  if (!is.null(seed)) set.seed(seed)
  
  # Initialize matrix
  n_params <- length(prior_list)
  samples <- matrix(NA, nrow = n_samples, ncol = n_params)
  colnames(samples) <- names(prior_list)
  
  for (i in seq_along(prior_list)) {
    prior <- prior_list[[i]]
    param_name <- names(prior_list)[i]
    
    if (is.null(prior) || is.prior.none(prior)) {
      samples[, i] <- 0  # No effect
    } else if (is.prior.point(prior)) {
      samples[, i] <- prior$parameters[["location"]]
    } else {
      samples[, i] <- rng(prior, n_samples)  # Uses existing rng.prior()
    }
  }
  
  # Reorder columns to match column_names if provided
  if (!is.null(column_names)) {
    # Match available columns
    available <- intersect(column_names, colnames(samples))
    samples <- samples[, available, drop = FALSE]
  }
  
  return(samples)
}
```

#### Step 2: Create `transform_prior_samples()` (exported function)

Location: `R/JAGS-formula.R`

```r
#' @title Transform prior samples to original scale
#' @description Generate prior samples and transform them using the same
#'   matrix transformation as posterior samples.
#' @param fit Fitted model with prior_list and formula_scale attributes
#' @param n_samples Number of samples to generate
#' @param seed Random seed for reproducibility
#' @return Matrix of prior samples on original scale
#' @export
transform_prior_samples <- function(fit, n_samples = 10000, seed = NULL) {
  
  prior_list <- attr(fit, "prior_list")
  formula_scale <- attr(fit, "formula_scale")
  
  if (is.null(prior_list)) {
    stop("'fit' must have 'prior_list' attribute")
  }
  
  # Get posterior column names for structure matching
  posterior <- as.matrix(coda::as.mcmc.list(fit))
  
  # Generate prior samples
  prior_samples <- .generate_prior_sample_matrix(
    prior_list, 
    n_samples = n_samples, 
    column_names = colnames(posterior),
    seed = seed
  )
  
  # Apply same transformation as posterior
  if (!is.null(formula_scale) && length(formula_scale) > 0) {
    prior_samples <- transform_scale_samples(prior_samples, formula_scale)
  }
  
  return(prior_samples)
}
```

#### Step 3: Extend `as_mixed_posteriors()`

Add `generate_prior_samples` argument:

```r
as_mixed_posteriors <- function(model, parameters, conditional = NULL, 
                                 conditional_rule = "AND", force_plots = FALSE,
                                 generate_prior_samples = FALSE,  # NEW
                                 n_prior_samples = 10000) {        # NEW
  # ... existing code ...
  
  # At the end, before return:
  if (generate_prior_samples) {
    prior_samples <- transform_prior_samples(model, n_samples = n_prior_samples)
    attr(out, "prior_samples") <- prior_samples
  }
  
  # Propagate formula_scale
  attr(out, "formula_scale") <- attr(model, "formula_scale")
  
  return(out)
}
```

#### Step 4: Modify `plot_posterior()`

Add `transform_scaled` argument:

```r
plot_posterior <- function(samples, parameter, plot_type = "base", prior = FALSE,
                           n_points = 1000, n_samples = 10000, force_samples = FALSE,
                           transform_scaled = FALSE,  # NEW
                           formula_scale = NULL,      # NEW
                           ...) {
  
  # Extract formula_scale from samples if not provided
  if (transform_scaled && is.null(formula_scale)) {
    formula_scale <- attr(samples, "formula_scale")
  }
  
  # ... existing parameter extraction ...
  
  # Transform posterior samples if requested
  if (transform_scaled && !is.null(formula_scale)) {
    # Get raw samples and transform
    raw_samples <- # extract from samples object
    transformed_samples <- transform_scale_samples(raw_samples, formula_scale)
    # Use transformed_samples for plotting
  }
  
  # For prior plotting with transform_scaled:
  if (prior && transform_scaled) {
    # Check if pre-generated prior_samples exist
    prior_samples <- attr(samples, "prior_samples")
    if (is.null(prior_samples)) {
      # Generate on the fly (need fit object)
      warning("Prior samples not pre-generated. For best results, use generate_prior_samples=TRUE in as_mixed_posteriors()")
      # Fall back to standard prior density (may be incorrect for intercept)
    }
    # Use prior_samples[, parameter] for density estimation
  }
  
  # ... rest of plotting logic ...
}
```

---

## Alternative: Simpler `plot_posterior` Enhancement

If modifying the mixed_posteriors infrastructure is too invasive, a simpler approach:

```r
plot_posterior <- function(samples, parameter, ...,
                           transform_scaled = FALSE,
                           fit = NULL) {  # Accept original fit object
  
  if (transform_scaled) {
    if (is.null(fit)) {
      stop("'fit' required when transform_scaled = TRUE")
    }
    
    # Transform posterior
    posterior_orig <- transform_scale_samples(fit)
    
    # For prior, generate and transform
    if (prior) {
      prior_samples_orig <- transform_prior_samples(fit, n_samples = n_samples)
      # Use density(prior_samples_orig[, parameter]) for prior overlay
    }
  }
}
```

This keeps all transformation logic in the plotting function, avoiding changes to `mix_posteriors` / `as_mixed_posteriors`.

---

## Recommendation

**Phase 1 (Immediate):** Implement the simpler approach - add `transform_scaled` and `fit` arguments to `plot_posterior()`. This:
- Requires minimal changes
- Keeps transformation logic centralized
- Works for single-model cases

**Phase 2 (Future):** Extend `as_mixed_posteriors()` with `generate_prior_samples` for:
- Pre-computed prior samples
- Multi-model averaging cases
- Better performance (compute once, plot multiple times)

---

## Factor Priors Consideration

For factor priors (orthonormal, meandif, treatment contrasts), the sample structure is more complex:
- Multiple columns per factor level
- Need to match column naming convention

The `rng.prior()` already handles this via `transform_factor_samples` argument.
Ensure `.generate_prior_sample_matrix()` correctly handles:
- `is.prior.factor()` → returns matrix, not vector
- Column naming: `parameter[1]`, `parameter[2]`, etc.

---

## Edge Cases to Handle

1. **Parameters not in prior_list**: Skip or use 0
2. **Parameters not affected by scaling**: Identity transformation (already handled by `transform_scale_samples`)
3. **Spike priors**: Generate constant samples at spike location
4. **Mixture priors**: Sample component, then sample from component (handled by `rng.prior`)
5. **Models without formula_scale**: Return untransformed samples

````

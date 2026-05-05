# Implemented finite-sample bounds for product-space inclusion BF table reporting

Timestamp: 2026-05-05T07:45:03.085764+00:00
Type: implementation-progress

Implemented BEAST/BSSVS-style reporting convention for product-space inclusion Bayes factors in `R/summary-tables.R`: when posterior inclusion probability is exactly 0 or 1 and prior inclusion probability is in (0, 1), table values now store/display finite-sample lower/upper bounds using 1/S or (S-1)/S, with '<'/'>' operator metadata and a table footnote. Bound metadata is preserved through `format_BF()` and `update.BayesTools_table()`, including inequality inversion for BF01. Updated BF Monte Carlo error column label to `error%(Inclusion BF)` and adjusted summary-table tests/reference outputs. Verification: unit profile passed; targeted summary/JAGS table fixture tests passed; full fixture profile stopped only in fixture-integrity due stale cached fit metadata.

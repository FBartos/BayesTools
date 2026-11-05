# Resolve the BayesTools installation for a dependent package: the CRAN release
# when it satisfies the package's `BayesTools (>= x)` import, otherwise the
# development source named in BAYESTOOLS_PAK_SPEC, and nothing for a package
# that does not import BayesTools (such as BayesTools itself). Run from the
# package's root; writes `pak_spec` to GITHUB_OUTPUT.

development_spec <- Sys.getenv("BAYESTOOLS_PAK_SPEC")
if (!nzchar(development_spec)) {
  stop("'BAYESTOOLS_PAK_SPEC' must not be empty.", call. = FALSE)
}

description <- read.dcf("DESCRIPTION")[1, ]
imports     <- if ("Imports" %in% names(description)) description[["Imports"]] else ""
imports_bayestools <- grepl("\\bBayesTools\\b", imports, perl = TRUE)

match <- regexec(
  "BayesTools\\s*\\(\\s*>=\\s*([^)[:space:]]+)\\s*\\)",
  imports,
  perl = TRUE
)

required_version <- regmatches(imports, match)[[1]]
required_version <- if (length(required_version) >= 2) required_version[[2]] else NA_character_

cran_version <- tryCatch(
  {
    repos <- getOption("repos")
    if (is.null(repos) ||
        is.na(repos[["CRAN"]]) ||
        identical(unname(repos[["CRAN"]]), "@CRAN@")) {
      repos <- c(CRAN = "https://cloud.r-project.org")
    }

    available <- utils::available.packages(
      contriburl = utils::contrib.url(repos, type = "source")
    )
    available["BayesTools", "Version"]
  },
  error = function(e) NA_character_
)

use_development <- !is.na(required_version) &&
  (is.na(cran_version) ||
     utils::compareVersion(required_version, cran_version) > 0)

pak_spec <- if (!imports_bayestools) {
  ""
} else if (use_development) {
  development_spec
} else {
  "any::BayesTools"
}

message(
  "BayesTools required: ", ifelse(is.na(required_version), "none", required_version),
  "; CRAN: ", ifelse(is.na(cran_version), "unavailable", cran_version),
  "; pak spec: ", pak_spec
)

output_file <- Sys.getenv("GITHUB_OUTPUT")
if (nzchar(output_file)) {
  cat("pak_spec=", pak_spec, "\n", sep = "", file = output_file, append = TRUE)
} else {
  cat(pak_spec, "\n", sep = "")
}

## Package-level documentation and namespace imports.
##
## fiber relies on base R's stats functions for its random draws and
## interpolation. Declaring them with @importFrom keeps R CMD check clean
## ("no visible global function definition for ...") and lets the functions
## resolve without stats:: qualification.

#' @keywords internal
"_PACKAGE"

#' @importFrom stats pgamma qgamma rbinom rgamma rnbinom runif
NULL

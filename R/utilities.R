#--------------------------------------#
#'@noRd
assert <- function (expr, error) {
  if (!expr) stop(error, call. = FALSE)
}
get_std <- function(object) {
  tmp <- lavaan::standardizedSolution(object,
                                      se = FALSE)
  pnames <- lavaan::lav_partable_labels(tmp)
  out <- tmp$est.std
  names(out) <- pnames
  out
}
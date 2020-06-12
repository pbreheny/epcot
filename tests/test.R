if (requireNamespace("tinytest", quietly=TRUE)) {
  home <- length(unclass(packageVersion("epcot"))[[1]]) == 4
  if (home) tinytest::test_package("epcot")
}

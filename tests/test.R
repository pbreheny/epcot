if (requireNamespace("tinytest", quietly=TRUE)) {
  tinytest::test_package("epcot", pattern="^[^_].*\\.[rR]$")
}

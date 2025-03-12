# semboottools 0.0.0.9010

* Major functions finalized. Ready to
  use. (0.0.0.9001)

* Renamed `plot_boot()` to
  `hist_qq_boot()` to avoid name
  conflict with functions with
  similar names. (0.0.0.9002)

* Initial setup of the `pkgdown` website.
  (0.0.0.9003)

* Added asymmetric bootstrap *p*-values
  to `standardizedSolution_boot()`.
  Enabled by default.
  (0.0.0.9004)

* Changed the default output format
  of `print.sbt_std_boot()` to `"text"`.
  (0.0.0.9004)

* Updated the maintainer email address.
  (0.0.0.9005)

* Added `parameterEstimates_boot()`, with
  a print method. Used `lavaan.printer`
  by default. Also updated the
  `pkgdown` site. (0.0.0.9006)

* Updated `store_boot()` to store more
  information. (0.0.0.9007)

* Updated related functions such that
  `store_boot()` will no longer be
  implicitly called to do bootstrapping.
  (0.0.0.9007)

* Enabled more tests to increase test
  coverage. (0.0.0.9008)

* Updated the print method of the output
  of `standardizedSolution_boot()` to
  use `lavaan.printer` by default.
  (0.0.0.9009)

* Updated `parameterEstimates_boot()`
  for a bug with multigroup models.
  (0.0.0.9009)

* Updated `hist_qq_boot()` and
  `scatter_boot()` to support the
  output of `parameterEstimates_boot()`.
  (0.0.0.9010)

# dsos 0.1.2

## New features

* Added functions for Bayesian test (`bf_from_os`), the corresponding S3 methods
(`print` and `plot`) and helper functions to convert bayes factors to p-values
(`as_pvalue` and `as_bf`).

## Breaking changes

* Deleted custom scoring functions (`score_*`) which depended on
`ranger` and `isotree`. Users can now pass a scoring function as an argument
to recreate this behaviour.

## Enhancements

* Removed dependency on `ranger` and `isotree` packages.

* Updated README with Bayesian example and new links.

# dsos 0.1.1

## Breaking changes

* Changed naming convention for functions. Used the prefixes `pt_*`, `at_*` and
`score_*` to denote permutation test, asymptotic test and scoring
functions.

## New features

* Exposed functions to test for no adverse shift from outlier scores
directly: `pt_from_os` and `at_from_os`

## Enhancements

* Added S3 method to print the result of (statistical) tests.

* Increased font size for p-value in S3 method to plot result of
(statistical) tests.

* Updated README with new examples and links.

* Removed dependency on `WeightedROC` package.

* Made dependencies on `ranger` and `isotree` optional, required only when
using specific scoring functions.

# dsos 0.1.0

## First Public Release

* Added a `NEWS.md` file to track changes to the package.
* Initial public release to GitHub and CRAN.

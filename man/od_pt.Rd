% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score-outlier-detection.R
\name{od_pt}
\alias{od_pt}
\title{Dataset Shift via Isolation Scores}
\usage{
od_pt(x_train, x_test, R = 1000, num_trees = 500, sub_ratio = 1/2)
}
\arguments{
\item{x_train}{Training sample.}

\item{x_test}{Test sample.}

\item{R}{The number of permutations. May be ignored.}

\item{num_trees}{The number of trees in random forests.}

\item{sub_ratio}{Subsampling ratio for sample splitting. May be ignored.}
}
\value{
A named list or object of class \code{outlier.test} containing:
\itemize{
   \item \code{statistic}: observed WAUC statistic
   \item \code{seq_mct}: sequential Monte Carlo test, if applicable
   \item \code{p_value}: p-value
   \item \code{outlier_scores}: outlier scores from training and test set
}
}
\description{
Test for no adverse shift via isolation scores for two-sample comparison.
The scores are predictions from extended isolation forest with the package
\pkg{isotree}. The prefix \emph{od} stands for outlier detection, the
relevant notion of outlyingness.
}
\details{
The empirical null distribution uses \code{R} permutations to estimate
the p-value. For speed, this is implemented as a sequential Monte Carlo test
with the \pkg{simctest} package. See Gandy (2009) for details. The suffix
\emph{pt} refers to permutation test. It does not use the asymptotic
(theoretical) null distribution for the weighted AUC (WAUC), the test
statistic. This is the recommended approach for small samples.
}
\section{Notes}{

Isolation forest detects \emph{isolated} points, instances that are typically
out-of-distribution relative to the high-density regions of the data distribution.
Any performant method for density-based out-of-distribution detection can
replace isolation forest, the default in this implementation.
}

\examples{
\donttest{
library(dsos)
set.seed(12345)
data(iris)
x_train <- iris[1:50, 1:4] # Training sample: Species == 'setosa'
x_test <- iris[51:100, 1:4] # Test sample: Species == 'versicolor'
iris_test <- od_pt(x_train, x_test) # Can also use: od_ss
str(iris_test)
}
}
\references{
Kamulete, V. M. (2021).
\emph{Test for non-negligible adverse shifts}.
arXiv preprint arXiv:2107.02990.

Liu, F. T., Ting, K. M., & Zhou, Z. H. (2008, December).
\emph{Isolation forest}.
In 2008 Eighth IEEE International Conference on Data Mining (pp. 413-422).
IEEE.

Gandy, A. (2009).
\emph{Sequential implementation of Monte Carlo tests with uniformly bounded resampling risk}.
Journal of the American Statistical Association, 104(488), 1504-1511.

Li, J., & Fine, J. P. (2010).
\emph{Weighted area under the receiver operating characteristic curve and its application to gene selection}.
Journal of the Royal Statistical Society: Series C (Applied Statistics), 59(4), 673-692.
}
\concept{anomalies}

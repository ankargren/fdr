fdr
======

[![Build Status](https://travis-ci.org/ankargren/fdr.svg?branch=master)](https://travis-ci.org/ankargren/fdr) [![](http://www.r-pkg.org/badges/version/fdr)](http://www.r-pkg.org/pkg/fdr) [![Coverage status](https://codecov.io/gh/ankargren/fdr/master/graph/badge.svg)](https://codecov.io/github/ankargren/fdr?branch=master)

Overview
--------

The `fdr` package implements random number generation and density evaluation for the multivariate normal distribution, the matrix normal and matrix t distributions as well as the inverse Wishart distribution. By avoiding argument checks and a C++ implementation, the functions are quick and useful in e.g. various Monte Carlo sampling approaches (such as importance sampling and MCMC).

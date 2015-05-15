# ndarray-blas-trsv-complex

[![Build Status](https://travis-ci.org/scijs/ndarray-blas-trsv-complex.svg?branch=master)](https://travis-ci.org/scijs/ndarray-blas-trsv-complex) [![npm version](https://badge.fury.io/js/ndarray-blas-trsv-complex.svg)](http://badge.fury.io/js/ndarray-blas-trsv-complex)

BLAS Level 2 TRSV (triangular solve) for complex [ndarrays](https://github.com/scijs/ndarray)

## Usage

##### `trsv( A_r, A_i, x_r, x_i [, uplo] )`
Calculate `x <- A^-1 x` for the real and complex parts `A_r` and `A_i` of the upper triangular matrix A using back-substitution. Data below the diagonal is ignored.  If `uplo` is 'lo', uses the lower triangular portion of A and performs forward-substitution instead. Result overwrites the x vectors.

## Credits
(c) 2015 Ricky Reusser. MIT License


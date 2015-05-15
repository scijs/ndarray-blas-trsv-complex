'use strict'

var blas1 = require('ndarray-blas-level1-complex'),
    div = require('complex-division')

function trsv (A_r, A_i, x_r, x_i, uplo) {
  var nm1, i, d, result,
      dotu = blas1.dotu,
      m = A_r.shape[0],
      n = A_r.shape[1],
      result

  if( A_r.dimension!==2 || A_i.dimension!==2 || m!==n || A_i.shape[0]!==m || A_i.shape[1]!==n) {
    throw new TypeError('trsv():: A_r and A_i must be square matrices of equal size')
  }

  if( x_r.dimension!==1 ||x_i.dimension!==1 || x_r.shape[0]!==A_r.shape[0] || x_r.shape[0]!==x_i.shape[0] ) {
    throw new TypeError('trsv():: x_r and x_i must be vectors with size equal to the width of A')
  }

  switch(uplo) {
    case 'lo':
      result = div( x_r.get(0), x_i.get(0), A_r.get(0,0), A_i.get(0,0) )
      x_r.set(0, result[0])
      x_i.set(0, result[1])
      for(i=1; i<n; i++) {
        result = dotu( A_r.pick(i,null).hi(i), A_i.pick(i,null).hi(i), x_r.hi(i), x_i.hi(i) )
        result = div( x_r.get(i) - result[0], x_i.get(i) - result[1], A_r.get(i,i), A_i.get(i,i) )
        x_r.set(i, result[0])
        x_i.set(i, result[1])
      }
      break
    case 'up':
    default:
      nm1 = n-1
      result = div( x_r.get(nm1), x_i.get(nm1), A_r.get(nm1,nm1), A_i.get(nm1,nm1) )
      x_r.set(nm1, result[0])
      x_i.set(nm1, result[1])
      for(i=n-2; i>=0; i--) {
        result = dotu( A_r.pick(i,null).lo(i+1), A_i.pick(i,null).lo(i+1), x_r.lo(i+1), x_i.lo(i+1) )
        result = div( x_r.get(i) - result[0], x_i.get(i) - result[1], A_r.get(i,i), A_i.get(i,i) )
        x_r.set(i, result[0])
        x_i.set(i, result[1])
      }
  }
  return true
}

module.exports = trsv

'use strict';

var assert = require('chai').assert,
    ndarray = require('ndarray'),
    trsv = require('../'),
    ndt = require('ndarray-tests')

describe("trsv()",function() {

  var A_r, A_i, A0_r, A0_i, b_r, b_i, b0_r, b0_i, c_r, c_i, c0_r, c0_i

  beforeEach(function() {
    var Ardata = [1,-3,0,7]
    var Aidata = [2,4,0,-2]
    var brdata = [12,-1]
    var bidata = [24,-30]
    var crdata = [-1,-18]
    var cidata = [8,-24]
    A_r = ndarray(new Float64Array(Ardata), [2,2])
    A_i = ndarray(new Float64Array(Aidata), [2,2])
    A0_r = ndarray(new Float64Array(Ardata), [2,2])
    A0_i = ndarray(new Float64Array(Aidata), [2,2])
    b_r = ndarray(new Float64Array(brdata))
    b_i = ndarray(new Float64Array(bidata))
    b0_r = ndarray(new Float64Array(brdata))
    b0_i = ndarray(new Float64Array(bidata))
    c_r = ndarray(new Float64Array(crdata))
    c_i = ndarray(new Float64Array(cidata))
    c0_r = ndarray(new Float64Array(crdata))
    c0_i = ndarray(new Float64Array(cidata))
  })

  it('upper triangular',function() {
    assert( trsv( A_r, A_i, b_r, b_i ), 'returns true on success' )
    assert( ndt.approximatelyEqual( ndarray([3,1]), b_r, 1e-8 ), 'Real component of solution is correct' )
    assert( ndt.approximatelyEqual( ndarray([2,-4]), b_i, 1e-8 ), 'Imag component of solution is correct' )
  })

  it('lower triangular',function() {
    assert( trsv( A_r.transpose(1,0), A_i.transpose(1,0), c_r, c_i, 'lo' ), 'returns true on success' )
    assert( ndt.approximatelyEqual( ndarray([3,1]), c_r, 1e-8 ), 'Real component of solution is correct' )
    assert( ndt.approximatelyEqual( ndarray([2,-4]), c_i, 1e-8 ), 'Imag component of solution is correct' )
  })

})


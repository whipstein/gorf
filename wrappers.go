package gorf

/*
#cgo LDFLAGS: libblas.a liblapack.a -L/usr/local/Cellar/gcc/11.2.0/lib/gcc/11 -lgfortran
#include <stdlib.h>
#include <stdbool.h>
#include <complex.h>
void zgeev_(char*, char*, int*, complex double*, int*, complex double*, complex double*, int*, complex double*, int*, complex double*, int*, double*, int*);
*/
import "C"

func _Zgeev(jobvl byte, jobvr byte, n int, a []complex128, lda int, w []complex128, vl []complex128, ldvl int, vr []complex128, ldvr int, work []complex128, lwork int, rwork []float64, info *int) {
	var _a, _w, _vl, _vr, _work *complex128
	var _rwork *float64
	_jobvl := (C.char)(jobvl)
	_jobvr := (C.char)(jobvr)
	_n := (C.int)(n)
	_a = &a[0]
	_lda := (C.int)(lda)
	_w = &w[0]
	_vl = &vl[0]
	_ldvl := (C.int)(ldvl)
	_vr = &vr[0]
	_ldvr := (C.int)(ldvr)
	_work = &work[0]
	_lwork := (C.int)(lwork)
	_rwork = &rwork[0]
	_info := (C.int)(*info)
	C.zgeev_(&_jobvl, &_jobvr, &_n, (*C.complexdouble)(_a), &_lda, (*C.complexdouble)(_w), (*C.complexdouble)(_vl), &_ldvl, (*C.complexdouble)(_vr), &_ldvr, (*C.complexdouble)(_work), &_lwork, (*C.double)(_rwork), &_info)
}

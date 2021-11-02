package gorf

import (
	"strconv"

	"github.com/whipstein/golinalg/golapack"
	"github.com/whipstein/golinalg/mat"
)

const (
	epsf64 = 2.2204460492503131e-016
	tiny   = 2.2250738585072014e-308
	huge   = 1.7976931348623157e+308
	radix  = 2
	digits = 53
	minexp = -1021
	maxexp = 1024
)

func absint(a int) int {
	if a < 0 {
		return -a
	}
	return a
}

func maxf64(a ...float64) float64 {
	maxval := a[0]
	for _, val := range a {
		if val > maxval {
			maxval = val
		}
	}
	return maxval
}

func maxint(a ...int) int {
	maxval := a[0]
	for _, val := range a {
		if val > maxval {
			maxval = val
		}
	}
	return maxval
}

func minf64(a ...float64) float64 {
	minval := a[0]
	for _, val := range a {
		if val < minval {
			minval = val
		}
	}
	return minval
}

func minint(a ...int) int {
	minval := a[0]
	for _, val := range a {
		if val < minval {
			minval = val
		}
	}
	return minval
}

func powint(a, b int) int {
	if b == 0 {
		return 1
	}

	result := a
	for i := 1; i < b; i++ {
		result *= a
	}
	return result
}

func matInv(m *mat.Matrix) {
	var colMajor int

	r, c := m.Rows, m.Cols
	work := mf(r, c)
	ipiv := make([]int, r)

	if m.Opts.Major == mat.Row {
		m.ToColMajor()
		colMajor = 1
	}
	if info, err := golapack.Dgetrf(r, c, m, &ipiv); err != nil || info != 0 {
		panic("golapack.Zgetrf error: " + strconv.Itoa(info))
	}
	if info, err := golapack.Dgetri(r, m, ipiv, work); err != nil || info != 0 {
		panic("golapack.Zgetri error: " + strconv.Itoa(info))
	}
	if colMajor == 1 {
		m.ToRowMajor()
		colMajor = 0
	}
}

func cmatInv(m *mat.CMatrix) *mat.CMatrix {
	var colMajor int

	r, c := m.Rows, m.Cols
	work := cvf(r * c)
	lwork := r * c
	ipiv := make([]int, r)

	if m.Opts.Major == mat.Row {
		m.ToColMajor()
		colMajor = 1
	}
	if info, err := golapack.Zgetrf(r, c, m, &ipiv); err != nil || info != 0 {
		panic("golapack.Zgetrf error: " + strconv.Itoa(info))
	}
	if info, err := golapack.Zgetri(r, m, &ipiv, work, lwork); err != nil || info != 0 {
		panic("golapack.Zgetri error: " + strconv.Itoa(info))
	}
	if colMajor == 1 {
		m.ToRowMajor()
	}

	return m
}

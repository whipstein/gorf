package gorf

import (
	"fmt"
	"math"
	"math/cmplx"

	"github.com/whipstein/golinalg/goblas"
	"github.com/whipstein/golinalg/golapack"
	"github.com/whipstein/golinalg/mat"
)

func convertNumber(a, b float64, encin, encout Encoding) (float64, float64) {
	if encin == encout {
		return a, b
	}

	switch encin {
	case RI:
		if encout == DB {
			return ri2db(a, b)
		} else if encout == MA {
			return ri2ma(a, b)
		}
	case MA:
		if encout == DB {
			return ma2db(a, b)
		} else if encout == RI {
			return ma2ri(a, b)
		}
	case DB:
		if encout == MA {
			return db2ma(a, b)
		} else if encout == RI {
			return db2ri(a, b)
		}
	}

	return 0, 0
}

func db2ma(m, a float64) (float64, float64) {
	return math.Pow(10, m/20), a
}

func db2ri(m, a float64) (float64, float64) {
	c := cmplx.Rect(math.Pow(10, m/20), a*math.Pi/180)
	return real(c), imag(c)
}

func ma2db(m, a float64) (float64, float64) {
	return 20 * math.Log10(m), a
}

func ma2ri(m, a float64) (float64, float64) {
	c := cmplx.Rect(m, a*math.Pi/180)
	return real(c), imag(c)
}

func ri2db(r, i float64) (float64, float64) {
	m, a := cmplx.Polar(complex(r, i))
	return 20 * math.Log10(m), a * 180 / math.Pi
}

func ri2ma(r, i float64) (float64, float64) {
	m, a := cmplx.Polar(complex(r, i))
	return m, a * 180. / math.Pi
}

func Passivity(m *mat.CMatrix) (p float64) {
	var info int
	var err error

	dummy := cmf(m.Rows, m.Cols, m.Opts)
	lambda := cvf(m.Rows)
	lwork := 3 * m.Rows
	work := cvf(lwork)
	rwork := vf(lwork)
	pcmplx := m.DeepCopy()
	golapack.Zlaset(mat.Full, m.Rows, m.Cols, 0, 1, pcmplx)
	if err = goblas.Zgemm(mat.ConjTrans, mat.NoTrans, m.Rows, m.Cols, m.Cols, -1, m, m, 1, pcmplx); err != nil {
		panic(err)
	}
	// if info, err = golapack.Zgeev('N', 'V', m.Rows, pcmplx, lambda, dummy, dummy, work, lwork, rwork); err != nil {
	// 	panic(err)
	// }
	_Zgeev('N', 'N', m.Rows, pcmplx.Data, m.Rows, lambda.Data, dummy.Data, m.Rows, dummy.Data, m.Rows, work.Data, lwork, rwork.Data, &info)
	if info != 0 {
		panic("Eigenvalue calculation failed!\n")
	}
	p = lambda.GetRe(0)
	for i := 1; i < lambda.Size; i++ {
		if lambda.GetRe(i) != math.NaN() && p > lambda.GetRe(i) {
			p = lambda.GetRe(i)
		}
	}
	return p
}

func AtoH(m *mat.CMatrix) *mat.CMatrix {
	h := m.DeepCopy()

	h.Set(0, 0, m.Get(0, 1)/m.Get(1, 1))
	h.Set(0, 1, (m.Get(0, 0)*m.Get(1, 1)-m.Get(0, 1)*m.Get(1, 0))/m.Get(1, 1))
	h.Set(1, 0, -1/m.Get(1, 1))
	h.Set(1, 1, m.Get(1, 0)/m.Get(1, 1))

	m.Data = h.Data

	return m
}

func AtoS(m *mat.CMatrix, z0 *mat.CVector) *mat.CMatrix {
	s := m.DeepCopy()

	r01 := z0.GetReCmplx(0)
	r02 := z0.GetReCmplx(1)
	z01 := z0.Get(0)
	z02 := z0.Get(1)
	z01conj := z0.GetConj(0)
	z02conj := z0.GetConj(1)

	denom := m.Get(0, 0)*z02 + m.Get(0, 1) + m.Get(1, 0)*z01*z02 + m.Get(1, 1)*z01

	s.Set(0, 0, (m.Get(0, 0)*z02+m.Get(0, 1)-m.Get(1, 0)*z01conj*z02-m.Get(1, 1)*z01conj)/denom)
	s.Set(0, 1, (2*(m.Get(0, 0)*m.Get(1, 1)-m.Get(0, 1)*m.Get(1, 0))*cmplx.Sqrt(r01*r02))/denom)
	s.Set(1, 0, (2*cmplx.Sqrt(r01*r02))/denom)
	s.Set(1, 1, (-m.Get(0, 0)*z02conj+m.Get(0, 1)-m.Get(1, 0)*z01*z02conj+m.Get(1, 1)*z01)/denom)
	m.Data = s.Data

	return m
}

func AtoT(m *mat.CMatrix, z0 *mat.CVector) *mat.CMatrix {
	m = StoT(AtoS(m, z0), z0)
	return m
}

func AtoY(m *mat.CMatrix) *mat.CMatrix {
	y := m.DeepCopy()

	y.Set(0, 0, m.Get(1, 1)/m.Get(0, 1))
	y.Set(0, 1, (m.Get(0, 1)*m.Get(1, 0)-m.Get(0, 0)*m.Get(1, 1))/m.Get(0, 1))
	y.Set(1, 0, -1/m.Get(0, 1))
	y.Set(1, 1, m.Get(0, 0)/m.Get(0, 1))
	m.Data = y.Data

	return m
}

func AtoZ(m *mat.CMatrix) *mat.CMatrix {
	z := m.DeepCopy()

	z.Set(0, 0, m.Get(0, 0)/m.Get(1, 0))
	z.Set(0, 1, (m.Get(0, 0)*m.Get(1, 1)-m.Get(0, 1)*m.Get(1, 0))/m.Get(1, 0))
	z.Set(1, 0, 1/m.Get(1, 0))
	z.Set(1, 1, m.Get(1, 1)/m.Get(1, 0))
	m.Data = z.Data

	return m
}

func HtoA(m *mat.CMatrix) *mat.CMatrix {
	a := m.DeepCopy()

	a.Set(0, 0, (m.Get(0, 1)*m.Get(1, 0)-m.Get(0, 0)*m.Get(1, 1))/m.Get(1, 0))
	a.Set(0, 1, -m.Get(0, 0)/m.Get(1, 0))
	a.Set(1, 0, -m.Get(1, 1)/m.Get(1, 0))
	a.Set(1, 1, -1/m.Get(1, 0))

	m.Data = a.Data

	return m
}

func HtoS(m *mat.CMatrix, z0 *mat.CVector) *mat.CMatrix {
	s := m.DeepCopy()

	r01 := z0.GetReCmplx(0)
	r02 := z0.GetReCmplx(1)
	z01 := z0.Get(0)
	z02 := z0.Get(1)
	z01conj := z0.GetConj(0)
	z02conj := z0.GetConj(1)

	denom := (z01+m.Get(0, 0))*(1+m.Get(1, 1)*z02) - m.Get(0, 1)*m.Get(1, 0)*z02

	s.Set(0, 0, ((m.Get(0, 0)-z01conj)*(1+m.Get(1, 1)*z02)-m.Get(0, 1)*m.Get(1, 0)*z02)/denom)
	s.Set(0, 1, (2*m.Get(0, 1)*cmplx.Sqrt(r01*r02))/denom)
	s.Set(1, 0, (-2*m.Get(1, 0)*cmplx.Sqrt(r01*r02))/denom)
	s.Set(1, 1, ((z01+m.Get(0, 0))*(1-m.Get(1, 1)*z02conj)+m.Get(0, 1)*m.Get(1, 0)*z02conj)/denom)
	m.Data = s.Data

	return m
}

func HtoT(m *mat.CMatrix, z0 *mat.CVector) *mat.CMatrix {
	m = StoT(HtoS(m, z0), z0)
	return m
}

func HtoY(m *mat.CMatrix) *mat.CMatrix {
	y := m.DeepCopy()

	y.Set(0, 0, 1/m.Get(0, 0))
	y.Set(0, 1, -m.Get(0, 1)/m.Get(0, 0))
	y.Set(1, 0, m.Get(1, 0)/m.Get(0, 0))
	y.Set(1, 1, (m.Get(0, 0)*m.Get(1, 1)-m.Get(0, 1)*m.Get(1, 0))/m.Get(0, 0))
	m.Data = y.Data

	return m
}

func HtoZ(m *mat.CMatrix) *mat.CMatrix {
	z := m.DeepCopy()

	z.Set(0, 0, (m.Get(0, 0)*m.Get(1, 1)-m.Get(0, 1)*m.Get(1, 0))/m.Get(1, 1))
	z.Set(0, 1, m.Get(0, 1)/m.Get(1, 1))
	z.Set(1, 0, -m.Get(1, 0)/m.Get(1, 1))
	z.Set(1, 1, 1/m.Get(1, 1))
	m.Data = z.Data

	return m
}

func StoA(s *mat.CMatrix, z0 *mat.CVector) *mat.CMatrix {
	var err error

	if (s.Rows != s.Cols) && (s.Rows != 2) {
		err = fmt.Errorf("Network must only have 2 ports\n")
		panic(err)
	}

	a := s.DeepCopy()

	z01 := z0.Get(0)
	z02 := z0.Get(1)
	z01conj := z0.GetConj(0)
	z02conj := z0.GetConj(1)

	denom := 2 * s.Get(1, 0) * cmplx.Sqrt(z0.GetReCmplx(0)*z0.GetReCmplx(1))

	a.Set(0, 0, ((z01conj+s.Get(0, 0)*z01)*(1-s.Get(1, 1))+s.Get(0, 1)*s.Get(1, 0)*z01)/denom)
	a.Set(0, 1, ((z01conj+s.Get(0, 0)*z01)*(z02conj+s.Get(1, 1)*z02)-s.Get(0, 1)*s.Get(1, 0)*z01*z02)/denom)
	a.Set(1, 0, ((1-s.Get(0, 0))*(1-s.Get(1, 1))-s.Get(0, 1)*s.Get(1, 0))/denom)
	a.Set(1, 1, ((1-s.Get(0, 0))*(z02conj+s.Get(1, 1)*z02)+s.Get(0, 1)*s.Get(1, 0)*z02)/denom)

	s.Data = a.Data

	return s
}

func StoH(s *mat.CMatrix, z0 *mat.CVector) *mat.CMatrix {
	var err error

	if (s.Rows != s.Cols) && (s.Rows != 2) {
		err = fmt.Errorf("Network must only have 2 ports\n")
		panic(err)
	}

	s = ZtoH(StoZ(s, z0))

	return s
}

func StoT(s *mat.CMatrix, z0 *mat.CVector) *mat.CMatrix {
	var err error

	if (s.Rows != s.Cols) && ((s.Rows % 2) != 0) {
		err = fmt.Errorf("Network does not have an even number of ports\n")
		panic(err)
	}

	rh, ch := s.Rows/2, s.Cols/2
	minv := cmf(rh, ch, s.Opts)
	mtmp := cmf(rh, ch, s.Opts)
	mtmp2 := cmf(rh, ch, s.Opts)
	t := s.DeepCopy()
	// Copy S matrix into minv for calculations
	golapack.Zlacpy(mat.Full, rh, ch, s.Off(rh, 0), minv)

	cmatInv(minv)

	// Calculate T[0:r/2, 0:c/2]
	//     mtmp=minv*S[r/2:r, c/2:c]
	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, rh, ch, ch, 1, minv, s.Off(rh, ch), 0, mtmp); err != nil {
		panic(err)
	}
	//     mtmp2=-1*S[0:r/2, 0:c/2]*mtmp + S[0:r/2, c/2:c]
	golapack.Zlacpy(mat.Full, rh, ch, s.Off(0, ch), mtmp2)
	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, rh, ch, ch, -1, s, mtmp, 1, mtmp2); err != nil {
		panic(err)
	}
	golapack.Zlacpy(mat.Full, rh, ch, mtmp2, t)

	// Calculate T[0:r/2, c/2:c]
	//     mtmp=S[0:r/2, 0:c/2]*minv
	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, rh, ch, ch, 1, s, minv, 0, t.Off(0, ch)); err != nil {
		panic(err)
	}

	// Calculate T[r/2:r, 0:c/2]
	//     mtmp=-1*minv*S[r/2:r, c/2:c]
	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, rh, ch, ch, -1, minv, s.Off(rh, ch), 0, t.Off(rh, 0)); err != nil {
		panic(err)
	}

	// Calculate T[r/2:r, c/2:c]
	golapack.Zlacpy(mat.Full, rh, ch, minv, t.Off(rh, ch))

	s.Data = t.Data

	return s
}

// y = sqrty0 * (Id - s) * (Id + s)**-1 * sqrty0
func StoY(s *mat.CMatrix, z0 *mat.CVector) *mat.CMatrix {
	var err error

	y := s.DeepCopy()
	yinv := s.DeepCopy()
	ytmp := s.DeepCopy()

	// Generate sqrt(Z0) matrix and identity matrix for calculation
	sqy0 := cmf(z0.Size, z0.Size, opts)
	for i, val := range z0.Data {
		sqy0.Set(i, i, 1/cmplx.Sqrt(val))
	}
	id := cmf(s.Rows, s.Cols, s.Opts)
	golapack.Zlaset(mat.Full, s.Rows, s.Cols, 0, 1, id)

	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, s.Rows, s.Cols, s.Cols, 1, id, id, -1, y); err != nil {
		panic(err)
	}
	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, s.Rows, s.Cols, s.Cols, 1, id, id, 1, yinv); err != nil {
		panic(err)
	}

	cmatInv(yinv)

	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, s.Rows, s.Cols, s.Cols, 1, sqy0, yinv, 0, ytmp); err != nil {
		panic(err)
	}
	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, s.Rows, s.Cols, s.Cols, 1, ytmp, y, 0, yinv); err != nil {
		panic(err)
	}
	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, s.Rows, s.Cols, s.Cols, 1, yinv, sqy0, 0, s); err != nil {
		panic(err)
	}

	return s
}

// z = sqrtz0 * (Id - s)**-1 * (Id + s) * sqrtz0
func StoZ(s *mat.CMatrix, z0 *mat.CVector) *mat.CMatrix {
	var err error

	z := s.DeepCopy()
	zinv := s.DeepCopy()
	ztmp := s.DeepCopy()

	// Generate sqrt(Z0) matrix and identity matrix for calculation
	sqz0 := cmf(z0.Size, z0.Size, opts)
	for i, val := range z0.Data {
		sqz0.Set(i, i, cmplx.Sqrt(val))
	}
	id := cmf(s.Rows, s.Cols, s.Opts)
	golapack.Zlaset(mat.Full, s.Rows, s.Cols, 0, 1, id)

	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, s.Rows, s.Cols, s.Cols, 1, id, id, 1, z); err != nil {
		panic(err)
	}
	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, s.Rows, s.Cols, s.Cols, 1, id, id, -1, zinv); err != nil {
		panic(err)
	}

	cmatInv(zinv)

	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, s.Rows, s.Cols, s.Cols, 1, sqz0, zinv, 0, ztmp); err != nil {
		panic(err)
	}
	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, s.Rows, s.Cols, s.Cols, 1, ztmp, z, 0, zinv); err != nil {
		panic(err)
	}
	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, s.Rows, s.Cols, s.Cols, 1, zinv, sqz0, 0, s); err != nil {
		panic(err)
	}

	return s
}

func TtoA(m *mat.CMatrix, z0 *mat.CVector) *mat.CMatrix {
	m = StoA(TtoS(m, z0), z0)
	return m
}

func TtoH(m *mat.CMatrix, z0 *mat.CVector) *mat.CMatrix {
	var err error

	if (m.Rows != m.Cols) && (m.Rows != 2) {
		err = fmt.Errorf("Network must only have 2 ports\n")
		panic(err)
	}

	m = StoH(TtoS(m, z0), z0)

	return m
}

func TtoS(m *mat.CMatrix, z0 *mat.CVector) *mat.CMatrix {
	var err error

	if (m.Rows != m.Cols) && ((m.Rows % 2) != 0) {
		err = fmt.Errorf("Network does not have an even number of ports\n")
		panic(err)
	}

	rh, ch := m.Rows/2, m.Cols/2
	minv := cmf(rh, ch, m.Opts)
	mtmp := cmf(rh, ch, m.Opts)
	s := m.DeepCopy()
	// Copy S matrix into minv for calculations
	golapack.Zlacpy(mat.Full, rh, ch, m.Off(rh, ch), minv)

	cmatInv(minv)

	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, rh, ch, ch, 1, m.Off(0, ch), minv, 0, s); err != nil {
		panic(err)
	}

	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, rh, ch, ch, 1, minv, m.Off(rh, 0), 0, mtmp); err != nil {
		panic(err)
	}
	golapack.Zlacpy(mat.Full, rh, ch, m, s.Off(0, ch))
	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, rh, ch, ch, -1, m.Off(0, ch), mtmp, 1, s.Off(0, ch)); err != nil {
		panic(err)
	}

	golapack.Zlacpy(mat.Full, rh, ch, minv, s.Off(rh, 0))

	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, rh, ch, ch, -1, minv, m.Off(rh, 0), 0, s.Off(rh, ch)); err != nil {
		panic(err)
	}

	m.Data = s.Data

	return m
}

func TtoY(m *mat.CMatrix, z0 *mat.CVector) *mat.CMatrix {
	m = StoY(TtoS(m, z0), z0)

	return m
}

func TtoZ(m *mat.CMatrix, z0 *mat.CVector) *mat.CMatrix {
	m = StoZ(TtoS(m, z0), z0)

	return m
}

func YtoA(m *mat.CMatrix) *mat.CMatrix {
	var err error

	if (m.Rows != m.Cols) && (m.Rows != 2) {
		err = fmt.Errorf("Network must only have 2 ports\n")
		panic(err)
	}

	a := m.DeepCopy()

	a.Set(0, 0, -m.Get(1, 1)/m.Get(1, 0))
	a.Set(0, 1, -1/m.Get(1, 0))
	a.Set(1, 0, (m.Get(0, 1)*m.Get(1, 0)-m.Get(0, 0)*m.Get(1, 1))/m.Get(1, 0))
	a.Set(1, 1, -m.Get(0, 0)/m.Get(1, 0))

	m.Data = a.Data

	return m
}

func YtoH(m *mat.CMatrix) *mat.CMatrix {
	var err error

	if (m.Rows != m.Cols) && (m.Rows != 2) {
		err = fmt.Errorf("Network must only have 2 ports\n")
		panic(err)
	}

	h := m.DeepCopy()

	h.Set(0, 0, 1/m.Get(0, 0))
	h.Set(0, 1, -m.Get(0, 1)/m.Get(0, 0))
	h.Set(1, 0, m.Get(1, 0)/m.Get(0, 0))
	h.Set(1, 1, (m.Get(0, 0)*m.Get(1, 1)-m.Get(0, 1)*m.Get(1, 0))/m.Get(0, 0))

	m.Data = h.Data

	return m
}

// s = (sqrty0 * z * sqrty0 - I) * (sqrty0 * z * sqrty0 + I)**-1
func YtoS(m *mat.CMatrix, z0 *mat.CVector) *mat.CMatrix {
	var err error

	zinv := m.DeepCopy()
	ztmp := m.DeepCopy()
	zwrk := m.DeepCopy()

	// Generate sqrt(Z0) matrix and identity matrix for calculation
	sqz0 := cmf(z0.Size, z0.Size, opts)
	for i, val := range z0.Data {
		sqz0.Set(i, i, cmplx.Sqrt(val))
	}
	golapack.Zlaset(mat.Full, m.Rows, m.Cols, 0, 1, ztmp)
	golapack.Zlaset(mat.Full, m.Rows, m.Cols, 0, 1, zinv)

	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, m.Rows, m.Cols, m.Cols, 1, sqz0, m, 0, zwrk); err != nil {
		panic(err)
	}
	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, m.Rows, m.Cols, m.Cols, 1, zwrk, sqz0, -1, ztmp); err != nil {
		panic(err)
	}

	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, m.Rows, m.Cols, m.Cols, 1, sqz0, m, 0, zwrk); err != nil {
		panic(err)
	}
	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, m.Rows, m.Cols, m.Cols, 1, zwrk, sqz0, 1, zinv); err != nil {
		panic(err)
	}

	cmatInv(zinv)

	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, m.Rows, m.Cols, m.Cols, 1, ztmp, zinv, 0, m); err != nil {
		panic(err)
	}

	return m
}

func YtoT(m *mat.CMatrix, z0 *mat.CVector) *mat.CMatrix {
	var err error

	if (m.Rows != m.Cols) && ((m.Rows % 2) != 0) {
		err = fmt.Errorf("Network does not have an even number of ports\n")
		panic(err)
	}

	m = StoT(YtoS(m, z0), z0)
	return m
}

func YtoZ(m *mat.CMatrix) *mat.CMatrix {
	return cmatInv(m)
}

func ZtoA(z *mat.CMatrix) *mat.CMatrix {
	var err error

	if (z.Rows != z.Cols) && (z.Rows != 2) {
		err = fmt.Errorf("Network must only have 2 ports\n")
		panic(err)
	}

	a := z.DeepCopy()

	a.Set(0, 0, z.Get(0, 0)/z.Get(1, 0))
	a.Set(0, 1, (z.Get(0, 0)*z.Get(1, 1)-z.Get(0, 1)*z.Get(1, 0))/z.Get(1, 0))
	a.Set(1, 0, 1/z.Get(1, 0))
	a.Set(1, 1, z.Get(1, 1)/z.Get(1, 0))

	z.Data = a.Data

	return z
}

func ZtoH(z *mat.CMatrix) *mat.CMatrix {
	var err error

	if (z.Rows != z.Cols) && (z.Rows != 2) {
		err = fmt.Errorf("Network must only have 2 ports\n")
		panic(err)
	}

	h := z.DeepCopy()

	h.Set(0, 0, (z.Get(0, 0)*z.Get(1, 1)-z.Get(1, 0)*z.Get(0, 1))/z.Get(1, 1))
	h.Set(0, 1, z.Get(0, 1)/z.Get(1, 1))
	h.Set(1, 0, -z.Get(1, 0)/z.Get(1, 1))
	h.Set(1, 1, 1/z.Get(1, 1))

	z.Data = h.Data

	return z
}

// s = (sqrty0 * z * sqrty0 - I) * (sqrty0 * z * sqrty0 + I)**-1
func ZtoS(z *mat.CMatrix, z0 *mat.CVector) *mat.CMatrix {
	var err error

	zinv := z.DeepCopy()
	ztmp := z.DeepCopy()
	zwrk := z.DeepCopy()

	// Generate sqrt(Z0) matrix and identity matrix for calculation
	sqy0 := cmf(z0.Size, z0.Size, opts)
	for i, val := range z0.Data {
		sqy0.Set(i, i, cmplx.Sqrt(1/val))
	}
	golapack.Zlaset(mat.Full, z.Rows, z.Cols, 0, 1, ztmp)
	golapack.Zlaset(mat.Full, z.Rows, z.Cols, 0, 1, zinv)

	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, z.Rows, z.Cols, z.Cols, 1, sqy0, z, 0, zwrk); err != nil {
		panic(err)
	}
	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, z.Rows, z.Cols, z.Cols, 1, zwrk, sqy0, -1, ztmp); err != nil {
		panic(err)
	}

	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, z.Rows, z.Cols, z.Cols, 1, sqy0, z, 0, zwrk); err != nil {
		panic(err)
	}
	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, z.Rows, z.Cols, z.Cols, 1, zwrk, sqy0, 1, zinv); err != nil {
		panic(err)
	}

	cmatInv(zinv)

	if err = goblas.Zgemm(mat.NoTrans, mat.NoTrans, z.Rows, z.Cols, z.Cols, 1, ztmp, zinv, 0, z); err != nil {
		panic(err)
	}

	return z
}

func ZtoT(z *mat.CMatrix, z0 *mat.CVector) *mat.CMatrix {
	var err error

	if (z.Rows != z.Cols) && ((z.Rows % 2) != 0) {
		err = fmt.Errorf("Network does not have an even number of ports\n")
		panic(err)
	}

	z = StoT(ZtoS(z, z0), z0)
	return z
}

func ZtoY(z *mat.CMatrix) *mat.CMatrix {
	return cmatInv(z)
}

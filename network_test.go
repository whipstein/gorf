package gorf

import (
	"math"
	"strconv"
	"testing"

	"github.com/whipstein/golinalg/golapack"
)

const (
	// eps = 2.2204460492503131e-016
	eps = 1e-09
)

func TestReadTouchstone(t *testing.T) {
	for _i, val := range []string{"./data/delay_short.s1p", "./data/line.s2p", "./data/tee.s3p", "./data/hfss_threeport_DB_50Ohm.s3p", "./data/ntwk_noise.s2p"} {
		net := NewNetwork()
		net.ReadTouchstone(val)

		for _j := range net.Data {
			if net.Freq.FreqScaled.Get(_j) != freq[_i][_j] {
				t.Errorf("Frequency doesn't match: got %v want %v\n", net.Freq.FreqScaled.Get(_j), freq[_i][_j])
			}
			for i := 0; i < net.NPorts; i++ {
				for j := 0; j < net.NPorts; j++ {
					if net.Data[_j].Get(i, j) != complex(convertNumber(s[_i][_j][i][j*2], s[_i][_j][i][j*2+1], enc[_i], RI)) {
						t.Errorf("Data doesn't match: got %v want %v\n", net.Data[_i].Get(i, j), complex(convertNumber(s[_i][_j][i][j*2], s[_i][_j][i][j*2+1], enc[_i], RI)))
					}
				}
			}
		}
	}
}

func TestConvertNumber(t *testing.T) {
	for i := range smatma {
		for j := 0; j < len(smatma[i])/2; j++ {
			val1, val2 := convertNumber(smatma[i][j*2], smatma[i][j*2+1], MA, MA)
			if math.Abs(val1-smatma[i][j*2]) > eps || math.Abs(val2-smatma[i][j*2+1]) > eps {
				t.Errorf("Value not correct: got %v, %v want %v, %v", val1, val2, smatma[i][j*2], smatma[i][j*2+1])
			}
			val1, val2 = convertNumber(smatma[i][j*2], smatma[i][j*2+1], MA, DB)
			if math.Abs(val1-smatdb[i][j*2]) > eps || math.Abs(val2-smatdb[i][j*2+1]) > eps {
				t.Errorf("Value not correct: got %v, %v want %v, %v", val1, val2, smatdb[i][j*2], smatdb[i][j*2+1])
			}
			val1, val2 = convertNumber(smatma[i][j*2], smatma[i][j*2+1], MA, RI)
			if math.Abs(val1-smatri[i][j*2]) > eps || math.Abs(val2-smatri[i][j*2+1]) > eps {
				t.Errorf("Value not correct: got %v, %v want %v, %v", val1, val2, smatri[i][j*2], smatri[i][j*2+1])
			}
			val1, val2 = convertNumber(smatri[i][j*2], smatri[i][j*2+1], RI, DB)
			if math.Abs(val1-smatdb[i][j*2]) > eps || math.Abs(val2-smatdb[i][j*2+1]) > eps {
				t.Errorf("Value not correct: got %v, %v want %v, %v", val1, val2, smatdb[i][j*2], smatdb[i][j*2+1])
			}
			val1, val2 = convertNumber(smatri[i][j*2], smatri[i][j*2+1], RI, MA)
			if math.Abs(val1-smatma[i][j*2]) > eps || math.Abs(val2-smatma[i][j*2+1]) > eps {
				t.Errorf("Value not correct: got %v, %v want %v, %v", val1, val2, smatma[i][j*2], smatma[i][j*2+1])
			}
		}
	}
}

func TestInvertMatrix(t *testing.T) {
	var info int
	var err error
	s := cmf(2, 2, opts)
	for i := 0; i < 2; i++ {
		for j := 0; j < 2; j++ {
			s.Set(i, j, sri2port_sto[i][j])
		}
	}
	lwork := s.Rows
	work := cvf(lwork)
	ipiv := make([]int, s.Rows)

	s.ToColMajor()
	if info, err = golapack.Zgetrf(s.Rows, s.Cols, s, &ipiv); err != nil || info != 0 {
		panic("golapack.Zgetrf error: " + strconv.Itoa(info))
	}
	if info, err = golapack.Zgetri(s.Rows, s, &ipiv, work, lwork); err != nil || info != 0 {
		panic("golapack.Zgetri error: " + strconv.Itoa(info))
	}
}

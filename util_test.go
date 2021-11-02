package gorf

import (
	"math"
	"testing"
)

func TestPassivity(t *testing.T) {
	var dataList = [][][]complex128{
		{
			{-7 + 8.2i, -1.7 - 2.9i},
			{-2 - 4i, -5 + 5.2i},
		},
		{
			{-2.9 + 4i, -8.8 + 4.3i, 2.7 - 5.9i, -2.7 + 6.3i, -9.2 + 8.5i, 3.6 - 0.8i},
			{-5.5 + 5.4i, -8.2 - 7.9i, -0.1 + 4.2i, 7.1 - 0.3i, 5.6 + 2.6i, 2.4 + 1.4i},
			{-5.2 - 1.2i, 7.3 + 5.6i, -4.9 + 2.2i, 8.2 - 2.7i, -5.4 + 4.6i, -6.1 + 3.1i},
			{-10.1 + 1.1i, -3.1 + 7.8i, 6.8 - 3i, -9 + 9.1i, -8.1 - 7.6i, 5.8 - 3i},
			{-8.2 + 6.4i, 4.5 + 4.3i, -9.7 + 4.3i, 0.03 - 3.1i, -8.2 + 7.9i, -9.5 + 3i},
			{-2.1 - 4.8i, 7.1 + 2.9i, -6.4 + 4.8i, 1.5 - 8.4i, -0.9 + 0.9i, -4.3 + 9.4i},
		},
	}
	var res = []float64{-141.89139789, -1254.587279543748}
	for idx, data := range dataList {
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := Passivity(mat)

		if math.Abs(dataOut-res[idx]) > eps {
			t.Errorf("result does not match: got %v want %v", dataOut, res[idx])
		}
	}
}

func BenchmarkPassivity(b *testing.B) {
	var data = [][]complex128{
		// {-7 + 8.2i, -1.7 - 2.9i},
		// {-2 - 4i, -5 + 5.2i},
		{-2.9 + 4i, -8.8 + 4.3i, 2.7 - 5.9i, -2.7 + 6.3i, -9.2 + 8.5i, 3.6 - 0.8i},
		{-5.5 + 5.4i, -8.2 - 7.9i, -0.1 + 4.2i, 7.1 - 0.3i, 5.6 + 2.6i, 2.4 + 1.4i},
		{-5.2 - 1.2i, 7.3 + 5.6i, -4.9 + 2.2i, 8.2 - 2.7i, -5.4 + 4.6i, -6.1 + 3.1i},
		{-10.1 + 1.1i, -3.1 + 7.8i, 6.8 - 3i, -9 + 9.1i, -8.1 - 7.6i, 5.8 - 3i},
		{-8.2 + 6.4i, 4.5 + 4.3i, -9.7 + 4.3i, 0.03 - 3.1i, -8.2 + 7.9i, -9.5 + 3i},
		{-2.1 - 4.8i, 7.1 + 2.9i, -6.4 + 4.8i, 1.5 - 8.4i, -0.9 + 0.9i, -4.3 + 9.4i},
	}
	size := len(data)
	z0 := cvf(size)
	z0.SetReAll(50)
	mat := cmf(size, size)
	for i := 0; i < size; i++ {
		for j := 0; j < size; j++ {
			mat.Set(i, j, data[i][j])
		}
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		dataOut := Passivity(mat)
		_ = dataOut
	}
}

func TestAtoH(t *testing.T) {
	for _, data := range [][][]complex128{ari2port_ato} {
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}
		outMat := mat.DeepCopy()
		outMat = StoH(AtoS(outMat, z0), z0)

		dataOut := AtoH(mat)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outMat.Get(i, j))) > eps || math.Abs(mat.GetIm(i, j)-imag(outMat.Get(i, j))) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outMat.Get(i, j))
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outMat.Get(i, j))) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outMat.Get(i, j))) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outMat.Get(i, j))
				}
			}
		}
	}
}
func TestAtoS(t *testing.T) {
	for idx, data := range [][][]complex128{ari2port_ato} {
		outData := [][][]complex128{sri2port_ato}
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := AtoS(mat, z0)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}
func TestAtoT(t *testing.T) {
	for idx, data := range [][][]complex128{ari2port_ato} {
		outData := [][][]complex128{tri2port_ato}
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := AtoT(mat, z0)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}
func TestAtoY(t *testing.T) {
	for _, data := range [][][]complex128{ari2port_ato} {
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}
		outMat := mat.DeepCopy()
		outMat = StoY(AtoS(outMat, z0), z0)

		dataOut := AtoY(mat)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outMat.Get(i, j))) > eps || math.Abs(mat.GetIm(i, j)-imag(outMat.Get(i, j))) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outMat.Get(i, j))
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outMat.Get(i, j))) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outMat.Get(i, j))) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outMat.Get(i, j))
				}
			}
		}
	}
}
func TestAtoZ(t *testing.T) {
	for _, data := range [][][]complex128{ari2port_ato} {
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}
		outMat := mat.DeepCopy()
		outMat = StoZ(AtoS(outMat, z0), z0)

		dataOut := AtoZ(mat)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outMat.Get(i, j))) > eps || math.Abs(mat.GetIm(i, j)-imag(outMat.Get(i, j))) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outMat.Get(i, j))
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outMat.Get(i, j))) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outMat.Get(i, j))) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outMat.Get(i, j))
				}
			}
		}
	}
}

func TestHtoA(t *testing.T) {
	for _, data := range [][][]complex128{hri2port_hto} {
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}
		outMat := mat.DeepCopy()
		outMat = StoA(HtoS(outMat, z0), z0)

		dataOut := HtoA(mat)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outMat.Get(i, j))) > eps || math.Abs(mat.GetIm(i, j)-imag(outMat.Get(i, j))) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outMat.Get(i, j))
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outMat.Get(i, j))) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outMat.Get(i, j))) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outMat.Get(i, j))
				}
			}
		}
	}
}
func TestHtoS(t *testing.T) {
	for idx, data := range [][][]complex128{hri2port_hto} {
		outData := [][][]complex128{sri2port_hto}
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := HtoS(mat, z0)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}
func TestHtoT(t *testing.T) {
	for idx, data := range [][][]complex128{hri2port_hto} {
		outData := [][][]complex128{tri2port_hto}
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := HtoT(mat, z0)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}
func TestHtoY(t *testing.T) {
	for _, data := range [][][]complex128{hri2port_hto} {
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}
		outMat := mat.DeepCopy()
		outMat = StoY(HtoS(outMat, z0), z0)

		dataOut := HtoY(mat)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outMat.Get(i, j))) > eps || math.Abs(mat.GetIm(i, j)-imag(outMat.Get(i, j))) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outMat.Get(i, j))
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outMat.Get(i, j))) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outMat.Get(i, j))) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outMat.Get(i, j))
				}
			}
		}
	}
}
func TestHtoZ(t *testing.T) {
	for _, data := range [][][]complex128{hri2port_hto} {
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}
		outMat := mat.DeepCopy()
		outMat = StoZ(HtoS(outMat, z0), z0)

		dataOut := HtoZ(mat)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outMat.Get(i, j))) > eps || math.Abs(mat.GetIm(i, j)-imag(outMat.Get(i, j))) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outMat.Get(i, j))
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outMat.Get(i, j))) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outMat.Get(i, j))) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outMat.Get(i, j))
				}
			}
		}
	}
}

func TestStoA(t *testing.T) {
	for idx, data := range [][][]complex128{sri2port_sto} {
		outData := [][][]complex128{ari2port_sto}
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := StoA(mat, z0)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}
func TestStoH(t *testing.T) {
	for idx, data := range [][][]complex128{sri2port_sto} {
		outData := [][][]complex128{hri2port_sto}
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := StoH(mat, z0)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}
func TestStoT(t *testing.T) {
	for idx, data := range [][][]complex128{sri2port_sto, sri_sto} {
		outData := [][][]complex128{tri2port_sto, tri_sto}
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := StoT(mat, z0)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}
func TestStoY(t *testing.T) {
	for idx, data := range [][][]complex128{sri2port_sto, sri_sto} {
		outData := [][][]complex128{yri2port_sto, yri_sto}
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := StoY(mat, z0)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}
func TestStoZ(t *testing.T) {
	for idx, data := range [][][]complex128{sri2port_sto, sri_sto} {
		outData := [][][]complex128{zri2port_sto, zri_sto}
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := StoZ(mat, z0)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}

func TestTtoA(t *testing.T) {
	for idx, data := range [][][]complex128{tri2port_tto} {
		outData := [][][]complex128{ari2port_tto}
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := TtoA(mat, z0)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}
func TestTtoH(t *testing.T) {
	for idx, data := range [][][]complex128{tri2port_tto} {
		outData := [][][]complex128{hri2port_tto}
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := TtoH(mat, z0)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}
func TestTtoS(t *testing.T) {
	for idx, data := range [][][]complex128{tri2port_tto, tri_tto} {
		outData := [][][]complex128{sri2port_tto, sri_tto}
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := TtoS(mat, z0)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}
func TestTtoY(t *testing.T) {
	for idx, data := range [][][]complex128{tri2port_tto, tri_tto} {
		outData := [][][]complex128{yri2port_tto, yri_tto}
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := TtoY(mat, z0)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}
func TestToZ(t *testing.T) {
	for idx, data := range [][][]complex128{tri2port_tto, tri_tto} {
		outData := [][][]complex128{zri2port_tto, zri_tto}
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := TtoZ(mat, z0)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}

func TestYtoA(t *testing.T) {
	for idx, data := range [][][]complex128{yri2port_yto} {
		outData := [][][]complex128{ari2port_yto}
		size := len(data)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := YtoA(mat)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}
func TestYtoH(t *testing.T) {
	for idx, data := range [][][]complex128{yri2port_yto} {
		outData := [][][]complex128{hri2port_yto}
		size := len(data)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := YtoH(mat)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}
func TestYtoS(t *testing.T) {
	for idx, data := range [][][]complex128{yri2port_yto, yri_yto} {
		outData := [][][]complex128{sri2port_yto, sri_yto}
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := YtoS(mat, z0)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}
func TestYtoT(t *testing.T) {
	for idx, data := range [][][]complex128{yri2port_yto, yri_yto} {
		outData := [][][]complex128{tri2port_yto, tri_yto}
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := YtoT(mat, z0)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}
func TestYtoZ(t *testing.T) {
	for idx, data := range [][][]complex128{yri2port_yto, yri_yto} {
		outData := [][][]complex128{zri2port_yto, zri_yto}
		size := len(data)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := YtoZ(mat)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}

func TestZtoA(t *testing.T) {
	for idx, data := range [][][]complex128{zri2port_zto} {
		outData := [][][]complex128{ari2port_zto}
		size := len(data)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := ZtoA(mat)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}
func TestZtoH(t *testing.T) {
	for idx, data := range [][][]complex128{zri2port_zto} {
		outData := [][][]complex128{hri2port_zto}
		size := len(data)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := ZtoH(mat)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}
func TestZtoS(t *testing.T) {
	for idx, data := range [][][]complex128{zri2port_zto, zri_zto} {
		outData := [][][]complex128{sri2port_zto, sri_zto}
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := ZtoS(mat, z0)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}
func TestZtoT(t *testing.T) {
	for idx, data := range [][][]complex128{zri2port_zto, zri_zto} {
		outData := [][][]complex128{tri2port_zto, tri_zto}
		size := len(data)
		z0 := cvf(size)
		z0.SetReAll(50)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := ZtoT(mat, z0)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}
func TestZtoY(t *testing.T) {
	for idx, data := range [][][]complex128{zri2port_zto, zri_zto} {
		outData := [][][]complex128{yri2port_zto, yri_zto}
		size := len(data)
		mat := cmf(size, size)
		for i := 0; i < size; i++ {
			for j := 0; j < size; j++ {
				mat.Set(i, j, data[i][j])
			}
		}

		dataOut := ZtoY(mat)

		for i := 0; i < mat.Rows; i++ {
			for j := 0; j < mat.Cols; j++ {
				if math.Abs(mat.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(mat.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, mat.Get(i, j), outData[idx][i][j])
				}
				if math.Abs(dataOut.GetRe(i, j)-real(outData[idx][i][j])) > eps || math.Abs(dataOut.GetIm(i, j)-imag(outData[idx][i][j])) > eps {
					t.Errorf("row %d col %d does not match: got %v want %v", i, j, dataOut.Get(i, j), outData[idx][i][j])
				}
			}
		}
	}
}

var enc = []Encoding{RI, RI, RI, DB, RI}
var freq = [][]float64{
	{75.0, 75.175, 75.35, 75.525, 75.7},
	{75.0, 75.175, 75.35, 75.525, 75.7, 75.875, 76.05, 76.225, 76.4, 76.575, 76.75, 76.925},
	{330.0, 330.85, 331.7, 332.55, 333.4, 334.25, 335.1, 335.95, 336.8},
	{2.9, 2.91022222222222, 2.92044444444444, 2.93066666666667, 2.94088888888889, 2.95111111111111},
	{1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0},
}
var s = [][][][]float64{
	{
		{
			{0.45345337996, 0.891279996524},
		},
		{
			{0.464543921496, 0.885550080459},
		},
		{
			{0.475521092896, 0.879704319764},
		},
		{
			{0.486384743723, 0.873744745949},
		},
		{
			{0.497134729488, 0.867673360624},
		},
	},
	{
		{
			{0.0, 0.0, 0.52275549736, -0.852482662568},
			{0.52275549736, -0.852482662568, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.517424428542, -0.855728906107},
			{0.517424428542, -0.855728906107, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.512093207875, -0.858929884477},
			{0.512093207875, -0.858929884477, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.506761904782, -0.862086058269},
			{0.506761904782, -0.862086058269, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.501430588672, -0.865197876063},
			{0.501430588672, -0.865197876063, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.496099328941, -0.868265774878},
			{0.496099328941, -0.868265774878, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.49076819497, -0.871290180598},
			{0.49076819497, -0.871290180598, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.485437256126, -0.87427150838},
			{0.485437256126, -0.87427150838, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.480106581757, -0.877210163048},
			{0.480106581757, -0.877210163048, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.474776241195, -0.880106539458},
			{0.474776241195, -0.880106539458, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.469446303753, -0.882961022862},
			{0.469446303753, -0.882961022862, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.464116838725, -0.885773989239},
			{0.464116838725, -0.885773989239, 0.0, 0.0},
		},
	},
	{
		{
			{-0.333333333333, 0.0, 0.666666666667, 0.0, 0.666666666667, 0.0},
			{0.666666666667, 0.0, -0.333333333333, 0.0, 0.666666666667, 0.0},
			{0.666666666667, 0.0, 0.666666666667, 0.0, -0.333333333333, 0.0},
		},
		{
			{-0.333333333333, 0.0, 0.666666666667, 0.0, 0.666666666667, 0.0},
			{0.666666666667, 0.0, -0.333333333333, 0.0, 0.666666666667, 0.0},
			{0.666666666667, 0.0, 0.666666666667, 0.0, -0.333333333333, 0.0},
		},
		{
			{-0.333333333333, 0.0, 0.666666666667, 0.0, 0.666666666667, 0.0},
			{0.666666666667, 0.0, -0.333333333333, 0.0, 0.666666666667, 0.0},
			{0.666666666667, 0.0, 0.666666666667, 0.0, -0.333333333333, 0.0},
		},
		{
			{-0.333333333333, 0.0, 0.666666666667, 0.0, 0.666666666667, 0.0},
			{0.666666666667, 0.0, -0.333333333333, 0.0, 0.666666666667, 0.0},
			{0.666666666667, 0.0, 0.666666666667, 0.0, -0.333333333333, 0.0},
		},
		{
			{-0.333333333333, 0.0, 0.666666666667, 0.0, 0.666666666667, 0.0},
			{0.666666666667, 0.0, -0.333333333333, 0.0, 0.666666666667, 0.0},
			{0.666666666667, 0.0, 0.666666666667, 0.0, -0.333333333333, 0.0},
		},
		{
			{-0.333333333333, 0.0, 0.666666666667, 0.0, 0.666666666667, 0.0},
			{0.666666666667, 0.0, -0.333333333333, 0.0, 0.666666666667, 0.0},
			{0.666666666667, 0.0, 0.666666666667, 0.0, -0.333333333333, 0.0},
		},
		{
			{-0.333333333333, 0.0, 0.666666666667, 0.0, 0.666666666667, 0.0},
			{0.666666666667, 0.0, -0.333333333333, 0.0, 0.666666666667, 0.0},
			{0.666666666667, 0.0, 0.666666666667, 0.0, -0.333333333333, 0.0},
		},
		{
			{-0.333333333333, 0.0, 0.666666666667, 0.0, 0.666666666667, 0.0},
			{0.666666666667, 0.0, -0.333333333333, 0.0, 0.666666666667, 0.0},
			{0.666666666667, 0.0, 0.666666666667, 0.0, -0.333333333333, 0.0},
		},
		{
			{-0.333333333333, 0.0, 0.666666666667, 0.0, 0.666666666667, 0.0},
			{0.666666666667, 0.0, -0.333333333333, 0.0, 0.666666666667, 0.0},
			{0.666666666667, 0.0, 0.666666666667, 0.0, -0.333333333333, 0.0},
		},
	},
	{
		{
			{-0.21148658816906, -8.66528404051443, -13.4849497033415, 83.5500112591331, -25.6759232493857, 64.8275204118573},
			{-13.4849497033415, 83.5500112591332, -0.377315745326041, -1.59879806058515, -14.1578618605019, -89.1053094741142},
			{-25.6759232493856, 64.8275204118573, -14.1578618605018, -89.1053094741143, -0.18224923110063, 6.45587163577479},
		},
		{
			{-0.23006525822869, -9.64943298495914, -13.1862404556219, 82.3960531750091, -24.4564518184145, 66.4298489263103},
			{-13.1862404556165, 82.3960531749996, -0.39422059998425, -2.55268523426717, -14.1162236359386, -90.2179228572975},
			{-24.4564518184146, 66.4298489263728, -14.1162236359403, -90.2179228572859, -0.187891831489039, 5.82307554364623},
		},
		{
			{-0.253309235365382, -10.7176894894606, -12.841308967834, 81.1741709223504, -23.2981060086464, 67.3616144060196},
			{-12.8413089678254, 81.1741709223365, -0.416829901061395, -3.54891280505034, -14.030500776937, -91.3896175020961},
			{-23.2981060086453, 67.3616144061121, -14.0305007769396, -91.3896175020766, -0.196365672889594, 5.19348463454599},
		},
		{
			{-0.282386367682559, -11.8882390557008, -12.4481635094292, 79.8640837927267, -22.1822990969097, 67.7732137468462},
			{-12.4481635094188, 79.864083792709, -0.446285748139503, -4.60189715702798, -13.8989889946607, -92.6404751217517},
			{-22.1822990969072, 67.7732137469483, -13.8989889946638, -92.6404751217265, -0.208129123449943, 4.56091373668791},
		},
		{
			{-0.318920930224311, -13.183771525354, -12.0042222895624, 78.4413289944383, -21.0925863507959, 67.7537609132864},
			{-12.0042222895515, 78.4413289944181, -0.48421030224758, -5.72863722888568, -13.7193992557416, -93.9947067802764},
			{-21.0925863507924, 67.7537609133809, -13.7193992557449, -93.9947067802486, -0.223855527820806, 3.91886378103244},
		},
		{
			{-0.365199591126864, -14.6330818460627, -11.5062559584243, 76.8755997303085, -20.0144038963704, 67.3508752903833},
			{-11.5062559584136, 76.8755997302882, -0.53290786979867, -6.94969804518525, -13.4888020125676, -95.4823116155129},
			{-20.0144038963663, 67.3508752904675, -13.4888020125706, -95.4823116154839, -0.244513329310037, 3.26031624998304},
		},
	},
	{
		{
			{0.0, 0.0, 0.0, 0.0},
			{10.0, 0.0, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.0, 0.0},
			{10.0, 0.0, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.0, 0.0},
			{10.0, 0.0, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.0, 0.0},
			{10.0, 0.0, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.0, 0.0},
			{10.0, 0.0, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.0, 0.0},
			{10.0, 0.0, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.0, 0.0},
			{10.0, 0.0, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.0, 0.0},
			{10.0, 0.0, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.0, 0.0},
			{10.0, 0.0, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.0, 0.0},
			{10.0, 0.0, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.0, 0.0},
			{10.0, 0.0, 0.0, 0.0},
		},
		{
			{0.0, 0.0, 0.0, 0.0},
			{10.0, 0.0, 0.0, 0.0},
		},
	},
}
var nfreq = [][]float64{
	{},
	{},
	{},
	{},
	{1, 2},
}
var nfmin = [][]float64{
	{},
	{},
	{},
	{},
	{0.5000, 1.0000},
}
var gopt = [][][]float64{
	{{}},
	{{}},
	{{}},
	{{}},
	{
		{0.00000, 134.27},
		{0.00000, 134.27},
	},
}
var rn = [][]float64{
	{},
	{},
	{},
	{},
	{0.1159, 0.1159},
}

var sfreq = []float64{30}

var ari2port_ato = [][]complex128{
	{1 + 0i, 1 - 2i},
	{8 - 2i, 4 + 2i},
}
var sri2port_ato = [][]complex128{
	{-0.995196836229336 + 0.0009651428968344965i, -0.022582961719726166 + 0.09329407543577613i},
	{0.0046647037717888066 + 0.0011291480859863083i, -0.9834610210859421 + 0.013681994698371022i},
}
var tri2port_ato = [][]complex128{
	{-197.51 + 51.02i, -201.49 + 48.98i},
	{198.49 - 50.98i, 202.51 - 49.02i},
}

var hri2port_hto = [][]complex128{
	{1 + 0i, 1 - 2i},
	{8 - 2i, 4 + 2i},
}
var sri2port_hto = [][]complex128{
	{-0.9122765010131432 + 0.14661814805281648i, -0.001422393083897159 - 0.019049412147708392i},
	{-0.04992461 + 0.04970129i, -0.99251800803285 - 0.004466416456362586i},
}
var tri2port_hto = [][]complex128{
	{10.530588235294118 + 7.632647058823529i, 10.645882352941177 + 7.661470588235293i},
	{-9.94 - 9.985i, -10.06 - 10.015i},
}

var ari2port_sto = [][]complex128{
	{0.20588235294117646 - 1.3235294117647058i, 1.4705882352941178 + 69.11764705882354i},
	{-0.01 + 0.02i, 0.5 - 1i},
}
var hri2port_sto = [][]complex128{
	{-54.70588235294119 + 28.823529411764707i, 0.235294118 + 0.0588235294i},
	{-0.4 - 0.8i, -0.02 - 2.11774071e-15i},
}
var sri_sto = [][]complex128{
	{-0.8 + 0.6i, -0.02 - 0.02i, -0.9 - 0.3i, 0.8 + 1i, 0.6 + 0.8i, -0.7 - 0.1i},
	{-0.02 - 0.03i, -1.8 - 2.6i, -0.02 - 0.03i, 0.6 + 0.8i, 0.7 + 0.9i, 0.2 + 0.2i},
	{-0.9 - 0.3i, -0.02 - 0.4i, -0.8 + 1.6i, -0.7 - 0.1i, 0.2 + 0.2i, 0.9 + 0.01i},
	{0.3 + 0.01i, 0.6 + 0.8i, -0.7 - 0.1i, 0.8 + 0.6i, -0.02 - 0.03i, -0.9 - 0.3i},
	{0.6 + 0.9i, 0.7 + 0.9i, 0.2 + 0.2i, -0.02 - 0.03i, 0.8 + 0.2i, -0.02 - 0.03i},
	{-0.7 - 0.1i, 0.2 + 0.2i, 0.9 + 0.01i, -0.9 - 0.3i, -0.02 - 0.03i, -2.8 + 1.6i},
}
var sri2port_sto = [][]complex128{
	{1 + 0i, 1 - 2i},
	{8 - 2i, 4 + 2i},
}
var tri_sto = [][]complex128{
	{-0.025430942526377898 + 0.3559578338787033i, 1.1690862009698817 + 0.17083167472109734i, -2.1563357070276195 - 2.8666181097329857i, 0.5393411091734047 - 0.7340112035893411i, -0.4250271909370256 + 0.8627441573803569i, -0.22321864146330608 - 0.9390943760769238i},
	{0.8918570770961208 + 1.4861928710434937i, 1.0588483768289327 + 0.760036017393786i, -6.609589346663594 + 0.7411920320794738i, -2.165155062178381 - 0.2476282774321544i, -0.4827531867455503 + 0.14101468696824682i, -1.5442513306254184 - 0.3734056523984794i},
	{-2.17252390637847 + 1.3796429658338654i, 0.7080654864358843 - 0.6722669368523164i, 1.0610293860475477 + 0.07007689638725823i, 7.76508196e-02 - 1.62461958i, -0.31153080252294846 + 1.126975106715508i, -0.3244768950237382 + 0.3452080936287402i},
	{0.003774767030086945 - 0.007865286691803502i, -0.4741289384013674 + 0.3204540411340853i, -1.2717167255587263 + 2.1546319373935328i, -0.33039380429028875 + 0.47327965633887276i, 0.4118109067869533 - 0.5142577514230819i, -0.5113543596739077 + 0.3598436572579409i},
	{-2.42887803e-01 - 0.05210235i, -0.018693333967734167 + 0.11257850279495207i, 1.0860206059729007 - 1.8915733553943657i, 0.4640392996612892 - 0.5275968875034696i, 4.15702712e-02 - 0.14135739i, 0.375104228001509 - 0.3407861221743924i},
	{1.0504456962223538 + 0.3815170554573251i, -0.35061336287589534 + 0.2129267053522648i, 1.2200425723015869 - 0.07779910246569244i, -0.5260191340541548 + 0.35136454213017737i, 3.33056007e-01 - 0.33574719i, 5.16651619e-01 + 0.20969443i},
}
var tri2port_sto = [][]complex128{
	{0.5882352941176471 - 2.3529411764705883i, 0.11764705882352941 + 0.029411764705882353i},
	{-0.4117647058823529 - 0.3529411764705882i, 0.11764705882352941 + 0.029411764705882353i},
}
var yri_sto = [][]complex128{
	{-0.03588858852008332 - 0.0658123745950294i, 0.0022073852498142015 + 0.02784480264998655i, -0.03166723934031997 - 0.0015259953789566194i, -0.0228230914728124 + 0.03661362157846373i, -0.006233717891954322 + 0.021471110846780227i, -0.006070543915571594 + 0.008842931938977613i},
	{-0.004085862560114975 + 0.018936020929956904i, -0.019671852248620874 + 0.002965454692701138i, 0.00654032 + 0.00247292i, 0.0125212207801607 - 0.011390576073113464i, 0.010167374974193975 - 0.007899856994114225i, -0.0005848040438370365 - 0.0024913530501663624i},
	{-0.04406712 - 0.01344734i, 0.015016916111850106 + 0.008844761549918296i, -0.031136473748422617 - 0.00707107803997322i, 0.014371960972328594 + 0.014712114561606097i, 0.009320677915364077 + 0.014620180109104925i, 0.005541446948157051 + 0.001588138879700725i},
	{-0.013991324364792997 + 0.003273378041821533i, 0.007756157726080032 - 0.004128771792023683i, -0.00021814769190118223 - 0.003714094598007246i, 0.00309600569813617 - 0.007285733274469971i, 0.001445113994712625 + 0.0033527139005295023i, -0.002175969518379138 - 0.004152901283770281i},
	{-0.010139012024407774 + 0.03163125907453644i, 0.012352006388905238 - 0.015589025042379549i, 0.010410910804243799 + 0.013029077209978656i, 0.015047033729639567 - 0.007098406643217624i, 0.006428772939992055 - 0.011439879692648229i, 0.005033745807298377 - 0.0003226231811509888i},
	{-0.018476894727085635 + 0.005246748160273587i, 0.004267032320675186 - 0.0015474477876960783i, 0.003106402982285195 + 0.004093166010069202i, 0.007132209423435622 + 0.0005522089638828575i, 0.006639057845670114 + 0.003281064955360188i, -0.028843511685587027 - 0.008120703622191007i},
}
var yri2port_sto = [][]complex128{
	{-0.014307692307692304 - 0.007538461538461537i, 0.0029230769230769245 + 0.002615384615384618i},
	{-0.0003076923076923076 + 0.014461538461538465i, -0.019076923076923078 - 0.003384615384615384i},
}
var zri_sto = [][]complex128{
	{-23.866727610349816 + 4.031123596868419i, -18.798452504897472 + 4.333531158923205i, 2.731091168119487 - 5.894807692620406i, -32.76854156386609 + 6.3282697088151i, -29.255750578630277 + 18.599393857818317i, 3.645426608194651 - 0.8243755508170214i},
	{-15.518379302325775 + 5.426632060893081i, -28.246798870911118 - 7.9179642874874805i, -0.11122682186843467 + 4.266557176622344i, 7.154530384986955 - 0.39869546177792486i, 5.668079167004378 + 32.664578742932726i, 2.4071463130761805 + 1.4735786836384956i},
	{-5.221298321118842 - 11.243592394823935i, 7.263988178149836 + 5.589254026431631i, -34.97467568056724 + 22.297540925813443i, 38.25067396212568 - 2.770141811623657i, -5.470265383487061 + 24.629643789503554i, -6.136282102543251 + 3.1717689971842766i},
	{-10.10972426292159 + 1.146252737106027i, -3.14713115443701 + 7.868938025521188i, 26.766639817760456 - 13.04369222819634i, -39.04653162629229 + 69.1322920594311i, -8.135833050521395 - 27.665084522821267i, 15.873502089417238 - 13.060216595585494i},
	{-38.18573219333717 + 26.401174173837454i, 14.499327398265512 + 24.266019231504433i, -19.771606346222597 + 24.370056898728635i, 10.033520941677368 - 63.102435305358036i, -8.244315571552017 + 127.93335910023917i, -9.580983556545187 + 13.017731541048288i},
	{-2.0898189233901636 - 4.8288898219031635i, 7.104483504459287 + 2.9014099571191054i, -6.492431543803547 + 4.82356580467751i, 21.526831121137665 - 8.424610273659319i, -0.9833045586803587 + 10.957686492906507i, -34.391342818653236 + 9.414464519938905i},
}
var zri2port_sto = [][]complex128{
	{-57.058823529411775 + 18.235294117647054i, -11.764705882352946 - 2.9411764705882373i},
	{-20 - 40i, -50 + 5.29435177e-12i},
}

var ari2port_tto = [][]complex128{
	{-7.85e+00 + 3.250e+00i, 5.75e+01 - 4.750e+01i},
	{1.70e-02 - 4.100e-02i, -4.15e+00 + 1.015e+01i},
}
var hri2port_tto = [][]complex128{
	{-5.994012225040564 - 3.2142708636534048i, -7.61631668676452 + 3.0588881034554465i},
	{3.45128696e-02 + 8.44109942e-02i, -4.04756955e-03 - 1.99592499e-05i},
}
var sri_tto = [][]complex128{
	{1.61954120e-01 - 6.94626890e-01i, 4.33290059e-01 - 3.16633240e-03i, -1.46410606e-01 - 8.74652950e-01i, 5.395232882022029 - 8.532443533251088i, -17.17643072073374 + 5.669142446365555i, 2.7365289578372556 - 7.479543435676878i},
	{-3.63453122e-01 - 4.73636626e-01i, -6.82339504e-01 - 1.18511140e-01i, 3.86524559e-01 - 8.01313704e-02i, -14.84919989625736 + 6.098291994573872i, -13.4368458349315 - 3.6179676825911877i, -1.2467733836938382 + 5.746712844956601i},
	{-1.12925870e-01 - 5.96466395e-01i, -4.92727454e-02 + 1.80378419e-01i, 8.55191810e-01 - 1.00785310e-02i, -4.401999208479723 - 1.2218877570434303i, -2.8062431148629186 + 1.6234472065027337i, 3.3798273088959 + 3.70931405789427i},
	{-1.87895131e-02 - 7.12868950e-02i, -7.79833910e-02 + 2.21370411e-02i, 8.46467846e-02 - 1.31807911e-03i, -5.81891382e-01 + 3.84824864e-01i, -7.72985335e-01 - 7.48365902e-02i, 2.15792400e-01 + 5.63699954e-01i},
	{-1.98107105e-02 + 4.49897301e-02i, -1.23796523e-02 - 5.02163581e-02i, -5.52824379e-02 + 5.59614174e-02i, -9.58205237e-01 - 4.19303246e-03i, 6.84078207e-01 + 3.36190828e-01i, -4.21462125e-01 - 1.75719692e-01i},
	{-2.54916314e-02 - 5.74512011e-02i, -5.85021596e-02 + 4.00927398e-02i, 3.59303054e-02 - 1.10554583e-01i, 6.23296692e-02 + 9.12587929e-02i, -6.67198386e-01 + 7.72617614e-01i, -3.50083514e-01 + 7.46373367e-02i},
}
var sri2port_tto = [][]complex128{
	{-0.1264411990776326 + 0.4485011529592621i, -9.04688701 + 8.59123751i},
	{-0.09607993850883935 - 0.09992313604919292i, 0.207532667179093 - 0.5841660261337432i},
}
var tri_tto = [][]complex128{
	{-2.9 + 4i, -8.8 + 4.3i, 2.7 - 5.9i, -2.7 + 6.3i, -9.2 + 8.5i, 3.6 - 0.8i},
	{-5.5 + 5.4i, -8.2 - 7.9i, -0.1 + 4.2i, 7.1 - 0.3i, 5.6 + 2.6i, 2.4 + 1.4i},
	{-5.2 - 1.2i, 7.3 + 5.6i, -4.9 + 2.2i, 8.2 - 2.7i, -5.4 + 4.6i, -6.1 + 3.1i},
	{-10.1 + 1.1i, -3.1 + 7.8i, 6.8 - 3i, -9 + 9.1i, -8.1 - 7.6i, 5.8 - 3i},
	{-8.2 + 6.4i, 4.5 + 4.3i, -9.7 + 4.3i, 0.03 - 3.1i, -8.2 + 7.9i, -9.5 + 3i},
	{-2.1 - 4.8i, 7.1 + 2.9i, -6.4 + 4.8i, 1.5 - 8.4i, -0.9 + 0.9i, -4.3 + 9.4i},
}
var tri2port_tto = [][]complex128{
	{-7 + 8.2i, -1.7 - 2.9i},
	{-2 - 4i, -5 + 5.2i},
}
var yri_tto = [][]complex128{
	{-3.18858922e-02 - 5.87469140e-02i, -1.97970970e-02 - 2.00773656e-02i, -1.02101818e-03 + 3.71108316e-02i, -1.5963774250081555 - 0.4357630353018983i, -6.47105790e-01 - 9.14542103e-01i, 7.91332669e-01 + 7.30327908e-01i},
	{-5.93573089e-02 - 2.30862089e-02i, -3.79720174e-02 - 4.36057820e-03i, 3.15974042e-02 + 1.41291350e-02i, -1.3691699014945682 + 0.7036585542691122i, -1.19139846e+00 + 6.26266705e-02i, 6.24125622e-01 - 8.84551877e-02i},
	{4.14542628e-02 - 4.28609410e-03i, 6.75386385e-03 - 8.91643616e-03i, -1.33127252e-02 - 1.42890464e-02i, 5.54911709e-01 - 7.25173122e-01i, 5.64778805e-01 - 2.96357572e-01i, -6.45728221e-01 + 2.48594769e-01i},
	{-1.19486661e-03 - 8.15137349e-04i, -1.84730221e-03 - 3.75725763e-04i, -3.16516524e-05 - 4.14543806e-04i, -2.71758649e-02 + 8.93167758e-03i, -3.26323186e-02 + 2.95574812e-03i, -5.95074573e-03 + 4.66666483e-03i},
	{-4.64223307e-03 - 1.66635948e-03i, -2.61748640e-03 + 2.71005906e-04i, 2.04582978e-03 + 8.81420795e-04i, -9.18631246e-02 + 5.53211864e-02i, -8.80337163e-02 - 2.16089940e-03i, 6.48718191e-02 - 4.94172051e-03i},
	{-9.12991727e-03 + 1.06196887e-02i, -1.56240032e-03 + 3.68517946e-03i, 5.60524678e-03 + 4.88737540e-04i, -4.46240269e-03 + 2.90741078e-01i, -7.57565942e-02 + 1.90174399e-01i, 9.49919385e-02 - 1.66462138e-01i},
}
var yri2port_tto = [][]complex128{
	{-0.12957303370786522 + 0.06948314606741576i, -0.7743280898876407 + 0.9255550561797757i},
	{-0.010337078651685394 - 0.00853932584269663i, -0.10889887640449439 - 0.03343820224719102i},
}
var zri_tto = [][]complex128{
	{-29.372620106578253 + 25.55076953622686i, 21.77538527951088 + 1.535670652813163i, -81.35103258266886 - 184.62019985852467i, 41.57003498203392 + 615.6147700986596i, -113.75624999969479 + 20.43023837085962i, 276.6831713344474 - 517.430272123246i},
	{-26.302854501997285 - 51.2272285443978i, 39.22359250892711 - 22.62958953730651i, 48.92759562977502 - 8.427818137200424i, -763.8828805572114 + 246.82777265740037i, -956.0490183661258 + 143.34209654789245i, -104.95466390065363 + 350.29002043557193i},
	{-40.143701974361136 + 6.338688913397791i, 31.886693129926076 - 9.332201957057064i, 78.89829676520411 + 20.407632065763444i, -913.7071779793919 - 207.09160501229516i, 412.8696986000733 + 475.9806664526273i, 94.55769831410966 + 173.12617958909044i},
	{-3.141376686249023 - 1.133656335089378i, -2.493717674293195 - 0.22026198846243006i, 0.8602777538070367 + 3.9607129230365947i, 6.790963961978615 - 24.554517845029153i, 32.399198801651124 + 53.79223491161721i, -2.7814045995934347 + 15.597187266835212i},
	{5.0432062352798885 + 2.9642198478962527i, -0.6272859719867453 - 1.0342821557347708i, -5.038142339172577 + 4.8798573875784585i, 23.037647711169612 + 0.5375348161530314i, 12.793540056512635 - 26.386015558668642i, -13.075865853449518 - 16.22602339332504i},
	{-1.7016467529919745 + 4.596882243427781i, -2.63422573e+00 - 3.57084951e-01i, -3.2611564914624944 - 6.480184433128504i, 0.3859395342298871 + 38.18128130559367i, 71.7309713862825 + 6.324908210730939i, 23.24322042869128 - 18.31830433914631i},
}
var zri2port_tto = [][]complex128{
	{-135.38071065989854 - 135.32994923857862i, 1877.9289340101516 - 764.9949238578687i},
	{8.629441624365489 + 20.81218274111675i, -247.05583756345177 + 1.2182741116751714i},
}

var ari2port_yto = [][]complex128{
	{0.54 + 1.52i, 0.1 - 0.2i},
	{-17.944 - 9.112i, 0.94 + 2.22i},
}
var hri2port_yto = [][]complex128{
	{-0.06022023399862354 - 0.07054370268410186i, 0.10220233998623536 - 0.29456297315898144i},
	{-0.16173434273916035 + 0.38196834136269786i, -6.382656572608397 + 5.380316586373022i},
}
var sri_yto = [][]complex128{
	{0.9996706300451921 - 0.0006932948472067046i, 0.0011445769021416127 + 0.0011607878524626525i, 0.001068942768334788 - 0.0028442084845683446i, 0.0022160393891765273 + 0.0011674228251476182i, 0.0003314616938696058 + 0.0008981662607645546i, -0.0008682673177250466 + 0.0016921388777036972i},
	{0.0012597326769750727 + 0.0002827914383497637i, 1.001116638933391 - 0.0004705598521810572i, -0.0006788870619611476 + 0.0012260567669577505i, -0.0005815100323264405 + 0.0009199858383121262i, 0.00045882555392623425 - 0.00012197070356456674i, 0.00011187512175320213 + 0.000357007050650171i},
	{-0.00023466160984567108 - 0.0011221485487272476i, 0.0013052409233297174 - 0.00041901742796324193i, 0.9947583133444131 + 0.0037492363715425636i, -0.0013269526292792322 - 0.002217949947646791i, 0.0014624294973328705 + 0.0014755634218526215i, 0.004930489677855865 - 0.003345115274062671i},
	{9.213254628798184e-05 - 0.00041013341371864964i, -0.0014950907784356399 + 0.0006845448123573017i, -0.0021545680670922995 - 0.0009705948103039264i, 1.0004878435891298 + 0.0011269124373110562i, 0.0016995936999227856 - 7.244510992976405e-05i, -0.0006379568262119273 + 0.0004915510081826557i},
	{0.0009794761012259096 + 0.001800173779029024i, -0.0007917550532868384 - 0.0011130691413567262i, 0.0009490700481819825 + 0.002670397572296479i, -0.0018576415234108912 - 0.0013269889461891005i, 1.000179407986476 + 0.000334025094826218i, 0.00034703052192053274 - 0.0030202395292078243i},
	{-0.0008490178962090311 + 0.0010852346445165972i, -0.0008479964602857715 + 0.0015925509385091852i, 0.0013962269273515293 - 0.0031213195644546055i, 0.0013484294042964695 + 0.002189855605134114i, 0.0014484499571329745 - 0.0014825906599242322i, 0.9965140925640767 + 0.005819543131564853i},
}
var sri2port_yto = [][]complex128{
	{1.0023259105903242 + 0.002211638860683582i, 0.0012900160608105282 - 0.0007624183250778405i},
	{0.0017798813484278392 - 0.0008978990102085915i, 1.0036663924368414 + 0.00309839620049962i},
}
var tri_yto = [][]complex128{
	{-13.916372483107537 - 57.55999387218317i, -120.57627846211307 + 279.28902732558515i, 125.39509703110821 + 227.38401713348628i, 14.24528734055563 + 57.74136496101629i, 120.95267856016052 - 279.23409623126577i, -126.31817213993888 - 226.94628164660512i},
	{365.62226737737745 + 34.40678957991363i, 139.65514046525297 - 18.232326795688483i, 47.59103910506014 + 161.96389267754407i, -366.0543523951661 - 33.80023428349163i, -138.6871762980821 + 18.468939367895267i, -48.96515410424441 - 162.51198907732845i},
	{190.81460515903845 + 7.407170795615723i, -0.6150131325525281 + 92.03876186625408i, -70.28962894687005 - 78.79585747725831i, -190.53160521962454 - 7.623040723891016i, 0.6901352764132788 - 92.03046427516041i, 71.15344973938774 + 78.77585170515883i},
	{-14.48569628 - 57.51643854i, -121.25379333890069 + 279.0535332624634i, 125.71404720562164 + 227.18615643292208i, 14.813883344023868 + 57.69778185391592i, 121.62901990294532 - 279.00041553352133i, -126.63679781853014 - 226.74590601769967i},
	{365.33899278256956 + 34.38587378141976i, 139.85175665130575 - 18.399303042602565i, 47.22287583565095 + 161.51750838848633i, -365.77073224978784 - 33.78193385068306i, -138.88584081478982 + 18.635970461931713i, -48.59346510038469 - 162.06736064013194i},
	{191.41466362019116 + 6.80585363115825i, -0.7900963945930141 + 92.5373345021611i, -71.31909195 - 78.93552587i, -191.12946923426736 - 7.022267609435145i, 0.8628325733195865 - 92.53021852770316i, 72.18374529583315 + 78.9147256305905i},
}
var tri2port_yto = [][]complex128{
	{-449.3390000000015 - 229.67200000000705i, 448.401 + 227.448i},
	{-448.801 - 228.148i, 447.85900000000146 + 225.93200000000695i},
}
var yri_yto = [][]complex128{
	{-2.9 + 4i, -8.8 + 4.3i, 2.7 - 5.9i, -2.7 + 6.3i, -9.2 + 8.5i, 3.6 - 0.8i},
	{-5.5 + 5.4i, -8.2 - 7.9i, -0.1 + 4.2i, 7.1 - 0.3i, 5.6 + 2.6i, 2.4 + 1.4i},
	{-5.2 - 1.2i, 7.3 + 5.6i, -4.9 + 2.2i, 8.2 - 2.7i, -5.4 + 4.6i, -6.1 + 3.1i},
	{-10.1 + 1.1i, -3.1 + 7.8i, 6.8 - 3i, -9 + 9.1i, -8.1 - 7.6i, 5.8 - 3i},
	{-8.2 + 6.4i, 4.5 + 4.3i, -9.7 + 4.3i, 0.03 - 3.1i, -8.2 + 7.9i, -9.5 + 3i},
	{-2.1 - 4.8i, 7.1 + 2.9i, -6.4 + 4.8i, 1.5 - 8.4i, -0.9 + 0.9i, -4.3 + 9.4i},
}
var yri2port_yto = [][]complex128{
	{-7 + 8.2i, -1.7 - 2.9i},
	{-2 - 4i, -5 + 5.2i},
}
var zri_yto = [][]complex128{
	{0.008178699943536104 + 0.017332257579755574i, -0.028650788338285275 - 0.02912849179075317i, -0.026736715041983806 + 0.07137903588311009i, -0.05557538728700337 - 0.02916980189672834i, -0.00813829672458126 - 0.022415189206463113i, 0.021591701371984687 - 0.04265803281673013i},
	{-0.031451073938142914 - 0.007068700449943486i, -0.02790540508071871 + 0.011768012613798293i, 0.017035909320615623 - 0.03078592597753147i, 0.014574237220864567 - 0.02295703757054806i, -0.011499491931066388 + 0.0030960741223153712i, -0.002837666896937609 - 0.008828862862907752i},
	{0.005913246682057675 + 0.028279141463048812i, -0.03260282401644721 + 0.010660108090979517i, 0.13109220189748128 - 0.09431566046566794i, 0.03356008 + 0.05550155i, -0.03671446550625466 - 0.0371177281664255i, -0.12327950346806653 + 0.08455869953505385i},
	{-0.0023116942570377812 + 0.010315289390809144i, 0.03728761747128428 - 0.017161913319429525i, 0.05407144881230645 + 0.024242426912816032i, -0.012253317910264856 - 0.02815063393462842i, -0.042500666784912364 + 0.001810124311804771i, 0.015762869264219465 - 0.01239529685292725i},
	{-0.024425767910427023 - 0.04500892099573167i, 0.019899467736838698 + 0.02793603316820571i, -0.0239011215827713 - 0.06690132525333634i, 0.04661588264334084 + 0.03307976409374254i, -0.004631463981149941 - 0.008356778040465886i, -0.00830731899539942 + 0.07577675000363765i},
	{0.021192753031232955 - 0.027215448240840403i, 0.021016773391608724 - 0.04000732818558646i, -0.03461478684418069 + 0.0785143657764065i, -0.03417184 - 0.05468922i, -0.03605997524215898 + 0.03725654650071418i, 0.08673053705591469 - 0.14634870059106778i},
}
var zri2port_yto = [][]complex128{
	{-0.05812120656661984 - 0.055193912492474365i, -0.03210344843932956 + 0.01908864367917805i},
	{-0.04430445754593996 + 0.022497894402508073i, -0.09159151566675147 - 0.07720787501362912i},
}

var ari2port_zto = [][]complex128{
	{-0.94 - 2.22i, 17.944 + 9.112i},
	{-0.1 + 0.2i, -0.54 - 1.52i},
}
var hri2port_zto = [][]complex128{
	{-9.04688701 + 8.59123751i, -0.1264411990776326 + 0.44850115295926213i},
	{0.20753266717909302 - 0.5841660261337432i, -0.09607993850883935 - 0.09992313604919292i},
}
var sri_zto = [][]complex128{
	{-1.028717961784003 + 0.440234536908726i, -0.3350839838141198 + 0.24196084684466915i, 0.06634110811149471 - 0.18173699188254513i, 0.022450243623397783 + 0.33973082472039234i, -0.25392437122785755 + 0.5815611786577279i, 0.1273681064490779 - 0.031671864009029334i},
	{-0.11252915519910302 + 0.25170607407421774i, -1.2612357154001426 - 0.5806373116161577i, 0.03251386025172889 + 0.21073215621722427i, 0.3672019131706386 + 0.06859952274958084i, 0.44560712440752015 + 0.31411834857451315i, 0.20929854123853542 + 0.11228430422132704i},
	{-0.1886218606393445 - 0.038417683901384905i, 0.3939168753650158 + 0.2279248851184109i, -1.1870498782445182 + 0.2559802952909334i, 0.1825609090528612 - 0.32010699704526224i, -0.19583829696902583 + 0.2687955010901915i, -0.23069491010117094 + 0.3270053580701752i},
	{-0.4366391779776497 + 0.3201313344397508i, -0.32614784018950943 + 0.6547063374853146i, 0.17326461890860254 - 0.32868097340989i, -1.1808633787789804 + 0.5314729791813929i, -0.5177795707338307 - 0.17544754696268733i, 0.1175289040312403 - 0.3501502129958069i},
	{-0.359467770550046 + 0.3893888135167872i, 0.34599261670739334 + 0.31116648566561617i, -0.38878058 + 0.3797565i, -0.1442041271335246 - 0.29728530382690954i, -1.2962733147602565 + 0.5530060829327172i, -0.3389830629800219 + 0.37684096670930795i},
	{-0.22957628256003712 - 0.2391815891849416i, 0.26896434099951605 + 0.06272279472765897i, -0.12983341538494125 + 0.3138506821337098i, -0.22384744 - 0.45526898i, -0.11229757242831377 - 0.029014723083481655i, -1.034778518798715 + 0.5169468127766439i},
}
var sri2port_zto = [][]complex128{
	{-1.2423848322180107 + 0.4107122125512722i, -0.1244248897260839 - 0.1146154627893207i},
	{-0.15502282151344582 - 0.16315547181999063i, -1.1901851305493953 + 0.236934048763864i},
}
var tri_zto = [][]complex128{
	{0.8978803453519895 + 0.06320162201239071i, 1.2113882951596415 + 1.3893169641089544i, -7.73253934e-01 - 1.16194842i, 0.34041269294230186 + 0.15287395760679262i, 6.31342539e-01 + 0.93141417i, -0.5552788299458957 - 1.5427190562408084i},
	{-0.1439229082336008 + 1.9606689010737968i, 1.7205644759854752 + 0.16016872767922452i, -2.685559775832344 + 1.1104795070434033i, -7.19927547e-01 + 1.91378582i, 1.8067888445198252 + 0.007912365137622834i, -3.0939604192959145 + 0.5353609937595109i},
	{-1.4525720719861277 + 0.5353686773680802i, 1.6148421215071065 + 0.5296957445134758i, 0.39821759625873576 + 0.5244339776550276i, -1.0063889027352695 - 0.07144812864061523i, 1.4415434575089945 + 0.9550318329041253i, -5.74672500e-01 + 0.44824721i},
	{-0.5808216326246263 + 0.07205246362408202i, -0.9324573181665216 - 1.1119771380076617i, -0.048846104371379595 + 1.782593842557639i, -3.69484436e-04 + 0.07458636i, -0.15864685432750725 - 1.1007051720449912i, -0.4684468537627223 + 1.8000595350464172i},
	{0.004567934690905884 - 1.3281463211953954i, -6.42208491e-01 + 0.06172396i, 1.3516821507023695 - 1.5298416826017882i, -0.004541964176935975 - 1.3358131718432218i, -0.95126482552232 + 0.24347493856402747i, 1.68313378e+00 - 1.21394747i},
	{1.889447303301564 - 0.7457389640922445i, -1.5415045756888712 - 0.44816160590334964i, 0.40903028986200896 - 0.6079213755687769i, 1.1377079967575168 - 0.1504654911875647i, -1.3670226807399821 - 1.021152405962514i, 1.495898752704417 - 0.40556582432283034i},
}
var tri2port_zto = [][]complex128{
	{1.58056 - 6.96112i, 2.47944 - 5.25888i},
	{-2.87944 + 4.55888i, -3.06056 + 3.22112i},
}
var yri_zto = [][]complex128{
	{0.008178699943536104 + 0.017332257579755574i, -0.028650788338285275 - 0.02912849179075317i, -0.026736715041983806 + 0.07137903588311009i, -0.05557538728700337 - 0.02916980189672834i, -0.00813829672458126 - 0.022415189206463113i, 0.021591701371984687 - 0.04265803281673013i},
	{-0.031451073938142914 - 0.007068700449943486i, -0.02790540508071871 + 0.011768012613798293i, 0.017035909320615623 - 0.03078592597753147i, 0.014574237220864567 - 0.02295703757054806i, -0.011499491931066388 + 0.0030960741223153712i, -0.002837666896937609 - 0.008828862862907752i},
	{0.005913246682057675 + 0.028279141463048812i, -0.03260282401644721 + 0.010660108090979517i, 0.13109220189748128 - 0.09431566046566794i, 0.03356008 + 0.05550155i, -0.03671446550625466 - 0.0371177281664255i, -0.12327950346806653 + 0.08455869953505385i},
	{-0.0023116942570377812 + 0.010315289390809144i, 0.03728761747128428 - 0.017161913319429525i, 0.05407144881230645 + 0.024242426912816032i, -0.012253317910264856 - 0.02815063393462842i, -0.042500666784912364 + 0.001810124311804771i, 0.015762869264219465 - 0.01239529685292725i},
	{-0.024425767910427023 - 0.04500892099573167i, 0.019899467736838698 + 0.02793603316820571i, -0.0239011215827713 - 0.06690132525333634i, 0.04661588264334084 + 0.03307976409374254i, -0.004631463981149941 - 0.008356778040465886i, -0.00830731899539942 + 0.07577675000363765i},
	{0.021192753031232955 - 0.027215448240840403i, 0.021016773391608724 - 0.04000732818558646i, -0.03461478684418069 + 0.0785143657764065i, -0.03417184 - 0.05468922i, -0.03605997524215898 + 0.03725654650071418i, 0.08673053705591469 - 0.14634870059106778i},
}
var yri2port_zto = [][]complex128{
	{-0.05812120656661984 - 0.055193912492474365i, -0.03210344843932956 + 0.01908864367917805i},
	{-0.04430445754593996 + 0.022497894402508073i, -0.09159151566675147 - 0.07720787501362912i},
}
var zri_zto = [][]complex128{
	{-2.9 + 4i, -8.8 + 4.3i, 2.7 - 5.9i, -2.7 + 6.3i, -9.2 + 8.5i, 3.6 - 0.8i},
	{-5.5 + 5.4i, -8.2 - 7.9i, -0.1 + 4.2i, 7.1 - 0.3i, 5.6 + 2.6i, 2.4 + 1.4i},
	{-5.2 - 1.2i, 7.3 + 5.6i, -4.9 + 2.2i, 8.2 - 2.7i, -5.4 + 4.6i, -6.1 + 3.1i},
	{-10.1 + 1.1i, -3.1 + 7.8i, 6.8 - 3i, -9 + 9.1i, -8.1 - 7.6i, 5.8 - 3i},
	{-8.2 + 6.4i, 4.5 + 4.3i, -9.7 + 4.3i, 0.03 - 3.1i, -8.2 + 7.9i, -9.5 + 3i},
	{-2.1 - 4.8i, 7.1 + 2.9i, -6.4 + 4.8i, 1.5 - 8.4i, -0.9 + 0.9i, -4.3 + 9.4i},
}
var zri2port_zto = [][]complex128{
	{-7 + 8.2i, -1.7 - 2.9i},
	{-2 - 4i, -5 + 5.2i},
}

var smatdb = [][]float64{
	{-0.011638199740124467, 143.333, -29.4815424208956, -126.062, -49.97273402183755, -107.176, -37.922079170687375, 51.6093, -60.02138490615103, 55.6223, -56.33693866310331, -118.165},
	{-29.50081527775497, -126.086, -0.017397877553168975, 142.397, -28.94208711722665, -124.776, -59.87410284037229, 55.5605, -38.856495024129934, 50.2964, -70.42737163025224, 47.4311},
	{-50.13923614464823, -107.006, -28.94177099566982, -124.774, -0.017162891498142555, 145.341, -56.4782906315203, -117.922, -70.75105170391471, 47.1432, -36.81617834038599, 52.3353},
	{-37.890613838568655, 51.6691, -59.88901159485402, 55.4923, -56.47441174864956, -117.908, -0.011107666268990988, 142.902, -29.37977680005382, -126.722, -50.62593624242556, -108.518},
	{-59.584655840941814, 56.9223, -38.85718035248474, 50.2897, -70.75314831243483, 47.1081, -29.37987909771291, -126.723, -0.016832181153131866, 141.968, -28.786256506244086, -125.418},
	{-56.65997798941029, -118.955, -70.44993468392997, 47.1776, -36.81617834038599, 52.3349, -50.62761902514525, -108.525, -28.786280391380686, -125.419, -0.016884397738919827, 144.864},
}
var smatma = [][]float64{
	{0.998661, 143.333, 0.0335678, -126.062, 0.00317222, -107.176, 0.0127027, 51.6093, 0.000997541, 55.6223, 0.00152459, -118.165},
	{0.0334934, -126.086, 0.997999, 142.397, 0.0357187, -124.776, 0.0010146, 55.5605, 0.0114071, 50.2964, 0.000301045, 47.4311},
	{0.00311199, -107.006, 0.03572, -124.774, 0.998026, 145.341, 0.00149998, -117.922, 0.000290033, 47.1432, 0.0144275, 52.3353},
	{0.0127488, 51.6691, 0.00101286, 55.4923, 0.00150065, -117.908, 0.998722, 142.902, 0.0339634, -126.722, 0.00294241, -108.518},
	{0.00104898, 56.9223, 0.0114062, 50.2897, 0.000289963, 47.1081, 0.033963, -126.723, 0.998064, 141.968, 0.0363653, -125.418},
	{0.00146893, -118.955, 0.000300264, 47.1776, 0.0144275, 52.3349, 0.00294184, -108.525, 0.0363652, -125.419, 0.998058, 144.864},
}
var smatri = [][]float64{
	{-0.8010456803103849, 0.5963636566534487, -0.01976003285054651, -0.027135554141850926, -0.000936781556903727, -0.0030307457899079937, 0.007888637948893376, 0.009956303460686567, 0.0005632573554477263, 0.000823303831046008, -0.0007196252615640324, -0.0013440662747866631},
	{-0.019727576112061034, -0.027067149541516906, -0.7906723874950947, 0.6089656637724387, -0.02037285904023424, -0.029338918593171528, 0.000573792527853044, 0.0008367647787640287, 0.007287039807950417, 0.00877615982348464, 0.00020364981327439363, 0.00022170892083609226},
	{-0.0009101694635168231, -0.002975915541103527, -0.020372576349513963, -0.02934069755277172, -0.8209274823200242, 0.5675684676298588, -0.0007023943167160737, -0.0013253611674728363, 0.00019727126739970267, 0.00021261041401473023, 0.008815771670156503, 0.01142081109683835},
	{0.007906833511942851, 0.01000069423862251, 0.0005738023948638905, 0.0008346473454390566, -0.0007023840445734197, -0.0013261248344061296, -0.7965856451003317, 0.6024092822193964, -0.020307836383985284, -0.0272232312733313, -0.0009345169644723741, -0.002790063556841195},
	{0.0005725079747132832, 0.000878973070855811, 0.007287491000943224, 0.008774615282231554, 0.00019735383454970294, 0.0002124382389249367, -0.020308072337130084, -0.027222556216305108, -0.7861418576018637, 0.6149087150319883, -0.021075045001074582, -0.029635747371927063},
	{-0.0007111421316101127, -0.0012853140524980754, 0.0002040978802522298, 0.00022023288803570294, 0.008815851402244421, 0.011420749550863333, -0.000934676728571004, -0.0027894088977179024, -0.02107550428408785, -0.029635298044922622, -0.8162001246780365, 0.5744015388554922},
}

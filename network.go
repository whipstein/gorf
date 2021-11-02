package gorf

import (
	"bytes"
	"encoding/gob"
	"errors"
	"fmt"
	"io"
	"io/fs"
	"regexp"
	"strconv"
	"strings"

	"github.com/whipstein/golinalg/mat"
)

type NoiseNetwork struct {
	Freq  *Frequency
	NFmin *mat.Vector
	Gopt  *mat.CVector
	Rn    *mat.Vector
}

func NewNoiseNetwork() *NoiseNetwork {
	return &NoiseNetwork{NewFrequency(), vf(0), cvf(0), vf(0)}
}

func (n *NoiseNetwork) ReadTouchstone(r *Reader) *NoiseNetwork {
	line := make([]byte, 500)
	format := "%f %f %f %f %f"
	freq := 0.
	m, a := 0., 0.

	for {
		if r.Readlnraw(&line); r.eof {
			break
		}
		line = bytes.TrimSpace(line)

		if len(line) == 0 || line[0] == '!' {
			continue
		} else {
			n.NFmin.Data = append(n.NFmin.Data, 0.)
			n.Gopt.Data = append(n.Gopt.Data, complex(0., 0.))
			n.Rn.Data = append(n.Rn.Data, 0.)

			if _, ok := fmt.Sscanf(string(line), format, &freq, n.NFmin.GetPtr(len(n.NFmin.Data)-1), &m, &a, n.Rn.GetPtr(len(n.Rn.Data)-1)); errors.Is(ok, fs.ErrExist) {
				if ok != io.EOF {
					panic(ok)
				}
				break
			}

			n.Freq.Append(freq)
			n.Gopt.Set(len(n.Gopt.Data)-1, complex(convertNumber(m, a, MA, RI)))
		}
	}

	return n
}

type Network struct {
	Name      string
	Comments  string
	NPorts    int
	PortNames []string
	Z0        *mat.CVector
	Freq      *Frequency
	Param     RFParam
	Data      []*mat.CMatrix
	Noise     *NoiseNetwork
}

func NewNetwork() *Network {
	return &Network{Freq: NewFrequency()}
}

func (n *Network) DeepCopy() *Network {
	b := bytes.Buffer{}
	e := gob.NewEncoder(&b)
	if err := e.Encode(n); err != nil {
		panic(err)
	}
	d := gob.NewDecoder(&b)
	result := &Network{}
	if len(n.Data) == 0 {
		result.Data = make([]*mat.CMatrix, 0)
	}
	if err := d.Decode(result); err != nil {
		panic(err)
	}
	return result
}

func (n *Network) Setup(p byte, z string) *Network {
	switch p {
	case 'a':
		n.Param = A
	case 'h':
		n.Param = H
	case 's':
		n.Param = S
	case 'y':
		n.Param = Y
	case 'z':
		n.Param = Z
	default:
		panic("parameter type not recogized")
	}

	i, ok := strconv.ParseFloat(z, 64)
	if errors.Is(ok, fs.ErrExist) {
		fmt.Printf("Number does not exist: %v", z)
		panic(ok)
	}
	n.Z0.SetReAll(float64(i))

	return n
}

func (n *Network) SetPorts(x int) *Network {
	n.NPorts = x
	n.PortNames = make([]string, x)
	n.Z0 = cvf(x)

	return n
}

func (n *Network) ReadTouchstone(f string) *Network {
	enc := RI

	p, ok := regexp.Compile(`s(\d+)p`)
	if errors.Is(ok, fs.ErrExist) {
		panic(ok)
	}

	x := strings.Split(f, ".")

	m := p.FindStringSubmatch(x[len(x)-1])
	if len(m) == 0 {
		panic("File is not a touchstone!")
	}

	i, ok := strconv.Atoi(m[1])
	if errors.Is(ok, fs.ErrExist) {
		fmt.Printf("Number does not exist: %v", m[1])
		panic(ok)
	}
	n.SetPorts(i)

	r := NewReader(f)
	line := make([]byte, 500)

	_i := 0
	freq := 0.
	sr := make([]float64, n.NPorts*n.NPorts)
	si := make([]float64, n.NPorts*n.NPorts)
	noisePattern, ok := regexp.Compile(`\!.*noise parameters`)
	if ok != nil {
		panic(ok)
	}

	// read lines until data starts
	for {
		if r.Readlnraw(&line); r.eof {
			break
		}
		line = bytes.TrimSpace(line)

		if len(line) == 0 {
			continue
		} else if noisePattern.Match(bytes.ToLower(line)) {
			n.Noise = NewNoiseNetwork()
			n.Noise.ReadTouchstone(r)
			break
		} else if line[0] == '!' {
			n.Comments += string(line) + "\n"
		} else if line[0] == '#' {
			res := strings.Split(strings.ToLower(string(line)), " ")
			n.Setup(res[2][0], res[5])
			n.Freq.Setup(res[1])
			switch res[3] {
			case "ri":
				enc = RI
			case "ma":
				enc = MA
			case "db":
				enc = DB
			}
		} else {
			if n.NPorts == 1 {
				format := "%f %f %f"

				if _, ok := fmt.Sscanf(string(line), format, &freq, &sr[0], &si[0]); errors.Is(ok, fs.ErrExist) {
					if ok != io.EOF {
						panic(ok)
					}
					break
				}

				if i != 0 {
					n.Data = append(n.Data, cmf(n.NPorts, n.NPorts, opts))
				}
				n.Freq.Append(freq)
				n.Data[_i].Set(0, 0, complex(convertNumber(sr[0], si[0], enc, RI)))
				_i++
			} else if n.NPorts == 2 {
				format := "%f %f %f %f %f %f %f %f %f"
				if _, ok := fmt.Sscanf(string(line), format, &freq, &sr[0], &si[0], &sr[2], &si[2], &sr[1], &si[1], &sr[3], &si[3]); errors.Is(ok, fs.ErrExist) {
					if ok != io.EOF {
						panic(ok)
					}
					break
				}

				n.Data = append(n.Data, cmf(n.NPorts, n.NPorts, opts))
				n.Freq.Append(freq)
				for i := 0; i < n.NPorts; i++ {
					for j := 0; j < n.NPorts; j++ {
						n.Data[_i].Set(i, j, complex(convertNumber(sr[i*2+j], si[i*2+j], enc, RI)))
					}
				}
				_i++
			} else {
				for i := 0; i < n.NPorts; i++ {
					addrs := make([]interface{}, 0)
					format := ""
					if i == 0 {
						format = "%f"
						n.Data = append(n.Data, cmf(n.NPorts, n.NPorts, opts))
						addrs = append(addrs, &freq)
					} else {
						if r.Readlnraw(&line); r.eof {
							break
						}
						line = bytes.TrimSpace(line)
					}
					for j := 0; j < n.NPorts; j++ {
						format += " %f %f"
						addrs = append(addrs, &sr[j], &si[j])
					}
					format = strings.TrimSpace(format)

					if _, ok := fmt.Sscanf(string(line), format, addrs...); errors.Is(ok, fs.ErrExist) {
						if ok != io.EOF {
							panic(ok)
						}
						break
					}

					if i == 0 {
						n.Freq.Append(freq)
					}
					for j := 0; j < n.NPorts; j++ {
						n.Data[_i].Set(i, j, complex(convertNumber(sr[j], si[j], enc, RI)))
					}
				}
				_i++
			}
		}
	}

	return n
}

func (n *Network) AtoH() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		AtoH(n.Data[i])
	}
	return n
}
func (n *Network) AtoS() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		AtoS(n.Data[i], n.Z0)
	}
	return n
}
func (n *Network) AtoT() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		AtoT(n.Data[i], n.Z0)
	}
	return n
}
func (n *Network) AtoY() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		AtoY(n.Data[i])
	}
	return n
}
func (n *Network) AtoZ() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		AtoZ(n.Data[i])
	}
	return n
}
func (n *Network) HtoA() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		HtoA(n.Data[i])
	}
	return n
}
func (n *Network) HtoS() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		HtoS(n.Data[i], n.Z0)
	}
	return n
}
func (n *Network) HtoT() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		HtoT(n.Data[i], n.Z0)
	}
	return n
}
func (n *Network) HtoY() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		HtoY(n.Data[i])
	}
	return n
}
func (n *Network) HtoZ() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		HtoZ(n.Data[i])
	}
	return n
}
func (n *Network) StoA() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		StoA(n.Data[i], n.Z0)
	}
	return n
}
func (n *Network) StoH() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		StoH(n.Data[i], n.Z0)
	}
	return n
}
func (n *Network) StoT() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		StoT(n.Data[i], n.Z0)
	}
	return n
}
func (n *Network) StoY() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		StoY(n.Data[i], n.Z0)
	}
	return n
}
func (n *Network) StoZ() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		StoZ(n.Data[i], n.Z0)
	}
	return n
}
func (n *Network) TtoA() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		TtoA(n.Data[i], n.Z0)
	}
	return n
}
func (n *Network) TtoH() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		TtoH(n.Data[i], n.Z0)
	}
	return n
}
func (n *Network) TtoS() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		TtoS(n.Data[i], n.Z0)
	}
	return n
}
func (n *Network) TtoY() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		TtoY(n.Data[i], n.Z0)
	}
	return n
}
func (n *Network) TtoZ() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		TtoZ(n.Data[i], n.Z0)
	}
	return n
}
func (n *Network) YtoA() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		YtoA(n.Data[i])
	}
	return n
}
func (n *Network) YtoH() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		YtoH(n.Data[i])
	}
	return n
}
func (n *Network) YtoS() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		YtoS(n.Data[i], n.Z0)
	}
	return n
}
func (n *Network) YtoT() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		YtoT(n.Data[i], n.Z0)
	}
	return n
}
func (n *Network) YtoZ() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		YtoZ(n.Data[i])
	}
	return n
}
func (n *Network) ZtoA() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		ZtoA(n.Data[i])
	}
	return n
}
func (n *Network) ZtoH() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		ZtoH(n.Data[i])
	}
	return n
}
func (n *Network) ZtoS() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		ZtoS(n.Data[i], n.Z0)
	}
	return n
}
func (n *Network) ZtoT() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		ZtoT(n.Data[i], n.Z0)
	}
	return n
}
func (n *Network) ZtoY() *Network {
	for i := 0; i < n.Freq.NPts; i++ {
		ZtoY(n.Data[i])
	}
	return n
}

func (n *Network) A() []*mat.CMatrix {
	net := n.DeepCopy()

	switch n.Param {
	case A:
	case H:
		net.HtoA()
	case S:
		net.StoA()
	case T:
		net.TtoA()
	case Y:
		net.YtoA()
	case Z:
		net.ZtoA()
	}

	return net.Data
}
func (n *Network) H() []*mat.CMatrix {
	net := n.DeepCopy()

	switch n.Param {
	case A:
		net.AtoH()
	case H:
	case S:
		net.StoH()
	case T:
		net.TtoH()
	case Y:
		net.YtoH()
	case Z:
		net.ZtoH()
	}

	return net.Data
}
func (n *Network) S() []*mat.CMatrix {
	net := n.DeepCopy()

	switch n.Param {
	case A:
		net.AtoS()
	case H:
		net.HtoS()
	case S:
	case T:
		net.TtoS()
	case Y:
		net.YtoS()
	case Z:
		net.ZtoS()
	}

	return net.Data
}
func (n *Network) T() []*mat.CMatrix {
	net := n.DeepCopy()

	switch n.Param {
	case A:
		net.AtoT()
	case H:
		net.HtoT()
	case S:
		net.StoT()
	case T:
	case Y:
		net.YtoT()
	case Z:
		net.ZtoT()
	}

	return net.Data
}
func (n *Network) Y() []*mat.CMatrix {
	net := n.DeepCopy()

	switch n.Param {
	case A:
		net.AtoY()
	case H:
		net.HtoY()
	case S:
		net.StoY()
	case T:
		net.TtoY()
	case Y:
	case Z:
		net.ZtoY()
	}

	return net.Data
}
func (n *Network) Z() []*mat.CMatrix {
	net := n.DeepCopy()

	switch n.Param {
	case A:
		net.AtoZ()
	case H:
		net.HtoZ()
	case S:
		net.StoZ()
	case T:
		net.TtoZ()
	case Y:
		net.YtoZ()
	case Z:
	}

	return net.Data
}

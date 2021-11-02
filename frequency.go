package gorf

import (
	"math"
	"strings"

	"github.com/whipstein/golinalg/mat"
)

type FreqUnit int

const (
	Hz  FreqUnit = iota
	KHz          = 1e3
	MHz          = 1e6
	GHz          = 1e9
	THz          = 1e12
)

type Sweep int

const (
	Lin Sweep = iota
	Log
)

type Frequency struct {
	Start      float64
	Stop       float64
	NPts       int
	Unit       FreqUnit
	SweepType  Sweep
	Freq       *mat.Vector
	FreqScaled *mat.Vector
	W          *mat.Vector
}

func NewFrequency() *Frequency {
	return &Frequency{Freq: vf(0), FreqScaled: vf(0), W: vf(0)}
}

func (f *Frequency) Setup(unit string) *Frequency {
	switch strings.ToLower(unit) {
	case "hz":
		f.Unit = Hz
	case "khz":
		f.Unit = KHz
	case "mhz":
		f.Unit = MHz
	case "ghz":
		f.Unit = GHz
	default:
		panic("frequency unit not recognized")
	}

	return f
}

func (f *Frequency) Append(x float64) *Frequency {
	if len(f.Freq.Data) == 0 {
		f.Start = x
	}
	f.NPts++
	f.Stop = x
	f.FreqScaled.Append(x)
	f.Freq.Append(x * float64(f.Unit))
	f.W.Append(x * float64(f.Unit) * 2 * math.Pi)
	return f
}

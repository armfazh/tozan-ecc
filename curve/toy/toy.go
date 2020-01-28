// Package toy provides small instances of elliptic curves.
package toy

import (
	"fmt"
	"math/big"

	C "github.com/armfazh/tozan-ecc/curve"
	GF "github.com/armfazh/tozan-ecc/field"
)

// ID is an identifier of a toy curve.
type ID string

const (
	W0    ID = "W0"
	W1    ID = "W1"
	W1ISO ID = "W1ISO"
	W2    ID = "W2"
	W3    ID = "W3"
	WC0   ID = "WC0"
	M0    ID = "M0"
	M1    ID = "M1"
	E0    ID = "E0"
	E1    ID = "E1"
)

type params struct {
	model C.Model
	p     int
	a, b  int
	h, r  int
	x, y  int
}

// Curves is a list of toy curves.
var Curves []ID
var toyCurves map[ID]*params

func init() {
	Curves = make([]ID, 0, 10)
	toyCurves = make(map[ID]*params)

	W0.register(&params{C.Weierstrass, 53, 3, 2, 51, 3, 46, 3})
	W1.register(&params{C.Weierstrass, 53, 0, 1, 54, 2, 13, 5})
	W1ISO.register(&params{C.Weierstrass, 53, 38, 22, 54, 2, 41, 45})
	W2.register(&params{C.Weierstrass, 53, 0, 2, 51, 3, 37, 27})
	W3.register(&params{C.Weierstrass, 59, 16, 0, 60, 4, 33, 11})
	WC0.register(&params{C.WeierstrassC, 53, 2, 3, 66, 6, 45, 4})
	M0.register(&params{C.Montgomery, 53, 4, 3, 44, 4, 16, 4})
	M1.register(&params{C.Montgomery, 53, 3, 1, 48, 4, 14, 22})
	E0.register(&params{C.TwistedEdwards, 53, 1, 3, 44, 4, 17, 49})
	E1.register(&params{C.TwistedEdwards, 53, -1, 12, 48, 4, 3, 19})
}

func (id ID) register(p *params) { toyCurves[id] = p; Curves = append(Curves, id) }

// New returns an elliptic curve and a generator point.
func (id ID) New() (C.EllCurve, C.Point, error) {
	if v, ok := toyCurves[id]; ok {
		F := GF.NewFp(fmt.Sprintf("%v", v.p), v.p)
		E := v.model.New(string(id), F,
			F.Elt(v.a), F.Elt(v.b),
			big.NewInt(int64(v.r)), big.NewInt(int64(v.h)))
		P := E.NewPoint(F.Elt(v.x), F.Elt(v.y))
		return E, P, nil
	}
	return nil, nil, fmt.Errorf("curve not supported")
}

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
	W4    ID = "W4"
	WC0   ID = "WC0"
	M0    ID = "M0"
	M1    ID = "M1"
	E0    ID = "E0"
	E1    ID = "E1"
)

type params struct {
	model C.Model
	p, m  int
	a, b  int
	h, r  int
	x, y  interface{}
}

// Curves is a list of toy curves.
var Curves []ID
var toyCurves map[ID]*params

func init() {
	Curves = make([]ID, 0, 10)
	toyCurves = make(map[ID]*params)

	W0.register(&params{model: C.Weierstrass, p: 53, m: 1, a: 3, b: 2, r: 51, h: 3, x: 46, y: 3})
	W1.register(&params{model: C.Weierstrass, p: 53, m: 1, a: 0, b: 1, r: 54, h: 2, x: 13, y: 5})
	W1ISO.register(&params{model: C.Weierstrass, p: 53, m: 1, a: 38, b: 22, r: 54, h: 2, x: 41, y: 45})
	W2.register(&params{model: C.Weierstrass, p: 53, m: 1, a: 0, b: 2, r: 51, h: 3, x: 37, y: 27})
	W3.register(&params{model: C.Weierstrass, p: 59, m: 1, a: 16, b: 0, r: 60, h: 4, x: 33, y: 11})
	WC0.register(&params{model: C.WeierstrassC, p: 53, m: 1, a: 2, b: 3, r: 66, h: 6, x: 45, y: 4})
	M0.register(&params{model: C.Montgomery, p: 53, m: 1, a: 4, b: 3, r: 44, h: 4, x: 16, y: 4})
	M1.register(&params{model: C.Montgomery, p: 53, m: 1, a: 3, b: 1, r: 48, h: 4, x: 14, y: 22})
	E0.register(&params{model: C.TwistedEdwards, p: 53, m: 1, a: 1, b: 3, r: 44, h: 4, x: 17, y: 49})
	E1.register(&params{model: C.TwistedEdwards, p: 53, m: 1, a: -1, b: 12, r: 48, h: 4, x: 3, y: 19})
	W4.register(&params{model: C.Weierstrass, p: 19, m: 2, a: 1, b: 4, r: 399, h: 3, x: []interface{}{0, 1}, y: 17})
}

func (id ID) register(p *params) { toyCurves[id] = p; Curves = append(Curves, id) }

// New returns an elliptic curve and a generator point.
func (id ID) New() (C.EllCurve, C.Point, error) {
	if v, ok := toyCurves[id]; ok {
		var F GF.Field
		if v.m == 1 {
			F = GF.NewFp(fmt.Sprintf("%v", v.p), v.p)
		} else if v.m == 2 {
			F = GF.NewFp2(fmt.Sprintf("%v", v.p), v.p)
		}
		E := v.model.New(string(id), F,
			F.Elt(v.a), F.Elt(v.b),
			big.NewInt(int64(v.r)), big.NewInt(int64(v.h)))
		P := E.NewPoint(F.Elt(v.x), F.Elt(v.y))
		return E, P, nil
	}
	return nil, nil, fmt.Errorf("curve not supported")
}

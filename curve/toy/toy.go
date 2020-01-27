// Package toy provides small instances of elliptic curves.
package toy

import (
	"math/big"

	C "github.com/armfazh/tozan-ecc/curve"
	GF "github.com/armfazh/tozan-ecc/field"
)

// EC is an elliptic curve
type EC struct {
	E C.EllCurve
	P C.Point
}

func init() {
	initCurves()
}

// ToyCurves is
var ToyCurves map[string]EC

func initCurves() {
	ToyCurves = make(map[string]EC)

	var P53, P59 GF.ID
	var f53 = GF.NewFp(P53, 53) // 1mod4, 2mod3
	var f59 = GF.NewFp(P59, 59) // 3mod4, 2mod3

	registerToyCurve("W0", C.NewWeierstrass(C.Custom, f53,
		f53.Elt(3), f53.Elt(2), big.NewInt(51), big.NewInt(3)),
		f53.Elt(46), f53.Elt(3))

	registerToyCurve("W1", C.NewWeierstrass(C.Custom, f53,
		f53.Zero(), f53.One(), big.NewInt(54), big.NewInt(2)),
		f53.Elt(13), f53.Elt(5))

	registerToyCurve("W1iso", C.NewWeierstrass(C.Custom, f53,
		f53.Elt(38), f53.Elt(22), big.NewInt(54), big.NewInt(2)),
		f53.Elt(41), f53.Elt(45))

	registerToyCurve("W2", C.NewWeierstrass(C.Custom, f53,
		f53.Zero(), f53.Elt(2), big.NewInt(51), big.NewInt(3)),
		f53.Elt(37), f53.Elt(27))

	registerToyCurve("W3", C.NewWeierstrass(C.Custom, f59,
		f59.Elt(16), f59.Zero(), big.NewInt(60), big.NewInt(4)),
		f59.Elt(33), f59.Elt(11))

	registerToyCurve("WC0", C.NewWeierstrassC(C.Custom, f53,
		f53.Elt(2), f53.Elt(3), big.NewInt(66), big.NewInt(6)),
		f53.Elt(45), f53.Elt(4))

	registerToyCurve("M0", C.NewMontgomery(C.Custom, f53,
		f53.Elt(4), f53.Elt(3), big.NewInt(44), big.NewInt(4)),
		f53.Elt(16), f53.Elt(4))

	registerToyCurve("M1", C.NewMontgomery(C.Custom, f53,
		f53.Elt(3), f53.Elt(1), big.NewInt(48), big.NewInt(4)),
		f53.Elt(14), f53.Elt(22))

	registerToyCurve("E0", C.NewEdwards(C.Custom, f53,
		f53.Elt(1), f53.Elt(3), big.NewInt(44), big.NewInt(4)),
		f53.Elt(17), f53.Elt(49))

	registerToyCurve("E1", C.NewEdwards(C.Custom, f53,
		f53.Elt(-1), f53.Elt(12), big.NewInt(48), big.NewInt(4)),
		f53.Elt(3), f53.Elt(19))

}

// registerToyCurve is
func registerToyCurve(name string, e C.EllCurve, x, y GF.Elt) {
	n := EC{E: e, P: e.NewPoint(x, y)}
	ToyCurves[name] = n
}

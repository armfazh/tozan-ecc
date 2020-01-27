package curve

import (
	"errors"
	"fmt"
	"math/big"

	GF "github.com/armfazh/tozan-ecc/field"
)

// WCCurve is a Weierstrass curve
type WCCurve struct {
	*params
	RationalMap
}

type WC = *WCCurve

func (e *WCCurve) String() string { return "y^2=x^3+Ax^2+Bx\n" + e.params.String() }

// NewWeierstrassC returns a Weierstrass curve
func NewWeierstrassC(id CurveID, f GF.Field, a, b GF.Elt, r, h *big.Int) *WCCurve {
	if e := (&WCCurve{params: &params{Id: id, F: f, A: a, B: b, R: r, H: h}}); e.IsValid() {
		e.RationalMap = e.ToWeierstrass()
		return e
	}
	panic(errors.New("can't instantiate a WeierstrassC curve"))
}

func (e *WCCurve) NewPoint(x, y GF.Elt) (P Point) {
	if P = (&ptWc{e, &afPoint{x: x, y: y}}); e.IsOnCurve(P) {
		return P
	}
	panic(fmt.Errorf("%v not on %v", P, e))
}

func (e *WCCurve) IsValid() bool {
	F := e.F
	t0 := F.Sqr(e.A)      // A^2
	t1 := F.Add(e.B, e.B) // 2B
	t1 = F.Add(t1, t1)    // 4B
	t0 = F.Sub(t0, t1)    // A^2-4B
	t0 = F.Mul(t0, e.B)   // B(A^2-4B)
	return !F.IsZero(t0)  // B(A^2-4B) != 0
}
func (e *WCCurve) IsEqual(ec EllCurve) bool {
	e0 := ec.(*WECurve)
	return e.F.IsEqual(e0.F) && e.F.AreEqual(e.A, e0.A) && e.F.AreEqual(e.B, e0.B)
}
func (e *WCCurve) Identity() Point             { return &infPoint{} }
func (e *WCCurve) IsOnCurve(p Point) bool      { return e.Codomain().IsOnCurve(e.Push(p)) }
func (e *WCCurve) Add(p, q Point) Point        { return e.Pull(e.Codomain().Add(e.Push(p), e.Push(q))) }
func (e *WCCurve) Double(p Point) Point        { return e.Pull(e.Codomain().Double(e.Push(p))) }
func (e *WCCurve) Neg(p Point) Point           { return e.Pull(e.Codomain().Neg(e.Push(p))) }
func (e *WCCurve) ClearCofactor(p Point) Point { return e.Pull(e.Codomain().ClearCofactor(e.Push(p))) }

// ptWc is an affine point on a WCCurve curve.
type ptWc struct {
	*WCCurve
	*afPoint
}

func (p *ptWc) String() string { return p.afPoint.String() }
func (p *ptWc) Copy() Point    { return &ptWc{p.WCCurve, p.copy()} }
func (p *ptWc) IsEqual(q Point) bool {
	qq := q.(*ptWc)
	return p.WCCurve.IsEqual(qq.WCCurve) && p.isEqual(p.F, qq.afPoint)
}
func (p *ptWc) IsIdentity() bool   { return false }
func (p *ptWc) IsTwoTorsion() bool { return p.F.IsZero(p.y) }

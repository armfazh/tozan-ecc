package field

import (
	"crypto/rand"
	"fmt"
	"io"
	"math/big"
)

// fpElt is a prime field element.
type fpElt struct{ n *big.Int }

func (e fpElt) String() string         { return "0x" + e.n.Text(16) }
func (e fpElt) Copy() Elt              { return &fpElt{new(big.Int).Set(e.n)} }
func (e fpElt) Polynomial() []*big.Int { return []*big.Int{new(big.Int).Set(e.n)} }

// fp implements a prime field.
type fp struct {
	p    *big.Int
	name string
	cte  struct {
		pMinus1div2 *big.Int
		pMinus2     *big.Int
	}
	hasSqrt
}

// NewFp creates a prime field as Z/pZ given p as an int, uint, *big.Int or string.
func NewFp(name string, p interface{}) Field {
	prime := FromType(p)
	if !prime.ProbablyPrime(4) {
		panic(fmt.Errorf("Modulus is not prime p:%v", prime))
	}
	f := fp{p: prime, name: name}
	f.precmp()
	return f
}

func (f *fp) precmp() {
	pMinus1div2 := big.NewInt(1)
	pMinus1div2.Sub(f.p, pMinus1div2)
	pMinus1div2.Rsh(pMinus1div2, 1)

	pMinus2 := big.NewInt(2)
	pMinus2.Sub(f.p, pMinus2)
	f.cte.pMinus1div2 = pMinus1div2
	f.cte.pMinus2 = pMinus2

	t := big.NewInt(16)
	pMod16 := t.Mod(f.p, t).Uint64()
	switch {
	case pMod16%4 == uint64(3):
		f.hasSqrt = generateSqrt3mod4(f)
	case pMod16%8 == uint64(5):
		f.hasSqrt = generateSqrt5mod8(f)
	case pMod16%16 == uint64(9):
		f.hasSqrt = generateSqrt9mod16(f)
	case pMod16%16 == uint64(1):
		f.hasSqrt = generateSqrt1mod16(f)
	}
}

func (f fp) String() string       { return fmt.Sprintf("GF(%v)", f.name) }
func (f fp) Zero() Elt            { return &fpElt{big.NewInt(0)} }
func (f fp) One() Elt             { return &fpElt{big.NewInt(1)} }
func (f fp) Rand(r io.Reader) Elt { e, _ := rand.Int(r, f.p); return &fpElt{e} }
func (f fp) P() *big.Int          { return new(big.Int).Set(f.p) }
func (f fp) Order() *big.Int      { return new(big.Int).Set(f.p) }
func (f fp) Ext() uint            { return uint(1) }
func (f fp) BitLen() int          { return f.p.BitLen() }
func (f fp) Elt(in interface{}) Elt {
	var n *big.Int
	if v, ok := in.([]interface{}); ok && len(v) == 1 {
		n = FromType(v[0])
	} else {
		n = FromType(in)
	}
	return f.mod(n)
}
func (f fp) mod(x *big.Int) Elt { return &fpElt{x.Mod(x, f.p)} }

// Implementing hasPredicates

func (f fp) IsZero(x Elt) bool      { return x.(*fpElt).n.Sign() == 0 }
func (f fp) AreEqual(x, y Elt) bool { return f.IsZero(f.Sub(x, y)) }
func (f fp) IsSquare(x Elt) bool    { return f.AreEqual(f.Exp(x, f.cte.pMinus1div2), f.One()) }
func (f fp) IsEqual(ff Field) bool  { return f.p.Cmp(ff.(fp).p) == 0 }

// Implementing hasArith

func (f fp) Neg(x Elt) Elt    { return f.mod(new(big.Int).Neg(x.(*fpElt).n)) }
func (f fp) Add(x, y Elt) Elt { return f.mod(new(big.Int).Add(x.(*fpElt).n, y.(*fpElt).n)) }
func (f fp) Sub(x, y Elt) Elt { return f.mod(new(big.Int).Sub(x.(*fpElt).n, y.(*fpElt).n)) }
func (f fp) Mul(x, y Elt) Elt { return f.mod(new(big.Int).Mul(x.(*fpElt).n, y.(*fpElt).n)) }
func (f fp) Sqr(x Elt) Elt    { return f.mod(new(big.Int).Mul(x.(*fpElt).n, x.(*fpElt).n)) }
func (f fp) Inv(x Elt) Elt    { return f.Exp(x, f.cte.pMinus2) }
func (f fp) Exp(x Elt, y *big.Int) Elt {
	return &fpElt{new(big.Int).Exp(x.(*fpElt).n, y, f.p)}
}

// Implementing extended operations
func (f fp) Generator() Elt { return f.One() }
func (f fp) Inv0(x Elt) Elt { return f.Inv(x) }
func (f fp) Sgn0(x Elt) int { return int(x.(*fpElt).n.Bit(0)) }
func (f fp) CMov(x, y Elt, b bool) Elt {
	var z big.Int
	if b {
		z.Set(y.(*fpElt).n)
	} else {
		z.Set(x.(*fpElt).n)
	}
	return &fpElt{&z}
}

type sqrt3mod4 struct {
	*fp
	exp *big.Int
}

func generateSqrt3mod4(f *fp) hasSqrt {
	e := big.NewInt(1)
	e.Add(f.p, e)
	e.Rsh(e, 2)
	return sqrt3mod4{exp: e, fp: f}
}

func (s sqrt3mod4) Sqrt(x Elt) Elt { return s.Exp(x, s.exp) }

type sqrt5mod8 struct {
	*fp
	sqrtOne *fpElt
	exp     *big.Int
}

func generateSqrt5mod8(f *fp) hasSqrt {
	// calculates s = sqrt(-1) for p=8*k+5
	// t = 2^k
	// s = 2*t^3+t
	k := big.NewInt(5)
	k.Sub(f.p, k)           // p-5
	k.Rsh(k, 3)             // k = (p-5)/8
	t := f.Exp(f.Elt(2), k) // t = 2^k
	s := f.Sqr(t)           // t^2
	s = f.Add(s, s)         // 2t^2
	s = f.Add(s, f.One())   // 2t^2+1
	s = f.Mul(s, t)         // t(2t^2+1)
	k.Add(k, big.NewInt(1)) // e = k+1 = (p+3)/8
	return sqrt5mod8{fp: f, exp: k, sqrtOne: s.(*fpElt)}
}

func (s sqrt5mod8) Sqrt(x Elt) Elt {
	t0 := s.Exp(x, s.exp)
	t1 := s.Sqr(t0)
	e := s.AreEqual(x, t1)
	t1 = s.Mul(t0, s.sqrtOne)
	return s.CMov(t1, t0, e)
}

type sqrt9mod16 struct {
	*fp
	c1 *fpElt   // c1 = sqrt(-1) in F, i.e., (c1^2) == -1 in F
	c2 *fpElt   // c2 = sqrt(c1) in F, i.e., (c2^2) == c1 in F
	c3 *fpElt   // c3 = sqrt(-c1) in F, i.e., (c3^2) == -c1 in F
	c4 *big.Int // c4 = (q + 7) / 16         # Integer arithmetic
}

func (s sqrt9mod16) Sqrt(x Elt) Elt {
	tv1 := s.Exp(x, s.c4)
	tv2 := s.Mul(s.c1, tv1)
	tv3 := s.Mul(s.c2, tv1)
	tv4 := s.Mul(s.c3, tv1)
	e1 := s.AreEqual(s.Sqr(tv2), x)
	e2 := s.AreEqual(s.Sqr(tv3), x)
	tv1 = s.CMov(tv1, tv2, e1)
	tv2 = s.CMov(tv4, tv3, e2)
	e3 := s.AreEqual(s.Sqr(tv2), x)
	return s.CMov(tv1, tv2, e3)
}

func generateSqrt9mod16(f *fp) hasSqrt {
	panic("not implemented")
}

type sqrt1mod16 struct {
	*fp
	c1 *big.Int
	c3 *big.Int
	c5 *fpElt
}

func (s sqrt1mod16) Sqrt(x Elt) Elt {
	z := s.Exp(x, s.c3)
	t := s.Mul(s.Sqr(z), x)
	z = s.Mul(z, x)
	b := t
	c := s.c5

	one := big.NewInt(1)
	endOuter := big.NewInt(2)
	endInner := new(big.Int)
	for i := new(big.Int).Set(s.c1); i.Cmp(endOuter) != 0; i.Sub(i, one) {
		endInner.Sub(i, endOuter)
		for j := new(big.Int).Set(one); j.Cmp(endInner) != 0; j.Add(j, one) {
			b = s.Sqr(b)
		}
		z = s.CMov(z, s.Mul(z, c), !s.AreEqual(b, s.One()))
		c = s.Sqr(c).(*fpElt)
		t = s.CMov(t, s.Mul(t, c), !s.AreEqual(b, s.One()))
		b = t
	}

	return z
}

func generateSqrt1mod16(f *fp) hasSqrt {
	one := big.NewInt(1)
	two := big.NewInt(2)
	c1 := findC1(f)

	c2 := f.Order()
	c2.Sub(c2, one)
	c2.Div(c2, (&big.Int{}).Exp(two, c1, nil))

	c3 := (&big.Int{}).Sub(c2, one)
	c3.Div(c3, two)

	c4 := findNonSquare(f)
	c5 := f.Exp(c4, c2).(*fpElt)

	fmt.Printf("c1: %v, c2: %v, c3: %v, c4: %v, c5: %v\n", c1, c2, c3, c4, c5)

	return sqrt1mod16{fp: f, c1: c1, c3: c3, c5: c5}
}

// Find the largest integer c1 such that 2^c1 divides q-1
func findC1(f *fp) *big.Int {
	zero := big.NewInt(0)
	one := big.NewInt(1)
	two := big.NewInt(2)

	c1 := big.NewInt(1)
	qMinus1 := (&big.Int{}).Sub(f.Order(), one)
	for ; ; c1.Add(c1, one) {
		qMinus1.Div(qMinus1, two)
		if (&big.Int{}).Mod(qMinus1, two).Cmp(zero) == 1 {
			return c1
		}
	}
	panic("no valid c1 found")
}

// Find a non square in the field
func findNonSquare(f *fp) Elt {
	two := &fpElt{big.NewInt(2)}
	for i := two.Copy(); !f.AreEqual(i, f.Zero()); i = f.Add(i, f.One()) {
		if !f.IsSquare(i) {
			return i
		}
	}
	panic("no non-squares found")
}

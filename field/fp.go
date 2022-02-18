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
	// Calculates c1 = sqrt(-1) mod p, for p=8*k+5
	// x(a) = -1 iff a^(4k+2) = -1
	//   sqrt(-1) = sqrt(a^(4k+2))
	//            = a^(2k+1)
	// Since x(2) = -1, 2 \in QNR, then
	//   sqrt(-1) = 2^(2k+1).
	k := big.NewInt(5)
	k.Sub(f.p, k)            // p-5
	k.Rsh(k, 3)              // k = (p-5)/8
	c1 := f.Exp(f.Elt(2), k) // c1 = 2^(k)
	c1 = f.Sqr(c1)           //    = 2^(2k)
	c1 = f.Add(c1, c1)       //    = 2^(2k+1)
	k.Add(k, big.NewInt(1))  // e = k+1 = (p+3)/8
	return sqrt5mod8{fp: f, exp: k, sqrtOne: c1.(*fpElt)}
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
	// Calculates c1 = sqrt(-1) mod p, for p=16*k+9
	// x(a) = -1 iff a^(8k+4) = -1
	//   c1 = sqrt(-1)
	//      = sqrt(a^(8k+4))
	//      = a^(4k+2)
	//   c2 = sqrt(sqrt(-1))
	//      = sqrt(a^(4k+2))
	//      = a^(2k+1)
	//   c3 = sqrt(-sqrt(-1))
	//      = sqrt(-1) * sqrt(-sqrt(-1))
	//      = c1*c2
	//
	// find a such that x(a) = -1.
	k := big.NewInt(9)
	k.Sub(f.p, k)                 // p-9
	k.Rsh(k, 4)                   // k = (p-9)/16
	a := findNonSquare(f)         // a is QNR
	c2 := f.Exp(a, k)             // c2 = a^(k)
	c2 = f.Sqr(c2)                //    = a^(2k)
	c2 = f.Mul(c2, a)             //    = a^(2k+1)
	c1 := f.Sqr(c2)               // c1 = a^(4k+2)
	c3 := f.Mul(c1, c2)           // c3 = c1*c2
	c4 := k.Add(k, big.NewInt(1)) // e = k+1 = (p+7)/16
	return sqrt9mod16{f, c1.(*fpElt), c2.(*fpElt), c3.(*fpElt), c4}
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
	b := t.Copy()
	c := s.c5.Copy()

	one := big.NewInt(1)
	endInner := new(big.Int)
	for i := new(big.Int).Set(s.c1); i.Cmp(one) > 0; i.Sub(i, one) {
		endInner.Sub(i, one)
		for j := new(big.Int).Set(one); j.Cmp(endInner) < 0; j.Add(j, one) {
			b = s.Sqr(b)
		}
		z = s.CMov(z, s.Mul(z, c), !s.AreEqual(b, s.One()))
		c = s.Sqr(c).(*fpElt)
		t = s.CMov(t, s.Mul(t, c), !s.AreEqual(b, s.One()))
		b = t
	}

	return z
}

// Tonelli-Shanks algorithm.
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

	return sqrt1mod16{fp: f, c1: c1, c3: c3, c5: c5}
}

// Find the largest integer c1 such that 2^c1 divides q-1
func findC1(f *fp) *big.Int {
	zero := big.NewInt(0)
	one := big.NewInt(1)
	two := big.NewInt(2)

	c1 := big.NewInt(1)
	qMinus1 := (&big.Int{}).Sub(f.Order(), one)
	// loops until a valid c1 is found.
	for ; ; c1.Add(c1, one) {
		qMinus1.Div(qMinus1, two)
		if (&big.Int{}).Mod(qMinus1, two).Cmp(zero) == 1 {
			return c1
		}
	}
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

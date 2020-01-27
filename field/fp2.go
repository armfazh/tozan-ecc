package field

import (
	"crypto/rand"
	"io"
	"math/big"
)

type fp2Elt struct {
	a, b *big.Int
}

func (e fp2Elt) String() string {
	return "\na: 0x" + e.a.Text(16) +
		"\nb: 0x" + e.b.Text(16) + " * i"
}

func (e fp2Elt) Copy() Elt { r := &fp2Elt{}; r.a.Set(e.a); r.b.Set(e.b); return r }

type fp2 struct {
	p    *big.Int
	name string
	cte  struct {
		pMinus1div2 *big.Int
	}
}

// NewFp2 creates a quadratic extension field Z/pZ[x] with irreducible polynomial x^2=-1 and given p as an int, uint, *big.Int or string.
func NewFp2(name string, p interface{}) Field {
	prime := FromType(p)
	if !prime.ProbablyPrime(4) {
		panic("p is not prime")
	}
	return fp2{p: prime, name: name}
}

func (f fp2) Elt(in interface{}) Elt {
	var a, b *big.Int
	if v, ok := in.([]interface{}); ok && len(v) == 2 {
		a = FromType(v[0])
		b = FromType(v[1])
	} else {
		a = FromType(in)
		b = big.NewInt(0)
	}
	return f.mod(a, b)
}
func (f fp2) P() *big.Int     { return new(big.Int).Set(f.p) }
func (f fp2) Order() *big.Int { return new(big.Int).Mul(f.p, f.p) }
func (f fp2) String() string  { return "GF(" + f.name + ") Irred: i^2+1" }
func (f fp2) Ext() uint       { return uint(2) }
func (f fp2) Zero() Elt       { return f.Elt(0) }
func (f fp2) One() Elt        { return f.Elt(1) }
func (f fp2) BitLen() int     { return f.p.BitLen() }

func (f fp2) AreEqual(x, y Elt) bool { return f.IsZero(f.Sub(x, y)) }
func (f fp2) IsEqual(ff Field) bool  { return f.p.Cmp(ff.(fp2).p) == 0 }
func (f fp2) IsZero(x Elt) bool {
	e := x.(*fp2Elt)
	return e.a.Mod(e.a, f.p).Sign() == 0 &&
		e.b.Mod(e.b, f.p).Sign() == 0
}

func (f fp2) Rand(r io.Reader) Elt {
	a, _ := rand.Int(r, f.p)
	b, _ := rand.Int(r, f.p)
	return &fp2Elt{a, b}
}

func (f fp2) mod(a, b *big.Int) Elt { return &fp2Elt{a: a.Mod(a, f.p), b: b.Mod(b, f.p)} }
func (f fp2) Add(x, y Elt) Elt {
	a := new(big.Int).Add(x.(fp2Elt).a, y.(fp2Elt).a)
	b := new(big.Int).Add(x.(fp2Elt).b, y.(fp2Elt).b)
	a.Mod(a, f.p)
	b.Mod(b, f.p)
	return fp2Elt{a, b}
}
func (f fp2) Sub(x, y Elt) Elt {
	a := new(big.Int).Sub(x.(fp2Elt).a, y.(fp2Elt).a)
	b := new(big.Int).Sub(x.(fp2Elt).b, y.(fp2Elt).b)
	a.Mod(a, f.p)
	b.Mod(b, f.p)
	return fp2Elt{a, b}
}

func (f fp2) Mul(x, y Elt) Elt          { return nil }
func (f fp2) Sqr(x Elt) Elt             { return nil }
func (f fp2) Inv(x Elt) Elt             { return nil }
func (f fp2) Neg(x Elt) Elt             { return nil }
func (f fp2) Sqrt(x Elt) Elt            { return nil }
func (f fp2) Exp(x Elt, e *big.Int) Elt { return nil }
func (f fp2) IsSquare(x Elt) bool       { return false }

func (f fp2) Generator() Elt { return f.Elt([]string{"0", "1"}) }
func (f fp2) Inv0(x Elt) Elt { return f.Inv(x) }
func (f fp2) CMov(x, y Elt, b bool) Elt {
	var za, zb big.Int
	if b {
		za.Set(y.(*fp2Elt).a)
		zb.Set(y.(*fp2Elt).b)
	} else {
		za.Set(x.(*fp2Elt).a)
		zb.Set(x.(*fp2Elt).b)
	}
	return &fp2Elt{&za, &zb}
}
func (f fp2) GetSgn0(id Sgn0ID) func(Elt) int {
	if id == SignBE {
		return f.Sgn0BE
	}
	if id == SignLE {
		return f.Sgn0LE
	}
	panic("Wrong signID")
}
func (f fp2) Sgn0BE(x Elt) int {
	/* [TODO] */
	sb := x.(*fp2Elt).b.Sign()
	cb := 2*(f.cte.pMinus1div2.Cmp(x.(*fp2Elt).b)&^1) - 1
	sa := x.(*fp2Elt).a.Sign()
	ca := 2*(f.cte.pMinus1div2.Cmp(x.(*fp2Elt).a)&^1) - 1
	return sb*cb + (1-sb)*(sa*ca+(1-sa)*1)
}

func (f fp2) Sgn0LE(x Elt) int {
	/* [TODO] */
	sa := x.(*fp2Elt).a.Sign()
	ca := 1 - 2*int(x.(*fp2Elt).a.Bit(0))
	sb := x.(*fp2Elt).b.Sign()
	cb := 1 - 2*int(x.(*fp2Elt).b.Bit(0))
	return sa*ca + (1-sa)*(sb*cb+(1-sb)*1)
}

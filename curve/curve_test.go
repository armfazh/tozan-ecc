package curve_test

import (
	"testing"

	C "github.com/armfazh/tozan-ecc/curve"
	"github.com/armfazh/tozan-ecc/internal/toy"
)

func TestCurves(t *testing.T) {
	for name, EC := range toy.ToyCurves {
		t.Run(name, func(t *testing.T) {
			testAdd(t, EC.E, EC.P)
		})
	}
}

func testAdd(t *testing.T, e C.EllCurve, g C.Point) {
	order := e.Order().Uint64()
	T := make([]C.Point, order)
	T[0] = e.Identity()
	for i := uint64(1); i < order; i++ {
		T[i] = e.Add(T[i-1], g)
		if !e.IsOnCurve(T[i]) {
			t.Fatalf("point not in the curve: %v\n", T[i])
		}
	}
	for _, P := range T {
		for _, Q := range T {
			R := e.Add(P, Q)
			if !e.IsOnCurve(R) {
				t.Fatalf("point not in the curve: %v\n", R)
			}
		}
	}
}

func BenchmarkCurve(b *testing.B) {
	ec := toy.ToyCurves["W0"]
	e := ec.E
	P := ec.P
	Q := e.Double(P)
	b.Run("double", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			P = e.Double(P)
		}
	})

	b.Run("add", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			P = e.Add(P, Q)
		}
	})
}

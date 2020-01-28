package curve_test

import (
	"testing"

	C "github.com/armfazh/tozan-ecc/curve"
	"github.com/armfazh/tozan-ecc/curve/toy"
)

func TestCurves(t *testing.T) {
	for _, curveID := range toy.Curves {
		e, g, _ := curveID.New()
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
}

func BenchmarkCurve(b *testing.B) {
	E, P, _ := toy.W0.New()
	Q := E.Double(P)
	b.Run("double", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			P = E.Double(P)
		}
	})
	b.Run("add", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			P = E.Add(P, Q)
		}
	})
}

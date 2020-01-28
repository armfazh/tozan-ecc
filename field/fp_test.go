package field_test

import (
	"fmt"
	"testing"

	GF "github.com/armfazh/tozan-ecc/field"
)

func TestSqrt(t *testing.T) {
	var primes = []int{
		607, // 3 mod 4
		613, // 5 mod 8
		// 617, // 9 mod 16
		// 641, // 1 mod 16
	}
	for _, p := range primes {
		testSqrt(t, p)
	}
}

func testSqrt(t *testing.T, p int) {
	F := GF.NewFp(fmt.Sprintf("%v", p), p)
	for i := 0; i < p; i++ {
		x := F.Elt(i)
		if F.IsSquare(x) {
			y := F.Sqrt(x)
			got := F.Sqr(y)
			want := x
			if !F.AreEqual(got, want) {
				t.Fatalf("got: %v\nwant: %v\nF:%v", got, want, F)
			}
		}
	}
}

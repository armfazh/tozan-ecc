package field_test

import (
	"fmt"
	"testing"

	GF "github.com/armfazh/tozan-ecc/field"
)

func TestSqrtF2(t *testing.T) {
	var primes3mod4 = []int{59, 67, 71, 79, 83}
	for _, p := range primes3mod4 {
		testSqrtF2(t, p)
	}
}

func testSqrtF2(t *testing.T, p int) {
	F := GF.NewFp2(fmt.Sprintf("%v", p), p)
	for i := 0; i < p; i++ {
		for j := 0; j < p; j++ {
			x := F.Elt([]interface{}{i, j})
			if F.IsSquare(x) {
				y := F.Sqrt(x)
				got := F.Sqr(y)
				want := x
				if !F.AreEqual(got, want) {
					t.Fatalf("got: %v\nwant: %v\n%v\nF:%v", got, want, i, F)
				}
			}
		}
	}
}

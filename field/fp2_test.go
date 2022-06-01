package field_test

import (
	"fmt"
	"testing"

	GF "github.com/armfazh/tozan-ecc/field"
)

func TestGeneratorF2(t *testing.T) {
	// Field taken from https://github.com/armfazh/h2c-go-ref/blob/master/field/fields.go
	F := GF.NewFp2("BN254G2", "0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47")
	F.Generator()
}

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

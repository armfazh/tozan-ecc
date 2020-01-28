package toy_test

import (
	"testing"

	"github.com/armfazh/tozan-ecc/curve/toy"
)

func TestGenerate(t *testing.T) {
	for _, curveId := range toy.Curves {
		if _, _, err := curveId.New(); err != nil {
			t.Fatalf("Curve: %v %v", curveId, err)
		}
	}
}

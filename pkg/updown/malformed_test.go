package updown

import "testing"

func TestNegativeSizesAndDistancesAreErrors(t *testing.T) {
	if _, _, err := checkArgs(0, -2, 0, 0, 0, 0, 0, 0, 0, 0); err == nil {
		t.Errorf("negative --size-up: expected an error")
	}
	if _, _, err := checkArgs(0, 0, 0, 0, 0, 0, 0, -1, 0, 5); err == nil {
		t.Errorf("negative --dist-down: expected an error")
	}
}

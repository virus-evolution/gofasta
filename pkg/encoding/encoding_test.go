package encoding

import(
  "testing"
)

func intersectionStringArrays(A []string, B []string) ([]string) {
  intersection := make([]string, 0)
  for i := 0; i < len(A); i++ {
    for j := 0; j < len(B); j++ {
      test := A[i] == B[j]
      if test {
        intersection = append(intersection, A[i])
      }
    }
  }
  return intersection
}

func TestEncoding(t *testing.T) {

  nucs := []string{"A", "G", "C", "T", "R", "M", "W", "S", "K", "Y", "V", "H", "D", "B", "N", "-", "?"}

  lookupChar := make(map[string][]string)

  lookupChar["A"] = []string{"A"}
  lookupChar["C"] = []string{"C"}
  lookupChar["G"] = []string{"G"}
  lookupChar["T"] = []string{"T"}
  lookupChar["R"] = []string{"A", "G"}
  lookupChar["Y"] = []string{"C", "T"}
  lookupChar["S"] = []string{"G", "C"}
  lookupChar["W"] = []string{"A", "T"}
  lookupChar["K"] = []string{"G", "T"}
  lookupChar["M"] = []string{"A", "C"}
  lookupChar["B"] = []string{"C", "G", "T"}
  lookupChar["D"] = []string{"A", "G", "T"}
  lookupChar["H"] = []string{"A", "C", "T"}
  lookupChar["V"] = []string{"A", "C", "G"}
  lookupChar["N"] = []string{"A", "C", "G", "T"}
  lookupChar["?"] = []string{"A", "C", "G", "T"}
  lookupChar["-"] = []string{"A", "C", "G", "T"}

  lookupByte := MakeByteDict()

  for i := 0; i < len(nucs); i++  {
    for j := 0; j < len(nucs); j++  {
      nuc1 := nucs[i]
      nuc2 := nucs[j]

      nuc1Chars := lookupChar[nuc1]
      nuc2Chars := lookupChar[nuc2]

      rune1 := []rune(nuc1)[0]
      rune2 := []rune(nuc2)[0]

      byte1 := lookupByte[rune1]
      byte2 := lookupByte[rune2]

      byteDifferent := (byte1 & byte2) < 16
      byteSame := (byte1 & 8 == 8) && byte1 == byte2

      nucDifferent := len(intersectionStringArrays(nuc1Chars, nuc2Chars)) == 0
      nucSame := len(intersectionStringArrays([]string{nuc1}, []string{"A", "C", "G", "T"})) == 1 && nuc1 == nuc2

      test := byteDifferent == nucDifferent && byteSame == nucSame

      if ! test {
        t.Errorf("problem in encoding test: %s %s", nuc1, nuc2)
      }
    }
  }
}

package sam


// getCigarOperationMapNoInsertions is a map of SAM CIGAR operation types to function literals
// that are used to build an aligned sequence. This version DISCARDS insertions relative
// to the reference.
func getCigarOperationMapNoInsertions() map[string]func(int, int, int, []byte) (int, int, []byte) {
	cigarOperationMap := map[string]func(int, int, int, []byte) (int, int, []byte){
		"M": func(query_start, ref_start, length int, seq []byte) (int, int, []byte) {
			return query_start + length, ref_start + length, seq[query_start : query_start+length]
		},
		"I": func(query_start, ref_start, length int, seq []byte) (int, int, []byte) {
			return query_start + length, ref_start, []byte{}
		},

		"D": func(query_start, ref_start, length int, seq []byte) (int, int, []byte) {
			gaps := make([]byte, length)
			for i, _ := range gaps {
				gaps[i] = '-'
			}
			return query_start, ref_start + length, gaps
		},

		"N": func(query_start, ref_start, length int, seq []byte) (int, int, []byte) {
			skip := make([]byte, length)
			for i, _ := range skip {
				skip[i] = '*'
			}
			return query_start, ref_start + length, skip
		},

		"S": func(query_start, ref_start, length int, seq []byte) (int, int, []byte) {
			return query_start + length, ref_start, []byte{}
		},
		"H": func(query_start, ref_start, length int, seq []byte) (int, int, []byte) {
			return query_start, ref_start, []byte{}
		},
		"P": func(query_start, ref_start, length int, seq []byte) (int, int, []byte) {
			return query_start, ref_start, []byte{}
		},
		"=": func(query_start, ref_start, length int, seq []byte) (int, int, []byte) {
			return query_start + length, ref_start + length, seq[query_start : query_start+length]
		},
		"X": func(query_start, ref_start, length int, seq []byte) (int, int, []byte) {
			return query_start + length, ref_start + length, seq[query_start : query_start+length]
		}}
	return cigarOperationMap
}

// getCigarOperationMapWithInsertions is a map of SAM CIGAR operation types to function literals
// that are used to build an aligned sequence. This version INCLUDES insertions relative
// to the reference.
func getCigarOperationMapWithInsertions() map[string]func(int, int, int, []byte) (int, int, []byte) {
	cigarOperationMap := map[string]func(int, int, int, []byte) (int, int, []byte){
		"M": func(query_start, ref_start, length int, seq []byte) (int, int, []byte) {
			return query_start + length, ref_start + length, seq[query_start : query_start+length]
		},
		"I": func(query_start, ref_start, length int, seq []byte) (int, int, []byte) {
			return query_start + length, ref_start, seq[query_start : query_start+length]
		},

		"D": func(query_start, ref_start, length int, seq []byte) (int, int, []byte) {
			gaps := make([]byte, length)
			for i, _ := range gaps {
				gaps[i] = '-'
			}
			return query_start, ref_start + length, gaps
		},

		"N": func(query_start, ref_start, length int, seq []byte) (int, int, []byte) {
			skip := make([]byte, length)
			for i, _ := range skip {
				skip[i] = '*'
			}
			return query_start, ref_start + length, skip
		},

		"S": func(query_start, ref_start, length int, seq []byte) (int, int, []byte) {
			return query_start + length, ref_start, []byte{}
		},
		"H": func(query_start, ref_start, length int, seq []byte) (int, int, []byte) {
			return query_start, ref_start, []byte{}
		},
		"P": func(query_start, ref_start, length int, seq []byte) (int, int, []byte) {
			return query_start, ref_start, []byte{}
		},
		"=": func(query_start, ref_start, length int, seq []byte) (int, int, []byte) {
			return query_start + length, ref_start + length, seq[query_start : query_start+length]
		},
		"X": func(query_start, ref_start, length int, seq []byte) (int, int, []byte) {
			return query_start + length, ref_start + length, seq[query_start : query_start+length]
		}}
	return cigarOperationMap
}

// getCigarOperationMapNoInsertionsWithRef is a map of SAM CIGAR operation types to function literals
// that are used to build an aligned sequence. This version DISCARDS insertions relative
// to the reference.
func getCigarOperationMapNoInsertionsWithRef() map[string]func(int, int, int, []byte, []byte) (int, int, []byte, []byte) {
	cigarOperationMap := map[string]func(int, int, int, []byte, []byte) (int, int, []byte, []byte){
		"M": func(query_start, ref_start, length int, seq []byte, refseq []byte) (int, int, []byte, []byte) {
			return query_start + length, ref_start + length, seq[query_start : query_start+length], refseq[ref_start : ref_start+length]
		},
		"I": func(query_start, ref_start, length int, seq []byte, refseq []byte) (int, int, []byte, []byte) {
			return query_start + length, ref_start, []byte{}, []byte{}
		},

		"D": func(query_start, ref_start, length int, seq []byte, refseq []byte) (int, int, []byte, []byte) {
			gaps := make([]byte, length)
			for i, _ := range gaps {
				gaps[i] = '-'
			}
			return query_start, ref_start + length, gaps, refseq[ref_start : ref_start+length]
		},

		"N": func(query_start, ref_start, length int, seq []byte, refseq []byte) (int, int, []byte, []byte) {
			skip := make([]byte, length)
			for i, _ := range skip {
				skip[i] = '*'
			}
			return query_start, ref_start + length, skip, refseq[ref_start : ref_start+length]
		},

		"S": func(query_start, ref_start, length int, seq []byte, refseq []byte) (int, int, []byte, []byte) {
			return query_start + length, ref_start, []byte{}, []byte{}
		},
		"H": func(query_start, ref_start, length int, seq []byte, refseq []byte) (int, int, []byte, []byte) {
			return query_start, ref_start, []byte{}, []byte{}
		},
		"P": func(query_start, ref_start, length int, seq []byte, refseq []byte) (int, int, []byte, []byte) {
			return query_start, ref_start, []byte{}, []byte{}
		},
		"=": func(query_start, ref_start, length int, seq []byte, refseq []byte) (int, int, []byte, []byte) {
			return query_start + length, ref_start + length, seq[query_start : query_start+length], refseq[ref_start : ref_start+length]
		},
		"X": func(query_start, ref_start, length int, seq []byte, refseq []byte) (int, int, []byte, []byte) {
			return query_start + length, ref_start + length, seq[query_start : query_start+length], refseq[ref_start : ref_start+length]
		}}
	return cigarOperationMap
}

// getCigarOperationMapWithInsertionsWithRef is a map of SAM CIGAR operation types to function literals
// that are used to build an aligned sequence. This version INCLUDES insertions relative
// to the reference.
func getCigarOperationMapWithInsertionsWithRef() map[string]func(int, int, int, []byte, []byte) (int, int, []byte, []byte) {
	cigarOperationMap := map[string]func(int, int, int, []byte, []byte) (int, int, []byte, []byte){
		"M": func(query_start, ref_start, length int, seq []byte, refseq []byte) (int, int, []byte, []byte) {
			return query_start + length, ref_start + length, seq[query_start : query_start+length], refseq[ref_start : ref_start+length]
		},
		"I": func(query_start, ref_start, length int, seq []byte, refseq []byte) (int, int, []byte, []byte) {
			gaps := make([]byte, length)
			for i, _ := range gaps {
				gaps[i] = '-'
			}
			return query_start + length, ref_start, seq[query_start : query_start+length], gaps
		},

		"D": func(query_start, ref_start, length int, seq []byte, refseq []byte) (int, int, []byte, []byte) {
			gaps := make([]byte, length)
			for i, _ := range gaps {
				gaps[i] = '-'
			}
			return query_start, ref_start + length, gaps, refseq[ref_start : ref_start+length]
		},

		"N": func(query_start, ref_start, length int, seq []byte, refseq []byte) (int, int, []byte, []byte) {
			skip := make([]byte, length)
			for i, _ := range skip {
				skip[i] = '*'
			}
			return query_start, ref_start + length, skip, refseq[ref_start : ref_start+length]
		},

		"S": func(query_start, ref_start, length int, seq []byte, refseq []byte) (int, int, []byte, []byte) {
			return query_start + length, ref_start, []byte{}, []byte{}
		},
		"H": func(query_start, ref_start, length int, seq []byte, refseq []byte) (int, int, []byte, []byte) {
			return query_start, ref_start, []byte{}, []byte{}
		},
		"P": func(query_start, ref_start, length int, seq []byte, refseq []byte) (int, int, []byte, []byte) {
			return query_start, ref_start, []byte{}, []byte{}
		},
		"=": func(query_start, ref_start, length int, seq []byte, refseq []byte) (int, int, []byte, []byte) {
			return query_start + length, ref_start + length, seq[query_start : query_start+length], refseq[ref_start : ref_start+length]
		},
		"X": func(query_start, ref_start, length int, seq []byte, refseq []byte) (int, int, []byte, []byte) {
			return query_start + length, ref_start + length, seq[query_start : query_start+length], refseq[ref_start : ref_start+length]
		}}
	return cigarOperationMap
}

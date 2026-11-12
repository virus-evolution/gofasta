#!/usr/bin/env python3
"""Write small synthetic fixtures for gofasta's end-to-end command tests into the directory given.

query.csv and target.csv are not written here: they are the output of `gofasta updown list -r ref.fasta`
for query.fasta and target.fasta,
and target_bad.csv is target.csv with two malformed lines added by hand.

Reference: 120 bp; geneA on the forward strand (11..61), geneB on the reverse strand (complement(70..108)).
"""
import os
import random
import sys

from Bio.Seq import Seq

random.seed(20261112)
out = sys.argv[1]
os.makedirs(out, exist_ok=True)

STOPS = {"TAA", "TAG", "TGA"}
CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOPS]


def rand_nt(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def orf(ncodons):
    return "ATG" + "".join(random.choice(CODONS) for _ in range(ncodons - 2)) + "TAA"


# layout (1-based): 1..10 intergenic, 11..61 geneA (17 codons), 62..69 intergenic,
# 70..108 geneB reverse (13 codons), 109..120 intergenic
geneA = orf(17)
geneB_rc = str(Seq(orf(13)).reverse_complement())
ref = rand_nt(10) + geneA + rand_nt(8) + geneB_rc + rand_nt(12)
assert len(ref) == 120
protA = str(Seq(geneA).translate())
protB = str(Seq(geneB_rc).reverse_complement().translate())
assert protA.endswith("*") and protB.endswith("*") and "*" not in protA[:-1] + protB[:-1]


def write(name, text):
    with open(os.path.join(out, name), "w") as fh:
        fh.write(text)


def fasta(records, wrap=None):
    lines = []
    for name, seq in records:
        lines.append(">" + name)
        if wrap:
            lines.extend(seq[i:i + wrap] for i in range(0, len(seq), wrap))
        else:
            lines.append(seq)
    return "\n".join(lines) + "\n"


def sub(seq, pos, base):
    """1-based substitution"""
    return seq[:pos - 1] + base + seq[pos:]


def other(base):
    return {"A": "G", "G": "A", "C": "T", "T": "C"}[base]


write("ref.fasta", fasta([("ref", ref)]))
write("ref_short.fasta", fasta([("ref_short", ref[:100])]))  # every alignment record is longer than this

# annotation
def genbank():
    lines = ["LOCUS       ref                      120 bp    RNA     linear   VRL 12-NOV-2026",
             "FEATURES             Location/Qualifiers",
             "     CDS             11..61",
             '                     /gene="geneA"',
             "                     /codon_start=1",
             '                     /translation="%s"' % protA[:-1],
             "     CDS             complement(70..108)",
             '                     /gene="geneB"',
             "                     /codon_start=1",
             '                     /translation="%s"' % protB[:-1],
             "ORIGIN"]
    low = ref.lower()
    for i in range(0, len(low), 60):
        chunk = low[i:i + 60]
        lines.append("%9d %s" % (i + 1, " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))))
    lines.append("//")
    return "\n".join(lines) + "\n"


write("ref.gb", genbank())
write("ref.gff", "\n".join([
    "##gff-version 3",
    "##sequence-region ref 1 120",
    "ref\tsynthetic\tgene\t11\t61\t.\t+\t.\tID=gene-A;Name=geneA",
    "ref\tsynthetic\tCDS\t11\t61\t.\t+\t0\tID=cds-A;Parent=gene-A;Name=geneA",
    "ref\tsynthetic\tCDS\t70\t108\t.\t-\t0\tID=cds-B;Name=geneB",
]) + "\n")

# queries in reference coordinates (aligned, no insertions)
q = {}
q["q1"] = sub(sub(ref, 5, other(ref[4])), 14, other(ref[13]))            # intergenic SNP + SNP in geneA codon 2
q["q2"] = sub(sub(ref, 80, other(ref[79])), 100, other(ref[99]))         # two SNPs in geneB
q["q3"] = ref[:20] + "---" + ref[23:]                                     # in-frame deletion in geneA
q["q4"] = sub(sub(ref, 30, "R"), 31, "N")                                # ambiguity codes
q["q5"] = "N" * 15 + sub(ref, 50, other(ref[49]))[15:110] + "-" * 10      # missing ends
q["q6"] = q["q1"]                                                         # identical to q1
q["q7"] = sub(sub(sub(ref, 14, other(ref[13])), 40, other(ref[39])), 115, other(ref[114]))
q["q8"] = ref
msa = [("ref", ref)] + [(k, q[k]) for k in sorted(q)]
write("msa.fasta", fasta(msa))
write("query.fasta", fasta([(k, q[k]) for k in ("q1", "q2", "q5")]))
write("target.fasta", fasta([(k, q[k]) for k in ("q3", "q4", "q6", "q7", "q8")], wrap=50))

# an msa where the reference has a gap column (an insertion in some queries)
ins_ref = ref[:44] + "---" + ref[44:]
ins_q1 = ref[:44] + "AAA" + ref[44:]
ins_q2 = sub(ref, 90, other(ref[89]))[:44] + "---" + sub(ref, 90, other(ref[89]))[44:]
write("msa_ins.fasta", fasta([("ref", ins_ref), ("i1", ins_q1), ("i2", ins_q2)]))

# malformed inputs
write("msa_badrow.fasta", fasta(msa[:3] + [("short", ref[:100])] + msa[3:]))
write("msa_badnuc.fasta", fasta(msa[:3] + [("badnuc", ref[:50] + "J" + ref[51:])] + msa[3:]))
write("empty.fasta", "")
write("ignore.txt", "q6\n")


# SAM
def cigar_seq(start, ops, subs=()):
    """Build SEQ for a record that starts at 1-based reference position `start` with CIGAR ops
    [(len, op)], applying 1-based reference-coordinate substitutions `subs` {pos: base}."""
    seq = []
    rpos = start
    subs = dict(subs)
    ins = iter("ACGTTGCA" * 4)
    for n, op in ops:
        if op in "M=X":
            for _ in range(n):
                seq.append(subs.get(rpos, ref[rpos - 1]))
                rpos += 1
        elif op == "I":
            seq.extend(next(ins) for _ in range(n))
        elif op == "S":
            seq.extend(rand_nt(n))
        elif op in "DN":
            rpos += n
    return "".join(seq)


def rec(qname, flag, pos, ops, subs=(), rname="ref"):
    cig = "".join("%d%s" % (n, op) for n, op in ops)
    seq = cigar_seq(pos, ops, subs)
    return "\t".join([qname, str(flag), rname, str(pos), "60", cig, "*", "0", "0", seq, "*"])


header = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ref\tLN:120\n"
records = [
    rec("s1", 0, 1, [(120, "M")], {14: other(ref[13]), 85: other(ref[84])}),
    rec("s2", 0, 1, [(30, "M"), (3, "D"), (40, "M"), (2, "I"), (47, "M")]),
    rec("s3", 0, 1, [(70, "M"), (50, "S")], {60: other(ref[59])}),
    rec("s3", 2048, 61, [(60, "S"), (60, "M")], {100: other(ref[99])}),
    rec("s4", 16, 5, [(5, "H"), (110, "M")], {95: other(ref[94])}),
    "s5\t4\t*\t0\t0\t*\t*\t0\t0\t" + ref[:50] + "\t*",
    rec("s6", 256, 1, [(120, "M")], {50: other(ref[49])}),
    rec("s7", 0, 30, [(60, "M")], {40: other(ref[39])}),
    rec("s8", 0, 1, [(44, "M"), (3, "I"), (76, "M")]),
]
write("aln.sam", header + "\n".join(records) + "\n")
bad = records[:4] + ["s9\t0\tref\t1\t60\t50M\t*\t0\t0\t" + ref[:40] + "\t*"] + records[4:]
write("aln_badrecord.sam", header + "\n".join(bad) + "\n")
two_sq = header.replace("@SQ\tSN:ref\tLN:120\n", "@SQ\tSN:ref\tLN:120\n@SQ\tSN:other\tLN:120\n")
write("aln_twosq.sam", two_sq + "\n".join(records[:2] + [rec("s10", 0, 1, [(120, "M")], rname="other")] + records[2:]) + "\n")
write("aln_headeronly.sam", header)
# s3's supplementary record runs past the end of the reference, so the whole of s3 is unusable
past_end = "s3\t2048\tref\t100\t60\t30M\t*\t0\t0\t" + ref[90:120] + "\t*"
write("aln_badquery.sam", header + "\n".join(records[:3] + [past_end] + records[4:]) + "\n")
# a blank line between records, and no newline at the end of the last record
write("aln_untidy.sam", header + "\n".join(records[:2]) + "\n\n" + "\n".join(records[2:]))
write("empty.sam", "")
print("wrote fixtures to", out)

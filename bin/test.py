import regex

r1_20nt = "TCGAGTCCCTCAAGTCCTTC"
r2_20nt = "CCGCCACCGCCGCCGCCGCC"
mismatch = 2

seq = "ATACGDCDTCGAGTCCCTCAAGTCCTTCCAGCAGCAGCAGCAGCAGCCTGCCACCGCCGCCGCCGCC"

r1_search = regex.search("(" + r1_20nt + ")" + "{e<=" + str(mismatch) + "}", str(seq))
r2_search = regex.search("(" + r2_20nt + ")" + "{e<=" + str(mismatch) + "}", str(seq))

import os
import sys
import assembly_statistics as asmstats

"""
Examples
Human:
chr - hs1 1 0 249250621 chr1
chr - hs2 2 0 243199373 chr2
chr - hs3 3 0 198022430 chr3
chr - hs4 4 0 191154276 chr4
chr - hs5 5 0 180915260 chr5

Droso:
chr - dmX x 0 22422827 set1-6-qual-1
chr - dm2l 2l 0 22570260 set1-6-qual-2
chr - dm2r 2r 0 21146708 set1-6-qual-3
chr - dm3l 3l 0 23849507 set1-6-qual-4
chr - dm3r 3r 0 27905053 set1-6-qual-5
chr - dm4 4 0 1351857 set1-6-qual-6

Column assignment as known/guessed
type    ?      id       shortid     ?       seq-len     annotation/highlight?
chr     -      hs1      1           0       249250621   chr1


mummer2Circos
#############
Ideogram (outer ring/track):
length and the contigs in the ideogram are derived from the reference



The query contigs are displayed on the outer ring

A.acidoterrestris reference has 2 contigs and is ~4,34Mb long (ref for BC02)
BC02 reference has 5 contigs and is ~4,37Mb long (ref is A.acidoterrestris)

A.acidocaldarius reference has 4 contigs and is ~3,2Mb long(ref for BC03)
BC03 has 2 contigs and is ~2,8Mb long (ref is A.acidocaldarius)
"""

def main(asm_fa):
    contigs, total_bp = asmstats.parse_assembly(asm_fa)
    for tig in contigs:
        ktf_out = f"{karyotype}.{"_".join(asm_fa.split(".")[:-1])}.txt"
        with open(ktf_out, "w") as out:
        


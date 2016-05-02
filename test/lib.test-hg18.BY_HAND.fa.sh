#!/bin/bash
ref_fa=${REF_FA:-$HOME/git/tegetype/data/ref.hg18.fa}

echo ">test-hg18_1"
samtools faidx $ref_fa chr1:1573331-1574330 | tail -n +2
samtools faidx $ref_fa chr1:1574632-1575631 | tail -n +2

echo ">test-hg18_2"
samtools faidx $ref_fa chr1:1597760-1598769 | tail -n +2
samtools faidx $ref_fa chr1:1599063-1600062 | tail -n +2

echo ">test-hg18_3"
samtools faidx $ref_fa chr1:814630-815688 | tail -n +2
echo "gggcgcggtggctcacgcctgtaatcccagcactttgggaggccgaggcgggcggatcacctgaggtcgggagttcgagaccagcctgaccaacatggagaaaccccgtctctactaaaaatacaaaaattagccgggcgtggtggcgcgcgcctgtaatcccagctactcgggaggctgaggcaggagaatcgcttgaacccgggaggcggaggttgcggtgagccgagatcgcgccattgcactccagcctgggcaacaagagcgaaactccgtctcaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
samtools faidx $ref_fa chr1:815630-816688 | tail -n +2

echo ">test-hg18_4"
samtools faidx $ref_fa chr1:99382117-99383116 | tail -n +2
echo "ggccgggcgcggtggctcacgcctgtaatcccagcactttgggaggccgaggcgggtggatcatgaggtcaggagatcgagaccatcctggctaacaaggtgaaaccccgtctctactaaaaatacaaaaaattagccgggcgcggtggcgggcgcctgtagtcccagctactggggaggctgaggcaggagaatggcgtgaacccgggaagcggagcttgcagtgagccgagattgcgccactgcagtccgcagtccggcctgggcgacagagcgagactccgtctcaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
samtools faidx $ref_fa chr1:99383138-99384137 | tail -n +2

echo ">test-hg18_5"
samtools faidx $ref_fa chr13:26930990-26931989 | tail -n +2
echo "tcacgcctgtaatcccagcactttgggaggccgaggcgggcggatcacgaggtcaggagatcgagaccacggtgaaaccccgtctctactaaaaatacaaaaaattagccgggcgcggtggcgggcgcctgtagtcccagctactcgggaggctgaggcaggagaatggcgtgaacccgggaggcggagcttgcagtgagccgagatcgcgccactgcactccagcctgggcgacagagcgagactccgtctcaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
samtools faidx $ref_fa chr13:26931990-26932989 | tail -n +2

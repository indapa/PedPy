#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
from itertools import  izip


def pairwise(iterable):
    "s -> (s0,s1), (s2,s3), (s4, s5), ..."
    a = iter(iterable)
    return izip(a, a)

""" take a ped file that is space-delimited for all fields and place tabs between the first 6 fields and tabs between each genotype, but white space between each allele ina genotype """
def main():
    usage = "usage: %prog [options] file\n\nconvert file.ped file.map to file.vcf"
    parser = OptionParser(usage)
    (options, args)=parser.parse_args()
    pedfile=args[0]
    pedfh=open(pedfile, 'r')
    for line in pedfh:
        fields= line.strip().split(' ')

        meta=fields[0:5]
        genos=fields[5::]
        gstrings=[]
        total= len(genos) + len(meta)
        for (x,y) in pairwise(genos):
            gstring=" ".join([x,y])
            gstrings.append(gstring)
        print "\t".join(meta) + "\t" + "\t".join(gstrings)
if __name__ == "__main__":
    main()

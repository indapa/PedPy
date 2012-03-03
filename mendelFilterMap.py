#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
from PedFile import *
import numpy as np

""" filter genotypes based on affected/unaffected status for dominant mendelian trait, given a ped file"""
# bioinformatics is so much fun!
def main():
    usage = "usage: %prog [options] file\n\nconvert file.ped file.map to file.vcf"
    parser = OptionParser(usage)
    parser.add_option("--model", type="string", dest="model", default = "dominant", help=" inheritance model (dominant recessive) ")
    parser.add_option("--affected", type="string", dest="affected", help="sample name of affecteds (one per line)")
    parser.add_option("--unaffected", type="string", dest="unaffected", help="sample name of unaffecteds (one per line)")
    (options, args)=parser.parse_args()

  
    pedfile=args[0]

  
    mapobjects=[]


    pedobj=Pedigree(pedfile)

    (affecteds, unaffecteds)=pedobj.getAffecteds()
    print affecteds
    print unaffecteds
    print len(affecteds), len(unaffecteds)

    # initially the rows are the individuals, column are the markers
    # we transpose on the fly and make the rows the genotypes and the columns the individuals, since this is the structure of the VCF
    genotype_matrix=np.array( pedobj.getGenotypeMatrix() ).transpose()
    names=genotype_matrix[0,:]

   
    (nrow, ncol) = np.shape(genotype_matrix)
    print nrow, ncol
    for i in range(nrow):
        if i==0: continue
        gls=zip(names, list(genotype_matrix[i,:]))
        alleles=[]
        for (sample, genotype) in gls:
            (a1, a2)=genotype.split(' ')
            alleles.append(a1)
            alleles.append(a2)
        seg_alleles=set(alleles)
        print seg_alleles
        print gls


if __name__ == "__main__":
    main()

#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
from PedFile import *
import numpy as np


def isSegregating( gstring):
    (a1, a2) = gstring.split(' ')
    if a1 != '0' and  a2 != '0':
        if a1 != '1' or a2 != '1':
            return True
        else:
            return False
    return False
        

""" filter genotypes based on affected/unaffected status for dominant mendelian trait, given a Mendelsoft pre  file"""

def main():
    usage = "usage: %prog [options] file.pre"
    parser = OptionParser(usage)
    parser.add_option("--model", type="string", dest="model", default = "dominant", help=" inheritance model (dominant recessive) ")
    parser.add_option("--affected", type="string", dest="affected", help="sample name of affecteds (one per line)")
    parser.add_option("--unaffected", type="string", dest="unaffected", help="sample name of unaffecteds (one per line)")
    parser.add_option("--v", action="store_true", dest="verbose",  help="Print additional info of filtering results to output")
    (options, args)=parser.parse_args()

    affecteds=[]
    unaffecteds=[]

    pedfile=args[0]
    fileName, fileExtension = os.path.splitext(pedfile)
    (chr, numb, pos, corrected)=fileName.split('.')
    
    if options.affected != None:
        affectedfh=open(options.affected, 'r')
        for line in affectedfh:
            affecteds.append(line.strip() )
            
    if options.unaffected != None:
        unaffectedfh=open (options.unaffected, 'r')
        for line in unaffectedfh:
            unaffecteds.append ( line.strip() )

    #check if any overlapping samples between unaffected and affected
    if len( list( set(unaffecteds).intersection( set(affecteds) ) ) ) != 0:
        sys.stderr.write("check list of affected and unaffecteds for overlapping samples!\n")
        exit(1)

    

    mapobjects=[]

    mendelsoft=True
    pedobj=Pedigree(pedfile,mendelsoft)

    #print affecteds
    #print unaffecteds
    #print len(affecteds), len(unaffecteds)

    # initially the rows are the individuals, column are the markers
    # we transpose on the fly and make the rows the genotypes and the columns the individuals, since this is the structure of the VCF
    genotype_matrix=np.array( pedobj.getGenotypeMatrix() ).transpose()
    #print genotype_matrix
    names=genotype_matrix[0,:]
    #print names
   
    (nrow, ncol) = np.shape(genotype_matrix)
    #print genotype_matrix
   # print nrow, ncol
   # print genotype_matrix
    for i in range(nrow):
        if i==0: continue
        gls=zip(names, list(genotype_matrix[i,:]))
        alleles=[]
        affected_genotypes=[]
        unaffected_genotypes=[]
        for (sample, genotype) in gls:
            if genotype == '0 0': continue # skip over samples not assigned genotype, even after mendelsoft correction
            (a1, a2)=genotype.split(' ')
            alleles.append(a1)
            alleles.append(a2)
        seg_alleles=set(alleles)
        #print seg_alleles
        #print gls
        #iterate thru and see if they are in affected or unaffected list
        for (sample, genotype) in gls: 
            
            if sample in affecteds: # if so ...
                if genotype == '0 0':
                    affecteds= filter((lambda x: x != sample), affecteds)
                    continue
                affected_genotypes.append( ( sample, isSegregating(genotype), genotype ) ) # are they segregating for a non-ref allele?
            if sample in unaffecteds:
                if genotype == '0 0':
                    unaffecteds= filter((lambda x: x != sample), unaffecteds)
                    continue
                unaffected_genotypes.append( (sample, isSegregating(genotype), genotype ) )
     
        shared_affected_segregating = filter( lambda x, segregating=True: segregating in x, affected_genotypes)
        shared_unaffected_segregating = filter( lambda x, segregating=False: segregating in x, unaffected_genotypes)
     
        
        if len(shared_affected_segregating) == len(affecteds): # all the affecteds have at least 1 segregating allele
            if len(shared_unaffected_segregating) == len(unaffecteds): # all the unaffecteds are homoz ref
                print  chr, numb, pos, corrected,"PASS"
                #print genotype_matrix
            elif len(shared_unaffected_segregating) < len(unaffecteds): # some but not all unaffectes are homoz ref
                print  chr, numb, pos, corrected,"MAYBE"
                if options.verbose == True:
                    print "unaffecteds  segregating for alt: ", filter( lambda x, segregating=True: segregating in x, unaffected_genotypes)
        else:
            print  chr, numb, pos, corrected,"NO" # since not all the affectes have a mutant alele, its not a cnadidate
            #print shared_affected_segregating, affected_genotypes
            if options.verbose == True:
                print "affecteds not segregating for alt: ", filter( lambda x, segregating=False: segregating in x, affected_genotypes), len(shared_affected_segregating), len(affecteds)
            
            

if __name__ == "__main__":
    main()

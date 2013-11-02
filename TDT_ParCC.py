"""
:File: TDT_ParCC.py
:Author: Jack A. Kosmicki & Kamil Slowikowski
:Last updated: 2013-10-27

find transmitted, untransmitted variants in families and count them up for each variant
find variants in cases and controls

Steps:
1)  read in the family relations file and save the family relations in a global dictionary
2)  read in the individuals in the VCF and save them in a global dictionary
3)  assign the optional parameters as global variables
4)  make a pool of workers
5)  read in the VCF file one line at a time and pass a set of lines to each worker
6)  for each line in the VCF, apply a series of filters to determine if it should be examined
7)  Upon passing the filters, examine each individual in the VCF
8)  For doTDT, only look at the probands (ignore cases/controls/parents)
9)  Apply filters to the proband
10) If the proband passes the filters, look at the parents and apply filters and determine transmission
11) After determining transmission, update the number of transmitted and untransmitted variants for that specific variant
12) print out the number of transmissions and untransmissions

Usage:
    TDT_ParCC.py family <family_Relations_File> [options]
    TDT_ParCC.py doTDT <vcf_File> <family_Relations_File> <outputFile_Name> <cores> <chunksize> [options]
    TDT_ParCC.py doCaseControl <vcf_File> <family_Relations_File> <outputFile_Name> <cores> <chunksize> [options]

Options:
    --pl=VALUE           Specify PL threshold [default: 30]
    --dp=DP_VAL          Specify minimum depth [default: 10]
    --ac=AC_VAL          Specify the minimum allele count value [default: 5]
    --gq_Par=GQ_PAR_VAL  GQ threshold for parents [default: 30]
    --gq_Kid=GQ_KID_VAL  GQ threshold for child [default: 30]
    --gq_CC=GQ_CC_VAL    GQ threshold for cases and controls [default: 30]
    --syn                Flag to examine all synonymous variants as opposed to
                         removing those above the minimum allele count
    -h, --help           Print this message and exit.
    -v, --version        Print the version and exit.
"""

__version__ = 1.0
__author__ = 'Jack A. Kosmicki <jkosmicki@fas.harvard.edu>'
__date__ = '10/27/2013'


import os
import gzip
import sys
import VCF
import csv
import Family
import filters
import multiprocessing as mp
from guppy import hpy
from functools import partial
from docopt import docopt
from collections import OrderedDict
#@profile

def find_VCFindivs(vcfFile):
    """ find all the individuals in the VCF and store them in an array"""
    
    # open up the vcf file and save the header line which includes the individuals
    with open(args['<vcf_File>'], 'rb') as vcf:
        for line in vcf:
            if line.startswith('#CHROM'):
                vcf_indivIDs = line.split('\t')[9:]
                break

    return vcf_indivIDs
    

def doTDT(v, family, AC_Thresh, GQ_Kid_Thresh, GQ_Parent_Thresh, DP_Thresh, 
          PL_Thresh, synFlag):
    """
    This function does the Transmission Disequilibrium Tests (TDT).
    family is a dictionary with:
       Key: individual id      Value: [father id, mother id]   
    """
    #hp.setrelheap()
    # apply filters to the line, if true look at the line otherwise move on
    if not filters.processLine(v.get('FILTER'), v.get('REF'), v.get('ALT'),
                               v.get('SNPEFF_FUNCTIONAL_CLASS'), v.get('AC'), 
                               AC_Thresh, synFlag):
        return None

    # array of [transmissions, untransmissions]
    TU = [0, 0]

    #loop through all the probands
    for indiv_id in family.keys():
        # indiv_data is their GT:AD:DP:GQ:PL stats
        indiv_data = v.get(indiv_id)

        if indiv_data == None:
            continue

        # check if the individual is a proband
        #if indiv_id not in family:
        #    continue

        #apply quality filters
        if not filters.passFilters(indiv_data, GQ_Kid_Thresh, DP_Thresh):
            continue

        # appply PL filter to child
        if not filters.PhredScaleFilter(indiv_data['GT'], indiv_data, PL_Thresh):
            continue

        father = v[family[indiv_id][0]]
        mother = v[family[indiv_id][1]]
        # determine if the parents have the alternate allele 
        # so they can pass it on AND apply filters
        if filters.TDT_Parent_Filters(indiv_data, father, mother, GQ_Parent_Thresh, DP_Thresh, PL_Thresh):
            TU = numberTransmissions(indiv_data['GT'], father['GT'], mother['GT'], TU)
            
    line = v.get('CHROM') +','+ v.get('POS')+','+ v.get('SNPEFF_GENE_NAME')+','+ v.get('SNPEFF_FUNCTIONAL_CLASS')+','+ v.get('AC')+','+ str(TU[0]) +','+str(TU[1])
    #print line
    #os._exit(os.EX_OK)
    #h = hp.heap()
    #print "memory usage from doTDT"
    #print h
    return line


def numberTransmissions(kid, dad, mom, TU):
    """ Series of cases:
    Kid Dad Mom Transmissions   Untransmissions
    Ref Het Het 0               2
    Ref Ref Het 0               1
    Het Het Het 1               1
    Het Ref Het 1               0
    Het Alt Het 0               1
    Alt Het Het 2               0
    Alt Het Alt 1               0
    """
    
    # Child Father Mother: [transmissions, untransmissions]
    transmissions = {
        'homoRef het het': [0, 2],
        'homoRef homoRef het': [0, 1],
        'homoRef het homoRef': [0,1],
        'het het het': [1, 1],
        'het homoRef het': [1, 0],
        'het het homoRef': [1, 0],
        'het homoAlt het': [0, 1],
        'het het homoAlt': [0, 1],
        'homoAlt het het': [2, 0],
        'homoAlt het homoAlt': [1, 0],
        'homoAlt homoAlt het': [1, 0]
    }

    key = '{} {} {}'.format(kid, dad, mom)

    counts = transmissions.get(key)

    if not counts:
        print "ERROR {}".format(key)
        return TU
    
    # update the total number of transmissions and untransmissions
    TU = [x+y for x,y in zip(TU, counts)]
    
    return TU


def doCaseControl(v, cases, controls, vcf_indivIDs, AC_Thresh, GQ_CC_Thresh, DP_Thresh, 
                  PL_Thresh, synFlag):
    """ Analyses case/control data, counting the number of reference and alternate alleles
        case and control are dictionaries with:
            Key: individual id      Value: gender    
    """
    
    # apply filters to the line, if true look at the line otherwise move on
    if not filters.processLine(v.get('FILTER'), v.get('REF'), v.get('ALT'),
                               v.get('SNPEFF_FUNCTIONAL_CLASS'), v.get('AC'), 
                               AC_Thresh, synFlag):
        return None
    
    # count the number of ref and alt alleles in cases and controls
    caseRefs = 0
    caseAlts = 0
    controlRefs = 0
    controlAlts = 0
    
    # loop through all the individuals in the vcf        
    for indiv_id in vcf_indivIDs:
        # indiv_data is their GT:AD:DP:GQ:PL stats
        indiv_data = v.get(indiv_id)
        
        if indiv_data == None:
            continue
        
        # check if the individual is a case
        if indiv_id in cases:
            # apply filters and update counts
            caseRefs, caseAlts = ProcessCC(indiv_data, GQ_CC_Thresh, DP_Thresh, PL_Thresh, caseRefs, caseAlts)

        elif indiv_id in controls:
            # apply filters and update counts
            controlRefs, controlAlts = ProcessCC(indiv_data, GQ_CC_Thresh, DP_Thresh, PL_Thresh, controlRefs, controlAlts)
                            
    line = v.get('CHROM') +','+ v.get('POS')+','+ v.get('SNPEFF_GENE_NAME')+','+ v.get('SNPEFF_FUNCTIONAL_CLASS')+','+ v.get('AC')+','+ str(caseRefs) +','+str(caseAlts) +','+str(controlRefs)+','+str(controlAlts)
    #print line

    return line


def ProcessCC(indivAttr, GQ_Thresh, DP_Thresh, PL_Thresh, Refcount, Altcount):
    """ apply filters to the individual
        update the number of Reference and Alternative alleles
    """
    
    #apply quality filters
    if not filters.passFilters(indivAttr, GQ_Thresh, DP_Thresh):
        return Refcount, Altcount
    
    # appply PL filter to individual
    if not filters.PhredScaleFilter(indivAttr['GT'], indivAttr, PL_Thresh):
        return Refcount, Altcount 

    # homoRef has 2 reference alleles
    if indivAttr['GT'] == 'homoRef':
        Refcount = Refcount + 2
        return Refcount, Altcount 
    
    # het has 1 reference allele and 1 alternate allele
    elif indivAttr['GT'] == 'het':
        Refcount += 1
        Altcount += 1
        return Refcount, Altcount 
    
    # homoAlt has 2 alternate alleles
    elif indivAttr['GT'] == 'homoAlt':
        Altcount = Altcount + 2
        return Refcount, Altcount 
    
    # if the individual has another genotype,
    # return the input
    else:
        return Refcount, Altcount
        

if __name__ == "__main__":
    args = docopt(__doc__, version='0.1')
    #hp = hpy()
    #hp.setref()
    # GLOBAL VARIABLES
    pool = mp.Pool(processes=int(args['<cores>']), maxtasksperchild=int(args['<chunksize>']) )
    print "size of pool is ", sys.getsizeof(pool)
    # case and control are dictionaries with:
    #   Key: individual id      Value: gender
    # family is a dictionary with:
    #   Key: individual id      Value: [father id, mother id]
    case, control, familyRelations = Family.readFamily(args['<family_Relations_File>'])

    PL_Threshold = float(args['--pl'])
    DP_Threshold = float(args['--dp'])
    AC_Threshold = float(args['--ac'])    
    GQ_Kid_Threshold = float(args['--gq_Kid'])
    GQ_Parent_Threshold = float(args['--gq_Par'])
    GQ_CC_Threshold = float(args['--gq_CC'])
    synonymousFlag = args['--syn']
  
    if args['family']:
        Family.readFamily(args['<family_Relations_File>'])
    elif args['doTDT']:
        work = []
        # create a file to write to
        writer = csv.writer(open(args['<outputFile_Name>'], 'wb'))
    
        # write out the header
        writer.writerow(["CHROM", "POSITION", "GENE_NAME", "SNPEFF_EFFECT", "AC",
                     "transmitted", "untransmitted"])  
        
        print(args)

        # use a partial function so that multiple parameters
        # can be passed into map()
        partial_TDT = partial(doTDT, family=familyRelations, 
                              AC_Thresh=AC_Threshold, GQ_Kid_Thresh=GQ_Kid_Threshold, 
                              GQ_Parent_Thresh=GQ_Parent_Threshold, 
                              DP_Thresh = DP_Threshold, PL_Thresh=PL_Threshold, 
                              synFlag=synonymousFlag)
        #print "initialization"
        #print hp.heap() 
        #hp.setref()
        
        for result in pool.map_async(partial_TDT, VCF.lines(args['<vcf_File>']), 
                                       chunksize=int(args['<chunksize>']) ):
            writer.writerow([result])        
        # for v in VCF.lines(args['<vcf_File>']):
            # work.append(v)

            # if len(work) == int(args['<chunksize>']):
                # for p in mp.active_children():
                    # print "memory of active children equal"
                    # print sys.getsizeof(p)                
                # for result in pool.map_async(partial_TDT, work):
                    # writer.writerow([result])
                # #sys.exit()
                # work = []
                
            # if len(work) > 0:
                # print "memory of active children less than"
                # for p in mp.active_children():
                    # print sys.getsizeof(p)
                # for result in pool.map_async(partial_TDT, work):
                    # writer.writerow([result])
        #h = hp.heap()
        #print "memory usage from parallel"
        #print h
        pool.close()
        pool.join()
#        results = map(partial_TDT, VCF.lines(args['<vcf_File>']))
        
#        print map(partial_TDT, VCF.lines(args['<vcf_File>']))        
#        print results

        # for result in pool.imap_unordered(partial_TDT, VCF.lines(args['<vcf_File>']), 
                                       # chunksize=int(args['<chunksize>']) ):
            # writer.writerow([result])
        # #h = hp.heap()
        # #print h
        # pool.close()
        # pool.join()
#        for result in results:
        
    elif args['doCaseControl']:
        # list of individuals in the VCF
        VCF_indivIDs = find_VCFindivs(args['<vcf_File>'])
        
        # create a file to write to
        writer = csv.writer(open(args['<outputFile_Name>'], 'wb'))
    
        # write out the header
        writer.writerow(["CHROM", "POSITION", "GENE_NAME", "SNPEFF_EFFECT", "AC",
                     "caseRefs", "caseAlts", "controlRefs", "controlAlts"]) 
        
        print(args)
        
        # use a partial function so that multiple parameters can be passed into map()
        partial_CC = partial(doCaseControl, cases=case, controls=control, 
                             vcf_indivIDs=VCF_indivIDs, AC_Thresh=AC_Threshold, 
                             GQ_CC_Thresh=GQ_CC_Threshold, DP_Thresh = DP_Threshold, 
                             PL_Thresh=PL_Threshold, synFlag=synonymousFlag)
        #print map(partial_TDT, VCF.lines(args['<vcf_File>']))

        #hp.setrelheap()
        results = pool.imap_unordered(partial_CC, VCF.lines(args['<vcf_File>']), 
                                      chunksize=int(args['<chunksize>']))

        #h = hp.heap()
        #print h
        for result in results:
            print sys.getsizeof(result)
            writer.writerow([result])

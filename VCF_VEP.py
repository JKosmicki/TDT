#!/usr/bin/env python
"""
:File:         VCF_VEP.py
:OriginalAuthor:   Kamil Slowikowski <kslowikowski@fas.harvard.edu>
:Edited:       Jack Kosmicki
:Last updated: January 21, 2015

Make life with VCF files easier.

Usage example:

    >>> import VCF
    >>> variants = list(VCF.lines('file.vcf.gz'))
    >>> print [v['CHROM'] for v in variants[:3]]
    [1, 1, 1]
"""


import gzip
import sys


def parse(line, hashTable):
    """Parse a single VCF line and return a dictionary. 
    
    Parameters
    ----------
    line: a line in the VCF (string)
    hashTables: hash tables of individuals to look up in the VCF
    """

    vcfLine = {}

    # the 1st 8 defined Columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT
    FIELDS = line.split('\t')[:9]

    vcfLine['FILTER'] = FIELDS[6]

    if not vcfLine['FILTER'].startswith('PASS'):
        return None

    # ignore multi-allelics
    if ',' in FIELDS[4]:
        return None

    # load dictionary
    vcfLine['CHROM'] = FIELDS[0]
    vcfLine['POS'] = int(FIELDS[1])
    vcfLine['ID'] = FIELDS[2]
    vcfLine['REF'] = FIELDS[3]
    vcfLine['ALT'] = FIELDS[4]
    vcfLine['QUAL'] = FIELDS[5]

    format = FIELDS[8].split(':')
    format_len = len(format)

    # apply filters to the line, if true look at the line otherwise move on
    if 'AD' not in format or 'PL' not in format:
        return None

    # INFO field consists of "key1=value;key2=value;...".
    infos = FIELDS[7].split(';')

    for i, info in enumerate(infos, 1):
        # It should be "key=value".
        try:
            key, value = info.split('=')
        # But sometimes it is just "value".
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Set the value to None if there is no value.
        if not value:
            value = None
        vcfLine[key] = value

    fields = line.split('\t')[9:]       # fields are all of the individuals in the vcf.

    for indivID in hashTable.keys():
        # individuals values based on the format field
        stats = {}

        # indivAttr (individual attributes) are GT:AD:DP:GQ:PL
        indivAttr = fields[hashTable[indivID]].strip('"').split(':')

        # If the sizes don't match, the individual has missing values 
        # in the format field.
        # i.e., format = GT:AD:DP and indivAttr = ./.
        if len(indivAttr) != format_len:
            vcfLine[indivID] = None
        else:
            # create a dictionary of GT:AD:DP:GQ:PL
            for i in range(format_len):
                if ',' in indivAttr[i]:
                    indivAttr[i] = indivAttr[i].split(',')
                stats[format[i]] = indivAttr[i]
            try:
                # sometimes AD is none
                if stats['AD'] != None:
                    stats['GT'] = find_gtype(stats)
                else:
                    vcfLine[indivID] = None
                    print 'AD is none. Does this even get called?'
            except KeyError:
                sys.stderr.write('Problem with format field: {}.\n'.format(stats))

            vcfLine[indivID] = stats

    return vcfLine


def find_gtype(stats):
    """ Determine the genotype of the individual by checking both 
        the genotype and AD.

    Parameters
    ----------
    stats = format field [GT,[AD_ref, AD_alt],DP,GQ,[PL_ref, PL_het, PL_alt]]
    """

    # Throw away individuals with 0 reference and alternate reads.
    if stats['AD'][0] == '0' and stats['AD'][1] == '0':
        return None
    elif stats['GT'] in ('0/0', '0|0'):
        return 'homoRef'
    elif stats['GT'] in ('0/1', '0|1', '1/0', '1|0'):
        return 'het'
    elif stats['GT'] in ('1/1', '1|1'):
        return 'homoAlt'
    else:
        return None

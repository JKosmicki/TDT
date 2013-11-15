"""
:File: filters.py
:Author: Jack A. Kosmicki
:Last updated: 2013-10-27

File of filters for quality control of VCF files

"""


def processLine(FILTER, REF, ALT, Class, AC, AC_Thresh, synFlag):
    """Processes the variant line and runs quality checks on the
       variant
    """
    
    if not FILTER.startswith('PASS'):
        return False

    # Check for the presence of multi-allelic sites
    #      but we want to allow AT but not A,T
    if ',' in REF or REF is None:
        return False
    if ',' in ALT or ALT is None:
        return False
    
    # process synonymous variants above a certain threshold?
    if synFlag == True:
        # process synonymous variants if the allele count
        # is less than the specified threshold
        #print v.get('AC')
        if Class == 'SILENT' and float(AC) > AC_Thresh:
            #print "synonymous variant with AC > the threshold"
            #print "{0} and {1}".format(AC), AC_Thresh)
            print "got inside check AC threshold"
            return False

    return True
    

def AllelicBalance(AD):
    """ calculate the allelic balance """
    ref = float(AD[0])
    alt = float(AD[1])
    
    # allelic balance = # non-reference calls/ total # of calls
    ABratio = alt / (ref + alt)
#     try:
#         ABratio = alt / (ref + alt)
#     except ZeroDivisionError:
#         ABratio = False
        
    return ABratio

    
def PhredScaleFilter(gtype, stats, PL_Thresh):
    """ This function determines what phred scale filter to use for
        a given individual as different genotypes require
        different filters
        Input: genotype, GT:AD:DP:GQ:PL values, the PL threshold
    """
    
    if gtype == 'het':
        return PhredScaleFilter_HET(stats, PL_Thresh)
    elif gtype == 'homoRef':
        return PhredScaleFilter_HOMOREF(stats, PL_Thresh)
    elif gtype == 'homoAlt':
        return PhredScaleFilter_HOMOALT(stats, PL_Thresh)
    else:
        return False


def PhredScaleFilter_HET(stats, PL_Thresh):
    """ apply filters for the normalized Phred-scaled likelihoods for 
        AA, AB, BB genotypes where A = reference allele, B = alternate allele
        stats is the GT:AD:DP:GQ:PL scores
        True = Passed Filter
        False = Failed Filter
    """
    
    homoRef = int(stats['PL'][0])
    het = int(stats['PL'][1])
    homoAlt = int(stats['PL'][2])
    
    if homoRef <= PL_Thresh:
        return False
    elif het != 0:
        return False
    elif homoAlt <= PL_Thresh:
        return False
    else:
        return True


def PhredScaleFilter_HOMOREF(stats, PL_Thresh):
    """ apply filters for the normalized Phred-scaled likelihoods for 
        AA, AB, BB genotypes where A = reference allele,
                                   B = alternate allele
        stats is the GT:AD:DP:GQ:PL scores
        True = Passed Filter
        False = Failed Filter
    """
    
    homoRef = int(stats['PL'][0])
    het = int(stats['PL'][1])
    homoAlt = int(stats['PL'][2])
    
    if homoRef != 0:
        return False
    elif het <= PL_Thresh:
        return False
    elif homoAlt <= PL_Thresh:
        return False
    else:
        return True
        
def PhredScaleFilter_HOMOALT(stats, PL_Thresh):
    """ apply filters for the normalized Phred-scaled likelihoods for 
        AA, AB, BB genotypes where A = reference allele,
                                   B = alternate allele
        stats is the GT:AD:DP:GQ:PL scores
        True = Passed Filter
        False = Failed Filter
    """
    
    homoRef = int(stats['PL'][0])
    het = int(stats['PL'][1])
    homoAlt = int(stats['PL'][2])
    
    if homoAlt != 0:
        return False
    elif het <= PL_Thresh:
        return False
    elif homoRef <= PL_Thresh:
        return False
    else:
        return True


def depthFilter(child_dp, dad_dp, mom_dp):
    """ filter comparing child's depth with parents' depths 
        child's depth must be > 10% of the parents' combined depths
    """

    #can be altered with option
    depth_ratio = 10  

    dpRatio_threshold = (dad_dp + mom_dp) / depth_ratio

    return child_dp >= dpRatio_threshold

	
def passFilters(stats, GQ_Thresh, DP_Thresh):
    """series of filters used to remove bad calls
        stats are the individuals values (GT:AD:DP:GQ:PL)
        AD: allelic depth
        DP: approximate number of reads that passed the filter
        GQ: genotype quality
        PL: phred scale score
    """
    
    genotype = stats['GT']
    
    # ensure individual has a genotype
    if genotype == None:
        return False

    # genotype quality score must be >= DP Threshold [default: 30]
    elif float(stats['GQ']) < GQ_Thresh:
        return False

    # read depth must be >= DP threshold [default: 10]
    elif float(stats['DP']) < DP_Thresh:
        return False

    # Allele balance of between 0.3 and 0.7 for hets
    elif genotype == 'het':
        # calculate the allelic balance from the allelic depth
        ABratio = AllelicBalance(stats['AD'])
        if ABratio < 0.3 or ABratio > 0.7:
            return False
        else:
            return True

    # Allele balance of <0.1 for homozygous Reference individuals
    elif genotype == 'homoRef':
        ABratio = AllelicBalance(stats['AD'])
        if ABratio >= 0.1:
            return False
        else:
            return True

    elif genotype == 'homoAlt':
        ABratio = AllelicBalance(stats['AD'])
        if ABratio <= 0.9:
            return False
        else:
            return True
            
            
def TDT_Parent_Filters(indiv, father, mother, GQ_Thresh,
                       DP_Thresh, PL_Thresh):
    """ determine if parents have an alternate allele to pass on
        indiv = individual's subject ID
        family = dictionary containing the families
            key = child's subject ID 
            value = [Dad ID, Mom ID]
        True = parents have alternate allele and passed filters
        False = parents don't have alternate allele or did not pass filters
    """

    # 1) find the parents' data 
    try:
        if father is None or mother is None:
            return False
        else:
            # Find the parent's GT:AD:DP:GQ:PL numbers
            dad_data, mom_data = father, mother
            
            dad_gtype, mom_gtype = dad_data['GT'], mom_data['GT']
            
            # determine the status of the phred Filter which will be either 
            # True (they passed) or False (they didn't pass)
            dad_Phred_Pass = PhredScaleFilter(dad_gtype, dad_data, PL_Thresh)
            mom_Phred_Pass = PhredScaleFilter(mom_gtype, mom_data, PL_Thresh)

            # 2) determine if the parents are heterozygous
            #   At least one parent must be heterozygous
            #   so one must be False
            if dad_gtype != 'het' and mom_gtype != 'het':
                #print "one parent isn't a het"
                return False
        
            # 3) apply allelic balance filters to parents
            #    Both parents must pass all the filters
            elif not passFilters(dad_data, GQ_Thresh, DP_Thresh) or not passFilters(mom_data, GQ_Thresh, DP_Thresh):
                #print "at least one parent didn't pass the filters"
                return False
        
            # 4) apply PL filter to parents
            elif not dad_Phred_Pass or not mom_Phred_Pass:
                #print "at least one parent didn't pass the phred filter"
                return False
        
            # 5) apply a filter for the depth of the family
            elif not depthFilter(int(indiv['DP']), int(dad_data['DP']),
                                 int(mom_data['DP'])):
                #print "failed depth filter"
                return False
                
            # 6) Check for Mendelian Errors
            elif not MendelianError(indiv['GT'], dad_gtype, mom_gtype):
                #print "Mendelian Error in Chrom {0} Pos {1} for individual {2}".format(variants.get('CHROM'), variants.get('POS'), indiv)
                return False
            else:
                return True
    except KeyError:
        print "KeyError with individual {}".format(indiv)
        return False
        

def MendelianError(proband_gtype, dad_gtype, mom_gtype):
    """ Checks for Mendelian errors in which the parents
        must have transmitted an allele but did not
        True = Passed
        False = Failed
    """
    
    # child is AA and either parent is BB
    if proband_gtype == 'homoAlt' and (dad_gtype == 'homoRef' or mom_gtype == 'homoRef'):
        return False
    # child is BB and either parent is AA
    elif proband_gtype == 'homoRef' and (dad_gtype == 'homoAlt' or mom_gtype == 'homoAlt'):
        return False
    # child is AB and both parents 
    elif proband_gtype == 'het' and ((dad_gtype == 'homoAlt' and mom_gtype == 'homoAlt') or (dad_gtype == 'homoRef' and mom_gtype == 'homoRef')):
        return False
    else:
        return True

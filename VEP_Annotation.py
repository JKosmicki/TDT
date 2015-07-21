"""
:File: VEP_Annotation.py
:Author: Jack A. Kosmicki
:Last updated: 2015-07-21

Extract Variant Effect Predictor (VEP) Annotations from a vcf.

Code graciously borrowed and modified from Konrad Karczewski
"""


import numpy as np


__version__ = 1.1
__author__ = 'Jack A. Kosmicki <jkosmicki@fas.harvard.edu>'
__date__ = '7/21/2015'


# Note that this list of VEP annotations is current as of v77 with 2 included for backwards compatibility (VEP <= 75)
csq_order = ['transcript_ablation',
'splice_donor_variant',
'splice_acceptor_variant',
'stop_gained',
'frameshift_variant',
'stop_lost',
'start_lost',
'initiator_codon_variant', # deprecated
'transcript_amplification',
'inframe_insertion',
'inframe_deletion',
'missense_variant',
'protein_altering_variant',
'splice_region_variant',
'incomplete_terminal_codon_variant',
'stop_retained_variant',
'synonymous_variant',
'coding_sequence_variant',
'mature_miRNA_variant',
'5_prime_UTR_variant',
'3_prime_UTR_variant',
'non_coding_transcript_exon_variant',
'non_coding_exon_variant', # deprecated
'intron_variant',
'NMD_transcript_variant',
'non_coding_transcript_variant',
'nc_transcript_variant', # deprecated
'upstream_gene_variant',
'downstream_gene_variant',
'TFBS_ablation',
'TFBS_amplification',
'TF_binding_site_variant',
'regulatory_region_ablation',
'regulatory_region_amplification',
'feature_elongation',
'regulatory_region_variant',
'feature_truncation',
'intergenic_variant',
'']
csq_order_dict = dict(zip(csq_order, range(len(csq_order))))
rev_csq_order_dict = dict(zip(range(len(csq_order)), csq_order))


def findVariantAnnotation(v, args, vepFieldNames):
    """ Find the gene, functional class, sift, polyphen, and loftee annotations.

    Parameters
    ----------
    v: line of the vcf in the form of a hash table
    args: optional command line arguments (hash table)
    vepFieldNames: Array of what is in each column of the VEP annotation.
                   e.g., (Allele, Gene, Feature, PolyPhen, SIFT, etc.)

    Returns
    -------
    (GENE_NAME, FUNCTIONAL_CLASS, SIFT, PolyPhen2, loftee)
    """

    if args['--vep']:
        gene, anno, pph2, sift, lof = VEP_annotate(v.get('CSQ'), vepFieldNames, v['ALT'])
    elif args['--ano']:
        # str() is called on some variables as they can be NONE producing a type error.
        gene = str(v.get('SEVERE_GENE_NAME'))
        anno = str(v.get('SEVERE_IMPACT'))
        pph2 = str(v.get('SIFT'))
        sift = str(v.get('POLYPHEN'))
        lof = ''
    else:
        gene, anno, pph2, sift, lof = ('', '', '', '', '')

    return gene, anno, pph2, sift, lof


def VEP_annotate(annotation, vepFieldNames, alt):
    """Take variant annotation and extract gene name and mutation type
    ['Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'ALLELE_NUM', 'DISTANCE', 'STRAND', 'SYMBOL', 'SYMBOL_SOURCE', 'HGNC_ID', 'BIOTYPE', 'CANONICAL', 'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL', 'UNIPARC', 'SIFT', 'PolyPhen', 'EXON', 'INTRON', 'DOMAINS', 'HGVSc', 'HGVSp', 'GMAF', 'AFR_MAF', 'AMR_MAF', 'ASN_MAF', 'EUR_MAF', 'AA_MAF', 'EA_MAF', 'CLIN_SIG', 'SOMATIC', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'LoF_info', 'LoF_flags', 'LoF_filter', 'LoF']


    Parameters
    ----------
    annotation: VEP annotation in |-delimited form (string)
    vepFieldNames: Array of what is in each column of the VEP annotation.
                   e.g., (Allele, Gene, Feature, PolyPhen, SIFT, etc.)
    alt: the alternate allele (string)

    Returns
    -------
    (GENE_NAME, FUNCTIONAL_CLASS, SIFT, PolyPhen2, loftee)
    """

    if not annotation:
        return ('', '', '', '', '')

    # array with dictionaries containing the information
    annotations = [dict(zip(vepFieldNames, x.split('|'))) for x in annotation.split(',') if len(vepFieldNames) == len(x.split('|'))]

    # loop through and choose the canonical annotation
    # check that alternative allele matches
    for entry in annotations:
        if entry['Allele'] != alt:
            continue

        if entry['CANONICAL'] != 'YES':
            continue

        gene_name = entry['SYMBOL']
        functional_class = worst_csq_from_csq(entry['Consequence'])

        sift = fixProteinPrediction(entry['SIFT'])
        pph2 = fixProteinPrediction(entry['PolyPhen'])
        lof = entry['LoF']

        return (gene_name, functional_class, sift, pph2, lof)

    # If there is no canonical transcript, return worst consequence
    worst_annotation = worst_csq_with_vep(annotations)

    if worst_annotation != None:
        gene_name = worst_annotation['SYMBOL']
        functional_class = worst_annotation['major_consequence']
        sift = fixProteinPrediction(worst_annotation['SIFT'])
        pph2 = fixProteinPrediction(worst_annotation['PolyPhen'])
        lof = worst_annotation['LoF']

        return (gene_name, functional_class, sift, pph2, lof)

    return ('', '', '', '', '')


def worst_csq_with_vep(annotation_list):
    """
    Takes list of VEP annotations [{'Consequence': 'frameshift', Feature: 'ENST'}, ...]
    Returns most severe annotation (as full VEP annotation [{'Consequence': 'frameshift', Feature: 'ENST'}])
    Also tacks on worst consequence for that annotation (i.e. worst_csq_from_csq)
    :param annotation_list:
    :return worst_annotation:
    """

    if len(annotation_list) == 0:
        return None

    worst = annotation_list[0]
    for annotation in annotation_list:
        if compare_two_consequences(annotation['Consequence'], worst['Consequence']) < 0:
            worst = annotation

        elif compare_two_consequences(annotation['Consequence'], worst['Consequence']) == 0 and annotation['CANONICAL'] == 'YES':
            worst = annotation

    worst['major_consequence'] = worst_csq_from_csq(worst['Consequence'])
    return worst


def compare_two_consequences(csq1, csq2):
    if csq_order_dict[worst_csq_from_csq(csq1)] < csq_order_dict[worst_csq_from_csq(csq2)]:
        return -1

    elif csq_order_dict[worst_csq_from_csq(csq1)] == csq_order_dict[worst_csq_from_csq(csq2)]:
        return 0

    return 1


def worst_csq_from_csq(csq):
    """
    Annotations occasionally consist of two annotations joined by '&'.
    In these situations, we return the most damaging of the two annotations.
    Example: 'non_coding_exon_variant&nc_transcript_variant'
            ['non_coding_exon_variant', 'nc_transcript_variant']
            Here we return the worst annotation, i.e., 'non_coding_exon_variant'.

    Parameters
    ----------
    csq: VEP consequence (string)

    Return
    ------
    most_severe_consequence
    """

    return rev_csq_order_dict[worst_csq_index(csq.split('&'))]


def fixProteinPrediction(annotation):
    """ PolyPhen and SIFT annotations have (#) in the string and we
        want to remove them.  However if there is no sift or polyphen
        annotation, split('(')[0] throws an error.
        E.g., probably_damaging(0.99) will become probably_damaging

    Parameters
    ----------
    annotation: SIFT or PolyPhen annotation (str)

    Returns
    -------
    annotation: SIFT or PolyPhen annotation devoid of (#) (str)
    """

    if len(annotation.split('(')) != 1:
        return annotation.split('(')[0]
    else:
        return annotation


def worst_csq_index(csq_list):
    """
    Input list of consequences (e.g. ['frameshift_variant', 'missense_variant'])
    Return index of the worst annotation (In this case, index of 'frameshift_variant', so 4)
    Works well with csqs = 'non_coding_exon_variant&nc_transcript_variant' by worst_csq_index(csqs.split('&'))
    :param annnotation:
    :return most_severe_consequence_index:
    """

    return min([csq_order_dict[ann] for ann in csq_list])

"""
:File: ProcessFile.py
:Author: Jack A. Kosmicki
:Last updated: 2013-09-10

read output file from TDT.py

Usage:
    ProcessFile.py genes <gene_File> [options]
    ProcessFile.py main <TDT_output_File> [options]
    ProcessFile.py doVarTDT <TDT_output_File> <outputFile_Name> [options]
    ProcessFile.py doTU_counts <TDT_output_File> <outputFile_Name> [options]

Options:
    --genes=G_FILE  file with names of the genes you want to look at
    -h, --help     Print this message and exit.
    -v, --version  Print the version and exit.
"""

__version__ = 1.0
__author__ = 'Jack A. Kosmicki <jkosmicki@fas.harvard.edu>'
__date__ = '10/09/2013'

import sys
import csv
import math
from collections import OrderedDict
from rpy2.robjects.packages import importr
from docopt import docopt

def findColumn(header, columnNames):
    """locates the index of where each column is
        header is the first row of the file containing the column names
        it is also in the format of an array
        columnNames is an array with the names of the columns to look up"""

    columnNameIndex = {}

    for i in range(len(columnNames)):
        try:
            columnNameIndex[columnNames[i]] = header.index(columnNames[i])
        except ValueError:
            print "Could not find column {} in header".format(columnNames[i])
            print "This is the header: ", header
            sys.exit()
            
    return columnNameIndex


def updateVariants(Dictionary, name, varType, variant):
    """ updates the counts for the observed variants for a specific gene and effect 
        variants = de novo, transmitted, untransmitted, case, control, and NHLBI
        
        Input Explanation:
            Dictionary: dictionary of stuff (e.g., variant, gene, etc.)
            Name: name of whatever we are looking for (e.g., variant, gene, etc.)
            varType: SILENT, MISSENSE, NONSENSE, NONE
            variant: variant's data
    """
    
    for i in range(len(variant)):
        Dictionary.get(name)[varType][i] += float(variant[i])
    
    return Dictionary


def TDT_stat(Dictionary):
    """ calculate TDT stat
        TDT statistic is (transmitted - untransmitted)^2 / (transmitted + untransmitted)
        Dictionary in the form: Variant: [T,U]
    """
    
    for key, value in Dictionary.items():
        for key2, value2 in value.items():
            T = float(value2[5])
            U = float(value2[6])
            
            #sometimes there will be 0 transmissions and untransmissions
            try:
                # calculate the TDT statistic which is a chi-squared statistic
                TDT = (math.pow((T-U),2))/(T+U)
                Dictionary[key][key2].append(TDT)
                
                # calculate the pvalue associated with this 
                # using the PDF with 1 degree of freedom
                # we do 1 - [value] to do lower.tail = False
                TDT_pval = stats.chi2.pdf(TDT,1)
                
                Dictionary[key][key2].append(TDT_pval)
            except ZeroDivisionError:
                Dictionary[key][key2].append('N/A')
                Dictionary[key][key2].append('N/A')

    
    return Dictionary
    

def writeFile(Dictionary):
    """ write out the dictionary """
    
    # create a file to write to
    writer = csv.writer(open('Processed_TDT.csv', 'wb'))

    writer.writerow(["GENE NAME", "SNPEFF_EFFECT", "transmitted", "untransmitted", "TDT_SCORE", "TDT_PVAL"])
    
    for key, value in Dictionary.items():
        valueOfValues = value
        for key2, value2 in valueOfValues.items():
            writer.writerow([key] + [key2] + value2)

    return None
    

def readGenes(geneFile):
    """read in the gene names from a file
    the file must contain a column called 'Gene_Names' """
    
    # Dictionary of gene names
    geneDictionary = OrderedDict()
    with open(geneFile) as inFile:
        header = inFile.next()  #header of the file
        header = header.strip().split('\t')
        
        # columns is a dictionary of the column indexes
        columns = findColumn(header, ["Gene_Names"])

        for line in inFile:
            field = line.strip().split('\t')
            
            #some columns have multiple genes so we must strip the " and split by ','
            geneNames = field[columns.get('Gene_Names')].strip('"').split(',')

            for i in range(len(geneNames)):
                # if the gene is not in the dictionary
                if((geneNames[i] in geneDictionary) == False):
                    # Dictionary is in the form: [de novo LOF, de novo missense, de novo synonymous, 
                    #                             case, control, transmitted, untransmitted, NHLBI]
                    # where each element in the array of values are the counts
                    geneDictionary[geneNames[i]] = [0, 0, 0, 0, 0, 0, 0, 0]

    print "number of genes in dictionary = ", len(geneDictionary)
    return geneDictionary
    
    
def read_calcTDT(TDT_file, outputFile, args):
    """ read in a file and 
        calculate TDT significance values
    """
    
    stats = importr('stats')
    
    # create a file to write to
    writer = open(outputFile, 'wb')
    
    
    with open(TDT_file) as inFile:
        header = inFile.next()  #header of the file
        header = header.strip('\n').strip('\r\n').strip('\r')
        writer.write(header + ',TDT_SCORE,TDT_PVAL\n')
        
        header = header.strip().split(',')
        
        # if no list of genes are specified
        if(args['--genes'] is None):
            # columns is a dictionary of the column indexes
            columns = findColumn(header, ['transmitted', 'untransmitted'])

            for line in inFile:
                field = line.strip().split(',')
            
                # find the number of transmissions (T) and untransmissions (U)
                T = float(field[columns.get('transmitted')])
                U = float(field[columns.get('untransmitted')])
            
                #sometimes there will be 0 transmissions and untransmissions
                try:
                    # calculate the TDT statistic which is a chi-squared statistic
                    TDT = (math.pow((T-U),2))/(T+U)
                    field.append(TDT)
                
                    # calculate the pvalue associated with this 
                    # using the PDF with 1 degree of freedom
                    # we do 1 - [value] to do lower.tail = False
                    TDT_val = stats.pchisq(TDT,1,lower_tail=False)[0]
                
                    field.append(TDT_val)
                
                    for i in range(len(field)):
                        writer.write(str(field[i]) + ',')
                    writer.write('\n')
            
                # when there are no transmission and untransmissions
                except ZeroDivisionError:
                    field.append('N/A')
                    field.append('N/A')
                
                    for i in range(len(field)):
                        writer.write(field[i] + ',')
                    writer.write('\n')
        
        else:
            # columns is a dictionary of the column indexes
            columns = findColumn(header, ['GENE_NAME', 'transmitted', 'untransmitted'])

            geneDictionary = readGenes(args['--genes'])
            
            for line in inFile:
                field = line.strip().split(',')
                
                if field[columns.get('GENE_NAME')] in geneDictionary:
                    # find the number of transmissions (T) and untransmissions (U)
                    T = float(field[columns.get('transmitted')])
                    U = float(field[columns.get('untransmitted')])
            
                    #sometimes there will be 0 transmissions and untransmissions
                    try:
                        # calculate the TDT statistic which is a chi-squared statistic
                        TDT = (math.pow((T-U),2))/(T+U)
                        field.append(TDT)
                
                        # calculate the pvalue associated with this 
                        # using the PDF with 1 degree of freedom
                        # we do 1 - [value] to do lower.tail = False
                        TDT_val = stats.pchisq(TDT,1,lower_tail=False)[0]
                
                        field.append(TDT_val)
                
                        for i in range(len(field)):
                            writer.write(str(field[i]) + ',')
                        writer.write('\n')
            
                    # when there are no transmission and untransmissions
                    except ZeroDivisionError:
                        field.append('N/A')
                        field.append('N/A')
                
                        for i in range(len(field)):
                            writer.write(field[i] + ',')
                        writer.write('\n')

    return None
    
      
def compressGenes(TDT_file):
    """read in the TDT output file
       and compress genes by their functional class
    """
    
    # Dictionary of gene names
    geneDictionary = OrderedDict()
    
    with open(TDT_file) as inFile:
        header = inFile.next()  #header of the file
        header = header.strip('\r').strip('\r\n').strip('\n').split(',')
        results = {}
        
        columnNames = ['GENE_NAME', 'SNPEFF_EFFECT', 'transmitted', 'untransmitted']
        # columns is a dictionary of the column indexes
        columns = findColumn(header, columnNames)

        for line in inFile:
            field = line.strip().split(',')
            
            # make a dictionary of the results for each line so that it is easy to search
            for i in range(len(columnNames)):
                results[columnNames[i]] = field[columns.get(columnNames[i])]

            # make an array of variant counts
            counts = [results.get('transmitted'), results.get('untransmitted')]
            
            # if the gene is not in the dictionary
            if((field[columns.get('GENE_NAME')] in geneDictionary) == False):
                geneDictionary[results.get('GENE_NAME')] = {'SILENT':[0, 0, 0, 0, 0, 0, 0, 0], 'MISSENSE':[0, 0, 0, 0, 0, 0, 0, 0], 'NONSENSE':[0, 0, 0, 0, 0, 0, 0, 0], 'NONE':[0, 0, 0, 0, 0, 0, 0, 0]}
                geneDictionary = updateVariants(geneDictionary, results.get('GENE_NAME'), results.get('SNPEFF_EFFECT'), counts)
            else:
                geneDictionary = updateVariants(geneDictionary, results.get('GENE_NAME'), results.get('SNPEFF_EFFECT'), counts)

    # calculate TDT statistic
    geneDictionary = TDT_stat(geneDictionary)
    print "number of genes in dictionary = ", len(geneDictionary)
    return geneDictionary

    
def find_TU_counts(TDT_file, outputFile):
    """read in the TDT output file
       and compress genes by their functional class
    """
    
    # Dictionary of gene names
    geneDictionary = {}
    
    with open(TDT_file) as inFile:
        header = inFile.next()  #header of the file
        header = header.strip('\r').strip('\r\n').strip('\n').split(',')
        results = {}
        
        columnNames = ['CHROM', 'POSITION', 'GENE_NAME', 'transmitted', 'untransmitted']
        # columns is a dictionary of the column indexes
        columns = findColumn(header, columnNames)

        for line in inFile:
            field = line.strip().split(',')
            
            # make a dictionary of the results for each line so that it is easy to search
            for i in range(len(columnNames)):
                results[columnNames[i]] = field[columns.get(columnNames[i])]

            # make an array of variant counts
            counts = [results.get('transmitted'), results.get('untransmitted')]
            
            # combine transmissions and untransmissions into a single variable
            # this is in the form T-U
            TU = field[columns.get('transmitted')] + '|' + field[columns.get('untransmitted')]
            
            # if the gene is not in the dictionary
            if((field[columns.get('GENE_NAME')] in geneDictionary) == False):
                geneDictionary[results.get('GENE_NAME')] = {results.get('SNPEFF_EFFECT'):{TU:1}}
            elif results.get('SNPEFF_EFFECT') not in geneDictionary[results.get('GENE_NAME')]:
                geneDictionary[results.get('GENE_NAME')][results.get('SNPEFF_EFFECT')] = {TU:1}
            elif TU not in geneDictionary[results.get('GENE_NAME')][results.get('SNPEFF_EFFECT')]:
                geneDictionary[results.get('GENE_NAME')][results.get('SNPEFF_EFFECT')][TU] = 1
            else:
                geneDictionary[results.get('GENE_NAME')][results.get('SNPEFF_EFFECT')][TU] += 1

            
    # create a file to write to
    writer = csv.writer(open(outputFile, 'wb'))

    writer.writerow(['GENE NAME', 'SNPEFF_EFFECT', 'TU', 'counts'])
    
    # for each name in the dictionary
    for key, value in geneDictionary.items():
        snpEffect = value
        for key2, value2 in snpEffect.items():
            TU = value2
            for key3, value3 in TU.items():
                writer.writerow([key] + [key2] + [key3] + [value3])

 
def main(TDT_file):
    geneDictionary = compressGenes(TDT_file)
    writeFile(geneDictionary)
            

if __name__ == "__main__":
    args = docopt(__doc__, version='0.1')
    if args['genes']:
        readGenes(args['<gene_File>'])
    elif args['main']:
        main(args['<TDT_output_File>'])
    elif args['doVarTDT']:
        read_calcTDT(args['<TDT_output_File>'], args['<outputFile_Name>'], args)
    elif args['doTU_counts']:
        find_TU_counts(args['<TDT_output_File>'], args['<outputFile_Name>'])

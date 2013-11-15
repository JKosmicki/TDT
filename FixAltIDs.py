"""
:File: FixAltIDs.py
:Author: Jack A. Kosmicki
:Last updated: 2013-16-9

search for schizophrenia GWAS loci in an autism vcf

Usage:
    FixAltIDs.py family <family_Relations_File> [options]

Options:
    -pl, --phred  Specify PL threshold
    -h, --help    Print this message and exit.
    --version     Print the version and exit.
"""

import gzip
import sys
import schizo
import csv
from docopt import docopt
from collections import OrderedDict

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
            sys.exit()
            
    return columnNameIndex
    
    
def findAlternateIDs(family_Relations_File):
    """ read in the family relations file """
    # dictionaries for cases, controls, families
    alt_ID = {}
    
    with open(family_Relations_File) as file:
        header = file.next()  #header of the file
        header = header.strip().split(',')
        columns = findColumn(header, ["SUBJECT_ID", "Alt_ID"])

        for line in file:
            field = line.strip().split(',')
            # add individuals to their respective dictionaries
            #print field[columns.get('SUBJECT_ID')]
            #print field[columns.get('Alt_ID')]
            if(field[columns.get('Alt_ID')] != ''):
                if(field[columns.get('SUBJECT_ID')] != field[columns.get('Alt_ID')] ):
                    alt_ID[field[columns.get('Alt_ID')]] = field[columns.get('SUBJECT_ID')]

    return alt_ID

    
def readFamily(family_Relations_File):
    """ read in the family relations file """
    # dictionaries for cases, controls, families
    case = {}
    control = {}
    family = {}
    father = {}
    mother = {}

    altID = findAlternateIDs(family_Relations_File)

    with open(family_Relations_File) as file:
        header = file.next()  #header of the file
        header = header.strip().split(',')
        columns = findColumn(header, ["SUBJECT_ID", "C/C Study", "Family Role", "Father ID", "Mother ID", "Gender"])

        for line in file:
            field = line.strip().split(',')
            # add individuals to their respective dictionaries
            if(field[columns.get('C/C Study')].strip().lower() == "case" or 
               field[columns.get('C/C Study')].strip().lower() == "swedish case" ):
                case[field[columns.get("SUBJECT_ID")]] = ["case", "N/A", "N/A", "N/A", field[columns.get("Gender")]]
            if field[columns.get('C/C Study')].strip().lower() == "control":
                control[field[columns.get("SUBJECT_ID")]] = ["control", "N/A", "N/A", "N/A", field[columns.get("Gender")]]
            # elif(field[columns.get('Family Role')] == "FATHER" or
                 # field[columns.get('Family Role')] == "father" or
                 # field[columns.get('Family Role')] == "Father" or                 
                 # field[columns.get('Family Role')] == "Trio father"):
                # father[field[columns.get("SUBJECT_ID")]] = ["N/A", "father", "N/A", "N/A", field[columns.get("Gender")]]
            elif(field[columns.get('Family Role')].strip().lower() == "proband" or 
                 field[columns.get('Family Role')].strip().lower() == "trio proband"):
                # if the proband is missing the mother or father, do not include them
                if(field[columns.get("Father ID")] == "0" or field[columns.get("Father ID")] == "N/A" or 
                   field[columns.get("Mother ID")] == "0" or field[columns.get("Mother ID")] == "N/A"):
                    print "Individual {0} has missing parents".format(field[columns.get("SUBJECT_ID")])
                    continue
                else:
                    # family dictionary is in the form {child_ID} = [Dad_ID, Mom_ID]
                    father = field[columns.get("Father ID")]
                    mother = field[columns.get("Mother ID")]
                    family[field[columns.get("SUBJECT_ID")]] = ["N/A", "proband", father, mother, field[columns.get("Gender")]]
                    # if the father has an alternate ID, set the ID to the alternate ID
                    if father in altID:
                        family[field[columns.get("SUBJECT_ID")]][2] = altID.get(father)
                    # if the mother has an alternate ID, set the ID to the alternate ID
                    if mother in altID:
                        family[field[columns.get("SUBJECT_ID")]][3] = altID.get(mother)
                        
    print "number of cases in dictionary = ", len(case)
    print "number of controls in dictionary = ", len(control)
    print "number of families in dictionary = ", len(family)
    return case, control, family

    
def main(family_Relations_File):
    case, control, family = readFamily(family_Relations_File)

    writer = csv.writer(open('dict.csv', 'wb'))
    writer.writerow(["SUBJECT_ID", "C/C Study", "Family Role", "Father ID", "Mother ID", "Gender"])
    for key, value in case.items():
        writer.writerow([key] + value)        
    for key, value in control.items():
        writer.writerow([key] + value)
    for key, value in family.items():
        writer.writerow([key] + value)  
        
if __name__ == "__main__":
    args = docopt(__doc__, version='0.1')
    if args['family']:
        main(args['<family_Relations_File>'])
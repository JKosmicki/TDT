import sys

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


def readFamily(family_Relations_File):
    """ read in the family relations file """
    # dictionaries for cases, controls, families
    case = {}
    control = {}
    family = {}

    with open(family_Relations_File) as file:
        header = file.next()  #header of the file
        header = header.strip().split(',')
        
        # columns is a dictionary of the column indexes
        columns = findColumn(header, ["SUBJECT_ID", "C/C Study", "Family Role", "Father ID", "Mother ID", "Gender"])

        for line in file:
            field = line.strip().split(',')
            indiv_id = field[columns.get("SUBJECT_ID")]
            father_id = field[columns.get("Father ID")]
            mother_id = field[columns.get("Mother ID")]
            family_role = field[columns.get('Family Role')].strip().lower()
            case_control = field[columns.get('C/C Study')].strip().lower()
            parental_nulls = ["0", "N/A", ""]
            
            # add individuals to their respective dictionaries
            if(case_control in ['case', 'swedish case'] or family_role == 'partial trio proband'):
                case[indiv_id] = field[columns.get('Gender')]

            if case_control == 'control':
                control[indiv_id] = field[columns.get('Gender')]

            # if the proband is missing the mother or father, do not include them
            elif family_role in ['proband', 'trio proband']:
                if father_id in parental_nulls or mother_id in parental_nulls:
                    print 'Individual {0} has missing parents'.format(indiv_id)
                    continue
                else:
                    # family dictionary is in the form {child_ID} = [Dad_ID, Mom_ID]
                    family[indiv_id] = [father_id, mother_id]
                        
    print 'number of cases in dictionary = ', len(case)
    print 'number of controls in dictionary = ', len(control)
    print 'number of families in dictionary = ', len(family)
    return case, control, family
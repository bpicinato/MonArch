"""
this script stores ensembles from treated libs, reads raw reads from other libs
from input, runs trimmomatic on these reads, and finally compares reads from the
stored ensemble with those who passed on trimmomatic's filter 
"""

def read_annotated_ensembles(full):
    """
    reads annotated ensembles of .full file for treated libs. returns a dict of
    dicts, containing all the infos needed from the ensembles.
    """
    annotated_ensembles = dict()
    with open(full, 'r') as file:
        # read header
        file.readline()
        for r in file.readlines():
            row = r.strip().split("\t")
            annotated_ensembles[row[0]] = {'chr': str(row[1]),
                    'start': int(row[2]),
                    'end': int(row[3]),
                    'id': str(row[4]),
                    'flag': int(row[5]),
                    'strand': str(row[6]),
                    'read' : str(row[7]),
                    'trans_size' : str(row[8]),
                    'halfA' : str(row[9]),
                    'halfB' : str(row[10]),
                    'emp_junc' : str(row[11]),
                    'real_junc' : str(row[12]),
                    'trans' : str(row[13])}
    return annotated_ensembles
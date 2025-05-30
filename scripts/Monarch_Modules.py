#! /usr/bin/python3
from datetime import datetime
import sys

# module to start a new ensemble with a given junction.
def BeginNewEnsemble(Ensembles, Junction, ID):
    Ensembles[ID] = {'id': "CircEnse_{:04d}".format(ID),
                     'chr': Junction['chr'],
                     'start': Junction['start'],
                     'end': Junction['end'],
                     'strand': Junction['strand'],
                     'junctions_start_ranges': [range(Junction['start']-3,
                                                      Junction['start']+4)],
                     'junctions_end_ranges': [range(Junction['end']-3,
                                                    Junction['end']+4)],
                     # counter of reads with each flag:
                     '97': 0,
                     '98': 0,
                     '99': 0,
                     '100': 0,
                     '101': 0,
                     '102': 0,
                     '103': 0}

# ----------------------------------------------------
# functions that will be used in the GetEnsembles step
# ----------------------------------------------------

# read junction from row (returns ThisJunction)
def ReadJunction(row, transcript_seqs = False):
    if(transcript_seqs == True):
        ThisJunction = {'chr': str(row[1]),
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
    else:
        # when called by Monarch-PostProcessing, the junctions do not
        # have sequencing information
        ThisJunction = {'chr': str(row[0]),
                        'start': int(row[1]),
                        'end': int(row[2]),
                        'id': str(row[3]),
                        'flag': int(row[4]),
                        'strand': str(row[5])}
    return ThisJunction

# write current ensemble to output .full
def WriteCurrentEnsembleToFull(WhichEnsemble, WhichJunction, WhichFile, MessageNumber):
    # print
    WhichFile.write(WhichEnsemble['id']+'\t')
    for v in WhichJunction.values():
        WhichFile.write(str(v)+'\t')
    WhichFile.write('\n')
    msg = "[ "+datetime.now().strftime("%m-%d-%y %H:%M:%S")+" | Step6 ] Annotated {} ensembles"
    sys.stdout.write(msg.format(MessageNumber-1)+'\r')

# check if junction belongs to any ensemble (returns BelongsTo)
def CheckIfBelongsToAnyEnsemble(Ensembles, ThisJunction, Sctrict = False):
    # (1) Check ensembles with same strand and chr
    SameStrandAndChr = [ID
                        for ID, ens in Ensembles.items()
                        if (ens['strand'] == ThisJunction['strand'])
                        and (ens['chr'] == ThisJunction['chr'])]
    #print(SameStrandAndChr)
    # (2) Check if start of this junction is in any junctions_start_ranges
    #SameStartRange = [(ID, ens)
    #                  for ID, ens in SameStrandAndChr
    #                  for range in ens['junctions_start_ranges']
    #                  if ThisJunction['start'] in range]
    SameStartRange = []
    for ID in SameStrandAndChr:
        for range in Ensembles[ID]['junctions_start_ranges']:
            if ThisJunction['start'] in range:
                SameStartRange.append(ID)
                break
    #print(SameStartRange)
    # (3) If exists any ensemble with same start range, check if the ens also
    # has same end range
    #BelongsTo = [ID
    #             for (ID, ens) in SameStartRange
    #             for range in ens['junctions_end_ranges']
    #             if ThisJunction['end'] in range
    #             and SameStartRange]
    BelongsTo = []
    for ID in SameStartRange:
        for range in Ensembles[ID]['junctions_end_ranges']:
            if ThisJunction['end'] in range:
                BelongsTo.append(ID)
                break
    #print(BelongsTo)
    return BelongsTo

# update the coordenates of the ensemble
def UpdateEnsembleCoordenates(WhichEnsemble, WhichJunction):
    # declare some variables
    ThisStart = WhichJunction['start']
    ThisEnd = WhichJunction['end']
    ThisStartRange = range(WhichJunction['start']-3,
                           WhichJunction['start']+4)
    ThisEndRange = range(WhichJunction['end']-3,
                         WhichJunction['end']+4)

    # update ensemble start (with the minimum)
    WhichEnsemble['start'] = min(WhichEnsemble['start'], ThisStart)
    # update ensemble end (with the maximum)
    WhichEnsemble['end'] = max(WhichEnsemble['end'], ThisEnd)
    # add coordenate to start and end ranges
    WhichEnsemble['junctions_start_ranges'].append(ThisStartRange)
    WhichEnsemble['junctions_end_ranges'].append(ThisEndRange)

# anotate the ensembles from a junction
def AnotateEnsembles(Ensembles, ThisJunction, EnsembleNumber, out_full):

    # if this is the first line of the bed file, create a new ensemble
    # containing data from current junction
    if EnsembleNumber == 1:
        BeginNewEnsemble(Ensembles, ThisJunction, 1)
        # sum one in right flag
        Ensembles[1][str(ThisJunction['flag'])] += 1
        # write to full
        WriteCurrentEnsembleToFull(WhichEnsemble = Ensembles[1],
                                   WhichJunction = ThisJunction,
                                   WhichFile = out_full,
                                   MessageNumber = EnsembleNumber)
        # increase EnsembleNumber: new ensemble was made
        return True

    else:
        # else, check if belongs to any
        BelongsTo = CheckIfBelongsToAnyEnsemble(Ensembles, ThisJunction)

        # if belongs to any ensemble...
        if BelongsTo:
            # Get the first one
            ID = sorted(BelongsTo)[0]
            # Update coordenates if needed
            UpdateEnsembleCoordenates(Ensembles[ID],ThisJunction)
            # add 1 to read counter
            Ensembles[ID][str(ThisJunction['flag'])] += 1

            # write to full
            WriteCurrentEnsembleToFull(WhichEnsemble = Ensembles[ID],
                                    WhichJunction = ThisJunction,
                                    WhichFile = out_full,
                                    MessageNumber = EnsembleNumber)
            # do not increase EnsembleNumber: junction in pre-existing ensemble
            return False
        # if this junction does not belong to any, begin new ensemble
        else:
            BeginNewEnsemble(Ensembles, ThisJunction, EnsembleNumber)
            # update read counter
            Ensembles[EnsembleNumber][str(ThisJunction['flag'])] += 1
            # print
            WriteCurrentEnsembleToFull(WhichEnsemble = Ensembles[EnsembleNumber],
                                       WhichJunction = ThisJunction,
                                       WhichFile = out_full,
                                       MessageNumber = EnsembleNumber)
            # increase EnsembleNumber: new ensemble was made
            return True


def GetEnsembles(input, full, transcript_seqs=False,
                 header=False, UpdateCoordenates=False):
    with open(input, 'r') as bed, open(full, 'a+') as out_full:
        EnsembleNumber=1
        Ensembles = dict()
        if header == True:
            bed.readline()
        for r in bed.readlines():
            row = r.strip().split('\t')
            # read junction from row
            ThisJunction = ReadJunction(row, transcript_seqs)

            # Anotate Ensembles
            add_ensemble_number = AnotateEnsembles(Ensembles, ThisJunction,
                                                    EnsembleNumber,
                                                    out_full)
            if add_ensemble_number:
                EnsembleNumber += 1

    return Ensembles

def PrintEnsemblesSimple(Ensembles, ens_simple):
    with open(ens_simple, 'a') as out:
        #out.write('\t'.join(['#read', 'chr', 'start', 'end', 'strand',
        #                        'flag97', 'flag98', 'flag99', 'flag100',
        #                        'flag101', 'flag102', 'flag103', '\n']))
        for k, v in Ensembles.items():
            for field in [str(x) for x in v.values()]:
                if 'range' not in field:
                    out.write(field+'\t')
            out.write('\n')

#! /usr/bin/python3

####### USAGE #######
# python3 CountEnsemblesFromMonarch.py <treated.full> <ctrl.full> <prefix STR> <output_dir>
#####################

import pandas as pd
import os
import sys

def main():
    
    treated_ensembles = sys.argv[1]
    ctrl = sys.argv[2]
    prefix = sys.argv[3]
    out_dir = sys.argv[4]
    
    TratEnsembles = dict()

    with open(treated_ensembles, 'r') as f:
        #skipping .full header
        f.readline()
        for line in f.readlines():
            row = line.strip().split("\t")
            if row[0] in TratEnsembles.keys():

                TratEnsembles[row[0]]['coordinates'].append((int(row[2]), int(row[3])))

            else:

                TratEnsembles[row[0]] = {'id': row[0],
                                         'chr': row[1],
                                         'strand': row[6],
                                         'coordinates': [(int(row[2]), int(row[3]))]
                                        }
                
    for ens, info in TratEnsembles.items():
        info['count'] = len(info['coordinates'])
    
    # initializing output
    out_df = pd.DataFrame(columns = ['CircEnse', 'count_treat', 'count_ctrl', 'CircEnse_lib_ctrl'])

    for ensemble,info in TratEnsembles.items():
        out_df = out_df.append({'CircEnse': ensemble,
                                'count_treat': info['count'],
                                'count_ctrl': 0, 
                                'CircEnse_lib_ctrl': ''}, ignore_index = True)

    # adding library    
    with open(ctrl, 'r') as f:

        #skipping .full header
        f.readline()

        for line in f.readlines():

            row = line.strip().split("\t")
            coordinates = (int(row[2]), int(row[3]))
            circense_lib = row[0] + ','


            for ensemble, info in TratEnsembles.items():
                # checking same strand and chr; and coordinates
                if (row[1] == info['chr']) and (row[6] == info['strand']) and (coordinates in info['coordinates']):

                        out_df.loc[out_df['CircEnse'] == info['id'], 'count_ctrl'] += 1
                        out_df.loc[out_df['CircEnse'] == info['id'], 'CircEnse_lib_ctrl'] += circense_lib
                        
    out_df.to_csv(os.path.join(out_dir, prefix+'_count.tsv'), sep = '\t', index = False)
    
if __name__ == '__main__':
    main()
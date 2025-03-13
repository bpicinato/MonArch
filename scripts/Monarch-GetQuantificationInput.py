#! /usr/bin/python3
import sys
import argparse
import os.path as op
import pandas as pd

def parsing_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
prepare input for the quantification of ensembles.

this script assumes you have run the Monarch-Seq pipeline for RNase R-treated
libraries and want to find/quantify in other libraries the annotated ensembles.

this script generates the junctions.tsv file to be used as input for the
quantification''')

    parser.add_argument('--simple', metavar='FILE',
    help='path to .simple file of ensembles of RNase R treated libraries')

    parser.add_argument('--full', metavar='FILE',
    help='path to .full file of ensembles of RNase R treated libraries')

    parser.add_argument('--top', metavar='X',
    help='the value of X, where only the top X ensembles (sorted by number of reads) will be used')

    parser.add_argument('--output', metavar='DIR',
    help='path to output dir where junctions.tsv file will be gerenated')

    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    return args


def read_slice_and_write(args):
    '''
    reads the two inputs, slice it to get only the top x ensembles and then
    writes the sliced files
    '''

    # reading
    cut = int(args.top)
    simple = pd.read_csv(args.simple, sep = '\t')
    full = pd.read_csv(args.full, sep = '\t')

    # slicing
    simple['total_reads'] = simple.iloc[:,5:11].sum(axis=1)
    simple_sliced = (simple.loc[simple['total_reads'] >= cut]
                    .sort_values('total_reads', ascending=False))
    full_sliced = full[full['#CircEnse'].isin(simple_sliced['#read'])]

    # writing
    simple_name = op.basename(args.simple) + '_top' + str(args.top)
    full_name = op.basename(args.full)  + '_top' + str(args.top)

    simple_sliced.to_csv(op.join(args.output, simple_name), sep='\t', index=False)
    full_sliced.to_csv(op.join(args.output, full_name), sep='\t', index=False)


def generate_junctions_tsv(args):
    '''
    reads the sliced full file, storing the information as a dict, and then
    writes the junctions.tsv file
    '''

    ensembles = {}
    full_name = op.basename(args.full)  + '_top' + str(args.top)

    # read .full file
    with open(op.join(args.output, full_name), 'r') as input_file:
        next(input_file)
        for line in input_file.readlines():
            row = line.strip().split('\t')
            ensemble = row[0]
            junction = row[11].replace('|', '')
            if ensemble in ensembles.keys():
                ensembles[ensemble].add(junction)
            else:
                ensembles[ensemble] = set([junction])

    # write junctions.tsv
    with open(op.join(args.output, 'junctions.tsv'), 'w') as out_file:
        for ens,s in ensembles.items():
            for j in s:
                out_file.write(ens+'\t'+j+'\n')


def main():
    args = parsing_args()
    sys.stderr.write('getting the top ensembles...\n')
    read_slice_and_write(args)
    sys.stderr.write('generating junctions.tsv file...\n')
    generate_junctions_tsv(args)
    sys.stderr.write('DONE!\n')


if __name__ == '__main__':
    main()

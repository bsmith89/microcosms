#!/usr/bin/env python3
"""Parse exported data files from LC Solutions.

Usage:

parse_lcs_export.py <FILE> > output.tsv


Warning:

Be sure to export WITHOUT summary statistics.

"""
from collections import namedtuple, OrderedDict, defaultdict
import pandas as pd
import sys
import os
import argparse

COL_RENAME = OrderedDict([
#'Unnamed: 0',
#('Sample Name', 'name'),
#('Sample ID', 'extract_id'),
#('Sample Type', 'type'),
#'Level',
('Ret. Time', 'tip'),
('Area', 'area'),
('Height', 'height'),
#'ISTD Area',
#'ISTD Height',
#'Area Ratio',
#'Height Ratio',
#('Conc.', 'concentration'),
#'Std. Conc.',
#('Deviation', 'stderr_maybe'),
#'%Dev',
#'Accuracy',
#'QC Check Results',
#'Mark',
('Peak Start', 'left'),
('Peak End', 'right'),
#'T.Plate#',
#'HETP',
#'meter',
#'Tailing F.',
#'Tailing F(10%)',
#'Resolution',
#'k'',
#'Separation F.',
#'Area%',
#'Height%',
#'USP Width',
('Width(5%)',  'width_05'),
('Width(10%)', 'width_10'),
('Width(50%)', 'width_50'),
#'Relative Retention Time',
#'Statistic'
])

COLS, NAMES = zip(*COL_RENAME.items())
COLS = list(COLS)
NAMES = list(NAMES)


def main():

    p = argparse.ArgumentParser()
    p.add_argument('--constant', '-c', type=str, metavar='KEY=VALUE',
                   help='fill an entire row with one value', action='append',
                   default=[])
    p.add_argument('paths', metavar='FILE', type=str, nargs='+')
    args = p.parse_args()

    const_cols = dict(pair.split('=', 1) for pair in args.constant)

    mol_lines = OrderedDict()
    Bounds = namedtuple('Bounds', ['start', 'stop'])

    tables = []


    for path in args.paths:
        with open(path) as handle:
            start = None
            name = ''
            stop = None
            for i, line in enumerate(handle):
                if line.startswith('ID'):
                    if name:
                        assert start != None
                        stop = i - 2
                        mol_lines[name] = Bounds(start, stop)
                    start = i + 2
                if line.startswith('Name'):
                    name = line[5:-1]
            else:
                if not start:  # Never found a table
                    print("No tables found in {}".format(path), file=sys.stderr)
                    continue
                else:
                    stop = i
                    mol_lines[name] = Bounds(start, stop)

        for mol in mol_lines:
            start, stop = mol_lines[mol]
            data = pd.read_table(path,
                                skiprows=start, nrows=stop - start,
                                index_col=1, na_values='-----')
            data = data[COLS]
            data.columns = NAMES
            data.index.name = 'injection_id'
            data['molecule_id'] = mol.lower()
            tables.append(data)

    columns = list(COL_RENAME.values()) + ['molecule_id']

    if tables:
        full_data = pd.concat(tables)
    else:
        full_data = pd.DataFrame(columns=columns)
        full_data.index.name = 'injection_id'
    for key in const_cols:
        full_data[key] = const_cols[key]

    full_data.to_csv(sys.stdout, sep='\t')


if __name__ == '__main__':
    main()

#!/usr/bin/env python

################################################################################
##
## substructure_search.py
## S. Takahama (satoshi.takahama@epfl.ch)
## Nov. 2014
##
################################################################################

import os
import pandas as pd
import numpy as np
# import operator
from collections import OrderedDict
from argparse import ArgumentParser, RawTextHelpFormatter
from util import searchgroups

## define arguments
parser = ArgumentParser(description='''
============================================================
Perform substructure searches. requires 2 files: one containing SMARTS patterns
and another containing SMILES strings; creates single output of substructure
counts. Example usage:

$ python substructure_search.py --groupfile SMARTSpatterns/FTIRgroups.csv --inputfile examples/example_main.csv --outputfile example_out.csv

''',formatter_class=RawTextHelpFormatter)

## Arguments
parser.add_argument('-g','--groupfile',type=str,
                    help='file of SMARTS patterns (substructure, pattern); csv format')
parser.add_argument('-i','--inputfile',type=str,
                    help='file of SMILES strings (label, SMILES); csv format')
parser.add_argument('-o','--outputfile',type=str,default='output.csv',
                    help='output file; csv format')
parser.add_argument('-e','--export',type=str,
                    help='text file with list of compounds to select in a single column')

## Flags (on/off):
parser.add_argument('-d','--default-directory',action='store_true',help='--groupfile exists in SMARTSpatterns/')

## parse arguments
args = parser.parse_args()

# for debugging
# from collections import namedtuple
# Args = namedtuple('Args',['default_directory','groupfile','inputfile','outputfile','export'])
# args = Args(True,'FTIRgroups.csv','example_main.csv','example_out.csv',None)

if __name__=='__main__':

    ## pattern file
    if args.default_directory:
        default_directory = 'SMARTSpatterns'
        groupfile = os.path.join(os.path.dirname(__file__),
                                 default_directory,
                                 args.groupfile)
    else:
        groupfile = args.groupfile

    ## output export
    if args.export:
        with open(args.export) as f:
            export = [x for x in f]
    else:
        export = None

    ## read SMARTS patterns        
    groups = pd.read_csv(groupfile).set_index('substructure')
    if not export and 'export' in groups.columns:
        export = groups.index[groups['export'].astype('bool')]

    ## read SMILES strings
    inp = pd.read_csv(args.inputfile).set_index('compound')

    ## apply search function
    search = searchgroups(groups.pattern, export)
    output = inp.SMILES.apply(search.count)

    ## export to output
    output.to_csv(args.outputfile,index_label='compound')

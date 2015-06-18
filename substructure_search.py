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
from collections import OrderedDict
from argparse import ArgumentParser, RawTextHelpFormatter
from util import searchgroups

###_* --- Define command-line arguments
parser = ArgumentParser(description='''
============================================================
Perform substructure searches. requires 2 files: one containing SMARTS patterns
and another containing SMILES strings; creates single output of substructure
counts. Example usage:

$ python substructure_search.py --groupfile SMARTSpatterns/FTIRgroups.csv --inputfile examples/example_main.csv --outputfile example_out.csv

''',formatter_class=RawTextHelpFormatter)

###_ . Arguments
parser.add_argument('-g','--groupfile',type=str,
                    help='file of SMARTS patterns (substructure, pattern); csv format')
parser.add_argument('-i','--inputfile',type=str,
                    help='file of SMILES strings (label, SMILES); csv format')
parser.add_argument('-o','--outputfile',type=str,default='output.csv',
                    help='output file; csv format')
parser.add_argument('-e','--export',type=str,
                    help='text file with list of compounds to select in a single column')

###_ . Flags (on/off):
parser.add_argument('-d','--default-directory',action='store_true',help='--groupfile exists in SMARTSpatterns/')

if __name__=='__main__':

###_* --- Parse arguments
    
    args = parser.parse_args()

    # for debugging
    # from collections import namedtuple
    # Args = namedtuple('Args',['default_directory','groupfile','inputfile','outputfile','export'])
    # args = Args(True,'FTIRgroups.csv','example_main.csv','example_out.csv',None)

    ## pattern directory
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

###_* --- Read patterns

###_ . SMARTS       
    groups = pd.read_csv(groupfile).set_index('substructure')

###_ . SMILES
    inp = pd.read_csv(args.inputfile).set_index('compound')

###_* --- Apply search function

    if not export and 'export' in groups.columns:
        export = groups.index[groups['export'].astype('bool')]
    search = searchgroups(groups.pattern, export)
    output = inp.SMILES.apply(search.count)

###_* --- Export to output
    
    output.to_csv(args.outputfile,index_label='compound')

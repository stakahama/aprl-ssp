#!/usr/bin/env python

################################################################################
##
## substructure_search.py
## Author: Satoshi Takahama (satoshi.takahama@epfl.ch)
## Nov. 2014
##
## -----------------------------------------------------------------------------
##
## This file is part of APRL-SSP
##
## APRL-SSP is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## APRL-SSP is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with APRL-SSP.  If not, see <http://www.gnu.org/licenses/>.
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
                    help='file of SMARTS patterns (substructure, pattern, [export]); csv format')
parser.add_argument('-i','--inputfile',type=str,
                    help='file of SMILES strings (compound, SMILES); csv format')
parser.add_argument('-o','--outputfile',type=str,default='output.csv',
                    help='output file; csv format')
parser.add_argument('-e','--export',type=str,
                    help='text file with list of compounds to select in a single column')

###_ . Flags (on/off):
parser.add_argument('-d','--default-directory',action='store_true',help='--groupfile exists in SMARTSpatterns/')

if __name__=='__main__':

###_* --- Parse arguments
    
    args = parser.parse_args()

    ## pattern directory
    if args.default_directory: 
        ddirectory = os.path.join(os.path.dirname(__file__),'SMARTSpatterns')
    else:
        ddirectory = ''

    ## output export
    if args.export:
        exportfile = os.path.join(ddirectory,args.export) ## looks in same directory
        with open(exportfile) as f:
            export = [x.strip('"\'\n') for x in f]
        print 'exporting ',\
              ', '.join('{:d}: {:s}'.format(*x) for x in zip(range(1,len(export)+1),export))
    else:
        export = None

###_* --- Read patterns

###_ . SMARTS
    groupfile = os.path.join(ddirectory,args.groupfile)
    groups = pd.read_csv(groupfile).drop_duplicates().set_index('substructure')

###_ . SMILES
    inp = pd.read_csv(args.inputfile).drop_duplicates().set_index('compound')

###_* --- Apply search function

    if not export and 'export' in groups.columns:
        export = groups.index[groups['export'].astype('bool')]

    search = searchgroups(groups.pattern, export)
    output = inp.SMILES.apply(search.count)

###_* --- Export to output
    
    output.to_csv(args.outputfile,index_label='compound')

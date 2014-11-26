#!/usr/bin/env python

################################################################################
##
## substructure_search.py
## S. Takahama (satoshi.takahama@epfl.ch)
## Nov. 2014
##
################################################################################

import os
import re
import pybel
import pandas as pd
from argparse import ArgumentParser, RawTextHelpFormatter

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

## create main class/function
class searchgroups:
    ##
    def __init__(self,groups):
        self.groups = groups
        self.brackets = re.compile('\{|\}')
    ##
    def count(self,smilesstr):
        ##
        groups = self.groups
        mol = pybel.readstring('smi',smilesstr)
        abundances = pd.Series([0]*len(groups),index=groups.index)
        hasbracket = groups.map(lambda x: bool(self.brackets.search(str(x))))
        ## SMARTS search
        for key in groups.index[~hasbracket]:
            abundances[key] = len(pybel.Smarts(groups[key]).findall(mol))
        ## evaluate expressions
        for key in groups.index[hasbracket]:
            abundances[key] = eval(groups[key].format(**abundances))
        ##
        return abundances

# for debugging
# from collections import namedtuple
# Args = namedtuple('Args',['default_directory','groupfile','inputfile','outputfile','export'])
# args = Args(True,'FTIRgroups.csv','example_main.csv','example_out.csv',None)

if __name__=='__main__':

    if args.default_directory:
        default_directory = 'SMARTSpatterns'
        groupfile = os.path.join(os.path.dirname(__file__),
                                 default_directory,
                                 args.groupfile)
    else:
        groupfile = args.groupfile

    ## read SMARTS patterns        
    groups = pd.read_csv(groupfile).set_index('substructure')

    ## read SMILES strings
    inp = pd.read_csv(args.inputfile).set_index('compound')

    ## apply search function
    output = inp.SMILES.apply(searchgroups(groups.pattern).count)

    ## export output
    select = args.export if args.export else \
             (groups.index[groups['export'].astype('bool')] \
              if 'export' in groups.columns else \
              [True]*len(output.columns))
    output[select].to_csv(args.outputfile,index_label='compound')

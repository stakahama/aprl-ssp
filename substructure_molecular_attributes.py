#!/usr/bin/env python

################################################################################
##
## substructure_atoms.py
## S. Takahama (satoshi.takahama@epfl.ch)
## Nov. 2014
##
################################################################################

import os
import re
import pybel
import pandas as pd
from operator import add, itemgetter
from collections import OrderedDict
from argparse import ArgumentParser, RawTextHelpFormatter

## define arguments
parser = ArgumentParser(description='''
============================================================
Perform substructure searches. requires 2 files: one containing SMARTS patterns
and another containing SMILES strings; creates single output of substructure
counts. Example usage:

$ python substructure_atoms.py --attributes molwt FTIRgroups_foratoms.csv --inputfile examples_2/example_main.csv --outputfile examples_2/example_out.csv

''',formatter_class=RawTextHelpFormatter)

## Arguments
parser.add_argument('-a','--attributes',type=str,
                    help='comma-separated list of attributes')
parser.add_argument('-i','--inputfile',type=str,
                    help='file of SMILES strings (label, SMILES); csv format')
parser.add_argument('-o','--outputfile',type=str,default='output.csv',
                    help='output file; csv format')

## parse arguments
args = parser.parse_args()

## create main class/function

## clean = lambda varStr: re.sub('\W|^(?=\d)','_', varStr)

class queryattr:
    ##
    def __init__(self,attributes):
        self.attributes = attributes

    ##
    def getattributes(self,smilesstr):
        ##
        attributes = self.attributes
        mol = pybel.readstring('smi',smilesstr)
        return pd.Series([getattr(mol,a) for a in attributes], index=attributes)
    ##

# for debugging
# from collections import namedtuple
# Args = namedtuple('Args',['default_directory','groupfile','inputfile','outputfile','export'])
# args = Args(True,'FTIRgroups_foratoms.csv','example_main.csv','example_out.csv',None)

if __name__=='__main__':

    ## read SMILES strings
    inp = pd.read_csv(args.inputfile).set_index('compound')

    ## apply search function
    query = queryattr(args.attributes.split(','))
    output = inp.SMILES.apply(query.getattributes)

    ## export output
    output.to_csv(args.outputfile,index_label='compound')

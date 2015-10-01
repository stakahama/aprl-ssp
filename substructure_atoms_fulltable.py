#!/usr/bin/env python

################################################################################
##
## substructure_atoms_fulltable.py
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
## along with APRL-SS.  If not, see <http://www.gnu.org/licenses/>.
##
################################################################################

import os
import re
import pybel
import pandas as pd
from collections import OrderedDict
from argparse import ArgumentParser, RawTextHelpFormatter
from util import searchgroups

###_* --- Define command-line arguments
parser = ArgumentParser(description='''
============================================================
Perform substructure searches. requires 2 files: one containing SMARTS patterns
and another containing SMILES strings; creates single output of substructure
counts. Example usage:

$ python substructure_atoms.py -d --groupfile FTIRgroups_foratoms.csv --inputfile examples_2/example_main.csv --outputfile examples_2/example_out.csv

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

    ## for debugging
    ## args = parser.parse_args('-g ../../SMARTSpatterns/FTIRgroups_foratoms_complete.csv -i ../../data/inputs/testcompounds.csv -o ../../example_out.csv'.split())

    ## pattern directory
    if args.default_directory:
        default_directory = 'SMARTSpatterns'
        groupfile = os.path.join(os.path.dirname(__file__),
                                 default_directory,
                                 args.groupfile)
    else:
        groupfile = args.groupfile

###_* --- Read patterns
## added .drop_duplicates() # 2015.08.13

###_ . SMARTS        
    groups = pd.read_csv(groupfile).drop_duplicates().set_index('substructure')

###_ . SMILES
    inp = pd.read_csv(args.inputfile).drop_duplicates().set_index('compound')

###_* --- Apply search function
    
    ## query tables
    search = searchgroups(groups.pattern,groups.export)
    dflist = []
    masslist = []
    for smiles in inp.SMILES.unique():
        indextable, masstable = search.matchatoms(smiles)
        indextable['SMILES'] = smiles
        dflist.append(indextable)
        masslist.append(masstable)
    master = pd.merge(inp.reset_index(),pd.concat(dflist),
                      on='SMILES',how='outer')
    atomicmass = pd.DataFrame(list(reduce(set.union,masslist)),
                              columns=['atomtype','atomicmass']).set_index('atomtype')

    ## produce same output as previous code
    ##   this method is better than alternative ('alt method') below
    ##   as we want to subset non-null groups
    ##   in the subtable rather than the full table
    def countatoms(df):
        return len(df['atom'].ix[df['group'].notnull()].unique())
    grouped = master.groupby(['compound','type'])
    counts = grouped.apply(countatoms).reset_index(name='count')
    output = counts.pivot_table(index='compound',columns='type',values='count').ix[inp.index]
    output.fillna(0,inplace=True)

    # #alt method
    # parameters = {
    #     'index':'compound',
    #     'columns':'type',
    #     'values':'atom',
    #     'aggfunc':lambda x: len(x.unique())
    #     }
    # output = master.ix[master['group'].notnull(),['compound','type','atom']].pivot_table(**parameters)
    # output.fillna(0,inplace=True)

###_* --- Export to output
    
    extension = os.path.splitext(args.outputfile)[1]
    output.to_csv(args.outputfile,index_label='compound')
    atomicmass.to_csv(args.outputfile.replace(extension,'_atomicmass'+extension),
                      index_label='atom')
    for var in ['atom','match']:
        master[var] = master[var].map('{:.0f}'.format)
    del master['SMILES'] # 2015.08.13
    master.to_csv(args.outputfile.replace(extension,'_fulltable'+extension),index=False)

import os
import re
import pybel
import pandas as pd
from collections import OrderedDict
from argparse import ArgumentParser, RawTextHelpFormatter
from util import searchgroups

## define arguments
parser = ArgumentParser(description='''
============================================================
Perform substructure searches. requires 2 files: one containing SMARTS patterns
and another containing SMILES strings; creates single output of substructure
counts. Example usage:

$ python substructure_atoms.py -d --groupfile FTIRgroups_foratoms.csv --inputfile examples_2/example_main.csv --outputfile examples_2/example_out.csv

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

## for debugging
## args = parser.parse_args('-g ../../SMARTSpatterns/FTIRgroups_foratoms_complete.csv -i ../../data/inputs/testcompounds.csv -o ../../example_out.csv'.split())

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
##
# groups = master[['compound','type','atom']].groupby(['compound','type'])
# counts = groups.aggregate(lambda x: len(x.unique())).reset_index()
# output = counts.pivot_table(index='compound',columns='type',values='atom')
parameters = {
    'index':'compound',
    'columns':'type',
    'values':'atom',
    'aggfunc':lambda x: len(x.unique())
    }
output = master.ix[master['group'].notnull(),['compound','type','atom']].pivot_table(**parameters)
output.fillna(0,inplace=True)

## export output
extension = os.path.splitext(args.outputfile)[1]
output.to_csv(args.outputfile,index_label='compound')
atomicmass.to_csv(args.outputfile.replace(extension,'_atomicmass'+extension),
                  index_label='atom')
master['atom'] = master['atom'].map('{:.0f}'.format)
master.to_csv(args.outputfile.replace(extension,'_fulltable'+extension),index=False)

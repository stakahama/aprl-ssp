#!/usr/bin/env python

################################################################################
##
## substructure_generate_fulltable.py
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
Perform substructure searches at atom level. Requires 2 files: one containing SMARTS patterns
and another containing SMILES strings; creates several output files of substructure
counts, atom counts, atomic masses of matched atoms, table of atoms and matched groups. Example usage:

$ python substructure_generate_fulltable.py -d -g MCMgroups.csv -i apinenemech.csv -o apinenemech

''',formatter_class=RawTextHelpFormatter)

###_ . Arguments

parser.add_argument('-g','--groupfile',type=str,
                    help='file of SMARTS patterns (substructure, pattern); csv format')
parser.add_argument('-i','--inputfile',type=str,
                    help='file of SMILES strings (label, SMILES); csv format')
parser.add_argument('-o','--outputprefix',type=str,default='output',
                    help='output prefix')
parser.add_argument('-e','--export',type=str,
                    help='text file with list of compounds to select in a single column')

###_ . Flags (on/off):
parser.add_argument('-d','--default-directory',action='store_true',
                    help='--groupfile exists in SMARTSpatterns/')

if __name__=='__main__':

###_* --- Parse arguments
    
    args = parser.parse_args()

    ## for debugging
    ## args = parser.parse_args('-g SMARTSpatterns/FTIRextra.csv -i examples/example_main.csv -o output'.split())

    ## pattern directory
    if args.default_directory: 
        ddirectory = os.path.join(os.path.dirname(__file__),'SMARTSpatterns')
    else:
        ddirectory = ''

###_* --- Read patterns
## added .drop_duplicates() # 2015.08.13

###_ . SMARTS
    groupfile = os.path.join(ddirectory,args.groupfile)        
    groups = pd.read_csv(groupfile).drop_duplicates().set_index('substructure')

###_ . SMILES
    inp = pd.read_csv(args.inputfile).drop_duplicates().set_index('compound')[['SMILES']]

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
    del master['SMILES']    
    atomicmass = pd.DataFrame(list(reduce(set.union,masslist)),
                              columns=['atomtype','atomicmass']).set_index('atomtype')

    ## create tables of counts
    def docount(var,req_uniq='group'):
        def fn(df):
            return len(df[var].ix[df[req_uniq].notnull()].unique())
        return fn
    
    param = {'atoms':('type','atom'), 'groups':('group','match')}
    tables = {}
    for k in param.keys():
        grouped = master.groupby(['compound',param[k][0]])
        counts = grouped.apply(docount(param[k][1])).reset_index(name='count')
        widef = counts.pivot_table(index='compound',columns=param[k][0],values='count').ix[inp.index]
        widef.fillna(0,inplace=True)
        tables[k] = widef
    
###_* --- Export to output

    def float2int(df,columns=None):
        if not columns:
            columns = df.columns
        for var in columns:
            df[var] = df[var].map('{:.0f}'.format)
        return df

    outpath = os.path.dirname(args.outputprefix)
    prefix = os.path.basename(args.outputprefix)    
    extension = '.csv'
    filename = '{}_{}{}'
    outputfiles = {
        k:os.path.join(outpath,filename.format(prefix,k,extension)) for k in
        ['atomcounts','groupcounts','atomicmass','atomfulltable']
        }
    
    float2int(tables['atoms']).to_csv(outputfiles['atomcounts'],index_label='compound')
    float2int(tables['groups']).to_csv(outputfiles['groupcounts'],index_label='compound')
    atomicmass.to_csv(outputfiles['atomicmass'],index_label='atom')
    master = float2int(master,['atom','match'])
    master.to_csv(outputfiles['atomfulltable'],index=False)

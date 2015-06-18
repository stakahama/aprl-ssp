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

## create main class/function

## clean = lambda varStr: re.sub('\W|^(?=\d)','_', varStr)

class searchgroups:
    ##
    def __init__(self,groups,include):
        self.groups = groups
        self.include = include
        self.brackets = re.compile('(\\{[^{}]+\\})')
        self.atomicmass = []
    ##
    def matchatoms(self,smilesstr):
        ##
        groups = self.groups
        include = self.include
        brackets = self.brackets
        mol = pybel.readstring('smi',smilesstr)
        pairs = OrderedDict(zip(groups.index,[None]*len(groups)))
        hasbracket = groups.map(lambda x: bool(brackets.search(str(x))))
        ## SMARTS search
        for key in groups.index[~hasbracket]:
            pairs[key] = set(pybel.Smarts(groups[key]).findall(mol))
        ## evaluate expressions
        for key in groups.index[hasbracket]:
            eq = groups[key]
            variables = brackets.findall(eq)
            br = {}
            for idx,var in enumerate(variables):
                newvar = '_{}_'.format(idx)
                eq = eq.replace(var,newvar)
                br[newvar] = pairs[var.strip('{}')]
            pairs[key] = eval(eq,br)
        ##
        usepairs = OrderedDict([(k,v) for (k,v) in pairs.items() if include.ix[k]])
        allpairs = reduce(set.union,usepairs.values())
        allatomidx = reduce(add,map(list,allpairs))
        atomicmass = [(atom.type,atom.atomicmass) for atom in mol.atoms if atom.idx in allatomidx]
        self.atomicmass += atomicmass
        allatoms = map(itemgetter(0),atomicmass)
        atomtypes = set(allatoms)
        return pd.Series(map(allatoms.count,atomtypes),index=atomtypes)
        # mass = sum([atom.atomicmass for atom in mol.atoms if atom.idx in allatomidx])
        # return mass
    ##
    def get_atomicmass(self):
        return set(self.atomicmass)
    

# for debugging
# from collections import namedtuple
# Args = namedtuple('Args',['default_directory','groupfile','inputfile','outputfile','export'])
# args = Args(True,'FTIRgroups_foratoms.csv','example_main.csv','example_out.csv',None)

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
    sg = searchgroups(groups.pattern,groups.export)
    output = inp.SMILES.apply(sg.matchatoms)
    atomicmass = sg.get_atomicmass()

    output = output.fillna(0)
    atomicmass = pd.DataFrame(map(itemgetter(1),atomicmass),
                              index=map(itemgetter(0),atomicmass),
                              columns=['atomicmass'])

    ## export output
    output.to_csv(args.outputfile,index_label='compound')
    extension = os.path.splitext(args.outputfile)[1]
    atomicmass.to_csv(args.outputfile.replace(extension,'_atomicmass'+extension),
                      index_label='atom')

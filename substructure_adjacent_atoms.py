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
import openbabel
import pandas as pd
from collections import OrderedDict
from argparse import ArgumentParser, RawTextHelpFormatter
from util import searchgroups
## import igraph ## didn't need

###_* --- Define command-line arguments

parser = ArgumentParser(description='''
============================================================
Find adjacent atoms and bond order for each atom.
Example usage:

$ python substructure_adjacent_atoms.py -i apinenemech.csv -o apinenemech_adjacent_atoms.csv

''',formatter_class=RawTextHelpFormatter)

###_ . Arguments

parser.add_argument('-i','--inputfile',type=str,
                    help='file of SMILES strings (label, SMILES); csv format')
parser.add_argument('-o','--outputfile',type=str,default='output',
                    help='output file name')


if __name__=='__main__':

###_* --- Parse arguments

    args = parser.parse_args()

    ## for debugging
    ## args = parser.parse_args('-i examples/example_main.csv -o output.csv'.split())

###_* --- Read SMILES file
    inp = pd.read_csv(args.inputfile).drop_duplicates().set_index('compound')[['SMILES']]

## http://openbabel.org/docs/dev/UseTheLibrary/Python_PybelAPI.html
## http://openbabel.org/docs/dev/UseTheLibrary/PythonExamples.html
## http://openbabel.org/dev-api/classOpenBabel_1_1OBAtom.shtml

## pyatom.idx == obatom.GetIdx()
## pyatom.idx != obatom.GetIndex()

    edgelist = []
    for compound in inp.index:
        mol = pybel.readstring('smi', inp.SMILES.ix[compound])
        mol.addh()
        for pyatom in mol.atoms:
            obatom = pyatom.OBAtom
            idx1 = obatom.GetIdx()
            atype1 = obatom.GetType()
            for neighbor in openbabel.OBAtomAtomIter(obatom):
                idx2 = neighbor.GetIdx()
                atype2 = neighbor.GetType()
                bond = obatom.GetBond(neighbor)
                bondorder = bond.GetBondOrder()
                edgelist.append((compound, idx1, idx2, atype1, atype2, bondorder))

    edgeframe = pd.DataFrame(edgelist, columns=['compound', 'atom1', 'atom2', 'atom1_type', 'atom2_type', 'bondorder'])

    edgeframe.to_csv(args.outputfile, index=False)

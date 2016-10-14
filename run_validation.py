#!/usr/bin/env python

################################################################################
##
## run_validation.py
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

import sys
import os
from subprocess import call

###_* ===== class definition

class command_args:

    def __init__(self,execpath):
        self.args = {
            'search':os.path.join(execpath,'substructure_search.py'),
            'generate':os.path.join(execpath,'substructure_generate_fulltable.py'),
            'validate_atoms':os.path.join(execpath,'validation_atoms.R'),
            'validate_groups':os.path.join(execpath,'validation_groups.R')
            }

    def add(self,**kwargs):
        return dict(self.args.items()+kwargs.items())

###_* ===== commands =====

commands = {}

commands['atoms'] = {
    'allatoms':'''
python {search} \
    -i {prefix}.csv \
    -o {prefix}_commonatoms.csv \
    -d -g common_atoms.csv
''',
    'matchedatoms':'''
python {generate} \
    -i {prefix}.csv \
    -o {prefix}_MCMgroups \
    -d -g MCMgroups.csv
''',
    'validation':'''
Rscript --vanilla {validate_atoms} \
    {prefix}_MCMgroups_atomfulltable.csv \
    {prefix}_commonatoms.csv \
    {prefix}
'''
    }

commands['groups'] = {
    'matchedgroups':'''
python {search} \
    -i compounds_{groupfile}.csv \
    -o matched_{groupfile}.csv \
    -d -g {groupfile}.csv
''',
    'validation':'''
Rscript --vanilla {validate_groups} \
    compounds_{groupfile}.csv \
    matched_{groupfile}.csv \
    {groupfile}
'''
    }

###_* ===== execute =====

execpath = os.path.dirname(os.path.abspath(__file__))

datapath = sys.argv[-1]
if not os.path.exists(datapath):
    datapath = os.path.join(execpath,'validation')

argslist = command_args(execpath)

here = os.getcwd()
os.chdir(datapath)

## read file lists
with open('filelist_atoms.csv') as f:
    prefixlist = f.read().strip().split()

with open('filelist_groups.csv') as f:
    grouplist = f.read().strip().split()

## perform atom-level validation
for prefix in prefixlist:
    print 'validating', prefix, '...'
    newargs = argslist.add(prefix=prefix)
    call(commands['atoms']['allatoms'].format(**newargs).split())
    call(commands['atoms']['matchedatoms'].format(**newargs).split())
    call(commands['atoms']['validation'].format(**newargs).split())

## perform group-level validation (less specific)
for groupfile in grouplist:
    print 'validating', groupfile, '...'
    newargs = argslist.add(groupfile=groupfile)
    export = ' -e SIMPOLexportlist.csv' if groupfile=='SIMPOLgroups' else ''
    call((commands['groups']['matchedgroups'].format(**newargs)+export).split())
    call(commands['groups']['validation'].format(**newargs).split())

os.chdir(here)

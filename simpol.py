#!/usr/bin/env python


################################################################################
##
## simpol.py
## Author: Satoshi Takahama (satoshi.takahama@epfl.ch)
## Jul. 2017
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

## Implements SIMPOL.1 from Pankow and Asher, doi:10.5194/acp-8-2773-2008, 2008.
## This code is nearly identical to https://github.com/stakahama/aprl-kpp-gp/blob/master/util/simpol2.py, except class is renamed.

from collections import OrderedDict
import pandas as pd
import numpy as np
import sys
from argparse import ArgumentParser, RawTextHelpFormatter

## -----------------------------------------------------------------------------

## Define command-line arguments

parser = ArgumentParser(description='''
============================================================
Calculate properties using SIMPOL.1. Example usage:

$ simpol.py -i apinene_SIMPOLgroups.csv -o apinene_props_298.csv -t 298.15

''',formatter_class=RawTextHelpFormatter)

## Arguments

parser.add_argument('-i','--inputfile',type=str,help='fragment table')
parser.add_argument('-o','--outfile',type=str,default='properties.csv',help='')
parser.add_argument('-t','--temp',type=float,default=298.15,help='')

## -----------------------------------------------------------------------------

class SIMPOL1:

    def __init__(self):
        table5 = [
            ('groups',[
                'zeroeth group',
                'carbon number','carbon number on the acid-side of an amide (asa)',
                'aromatic ring','non-aromatic ring','C=C (non-aromatic)',
                'C=C-C=O in non-aromatic ring',
                'hydroxyl (alkyl)','aldehyde','ketone','carboxylic acid',
                'ester','ether','ether (alicyclic)','ether, aromatic',
                'nitrate','nitro','aromatic hydroxyl',
                'amine, primary','amine, secondary','amine, tertiary','amine, aromatic',
                'amide, primary','amide, secondary','amide, tertiary',
                'carbonylperoxynitrate','peroxide','hydroperoxide','carbonylperoxyacid',
                'nitrophenol','nitroester'
                ]),
            ('k',[
                0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
                21,22,23,24,25,26,27,28,29,30
                ]),
            ('coefficient',[
                'b0','b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','b11',
                'b12','b13','b14','b15','b16','b17','b18','b19','b20','b21','b22',
                'b23','b24','b25','b26','b27','b28','b29','b30'
                ]),
            ('footnote comment',[
                'a','b','c','d','e','f','g','h','i','j','k','L','m','m','m','n',
                'o','p','q','q','q','q','c','c','c','r','r','r','r','p','L'
                ]),
            ('Bk,1',[
                -4.26938E+02,-4.11248E+02,-1.46442E+02,3.50262E+01,
                -8.72770E+01,5.73335E+00,-2.61268E+02,-7.25373E+02,
                -7.29501E+02,-1.37456E+01,-7.98796E+02,-3.93345E+02,
                -1.44334E+02,4.05265E+01,-7.07406E+01,-7.83648E+02,
                -5.63872E+02,-4.53961E+02,3.71375E+01,-5.03710E+02,
                -3.59763E+01,-6.09432E+02,-1.02367E+02,-1.93802E+03,
                -5.26919E+00,-2.84042E+02,1.50093E+02,-2.03387E+01,
                -8.38064E+02,-5.27934E+01,-1.61520E+03
                ]),
            ('Bk,2',[
                2.89223E-01,8.96919E-01,1.54528E+00,-9.20839E-01,1.78059E+00,
                1.69764E-02,-7.63282E-01,8.26326E-01,9.86017E-01,5.23486E-01,
                -1.09436E+00,-9.51778E-01,-1.85617E+00,-2.43780E+00,-1.06674E+00,
                -1.03439E+00,-7.18416E-01,-3.26105E-01,-2.66753E+00,1.04092E+00,
                -4.08458E-01,1.50436E+00,-7.16253E-01,6.48262E-01,3.06435E-01,
                -6.25424E-01,2.39875E-02,-5.48718E+00,-1.09600E+00,-4.63689E-01,
                9.01669E-01
                ]),
            ('Bk,3',[
                4.42057E-03,-2.48607E-03,1.71021E-03,2.24399E-03,-3.07187E-03,
                -6.28957E-04,-1.68213E-03,2.50957E-03,-2.92664E-03,5.50298E-04,
                5.24132E-03,-2.19071E-03,-2.37491E-05,3.60133E-03,3.73104E-03,
                -1.07148E-03,2.63016E-03,-1.39780E-04,1.01483E-03,-4.12746E-03,
                1.67264E-03,-9.09024E-04,-2.90670E-04,1.73245E-03,3.25397E-03,
                -8.22474E-04,-3.37969E-03,8.39075E-03,-4.24385E-04,-5.11647E-03,
                1.44536E-03
                ]),
            ('Bk,4',[
                2.92846E-01,1.40312E-01,-2.78291E-01,-9.36300E-02,-1.04341E-01,
                7.55434E-03,2.89038E-01,-2.32304E-01,1.78077E-01,-2.76950E-01,
                -2.28040E-01,3.05843E-01,2.88290E-01,9.86422E-02,-1.44003E-01,
                3.15535E-01,-4.99470E-02,-3.93916E-02,2.14233E-01,1.82790E-01,
                -9.98919E-02,-1.35495E-01,-5.88556E-01,3.47940E-02,-6.81506E-01,
                -8.80549E-02,1.52789E-02,1.07884E-01,2.81812E-01,3.84965E-01,
                2.66889E-01
                ])
            ]

        self.table = pd.DataFrame.from_dict(OrderedDict(table5))

    def get_groupnames(self):
        return self.table['groups']

    def p0_atm(self,nuk,temp):
        ## calculates vapor pressure from vector of abundances (nuk) and
        ##   temperature (temp in K)
        Bk1 = self.table['Bk,1']
        Bk2 = self.table['Bk,2']
        Bk3 = self.table['Bk,3']
        Bk4 = self.table['Bk,4']
        bk = Bk1/temp + Bk2 + Bk3*temp + Bk4*np.log(temp)
        return 10**np.sum(nuk*bk)

    def deltaHvap_kJpermol(self,nuk,temp):
        ## calculates enthalpy of vaporization from vector of abundances (nuk) and
        ##   temperature (temp in K)
        ## R is the gas constant
        ## db is db(T)/d(1/T)
        Bk1 = self.table['Bk,1']
        Bk2 = self.table['Bk,2']
        Bk3 = self.table['Bk,3']
        Bk4 = self.table['Bk,4']
        R = 8.3144622e-3 ## kJ K^-1 mol^-1
        dbT = Bk1 - Bk3*temp**2 - Bk4*temp
        return -2.303*R*np.sum(nuk*dbT)

    def calc_properties(self,temp,index=None):
        compounds = self.compounds
        p0_atm = self.p0_atm
        deltaHvap_kJpermol = self.deltaHvap_kJpermol
        if not index:
            index = compounds.index
        props = pd.DataFrame(np.tile(np.nan,(len(index),2)),
                             index=index,
                             columns=['p0','DeltaH'])
        for idx in index:
            abundances = compounds.ix[idx,self.table['groups']].values # preserve order of abundances
            props.ix[idx,'p0'] = p0_atm(abundances,temp)
            props.ix[idx,'DeltaH'] = deltaHvap_kJpermol(abundances,temp)
        return props

    def read_compounds(self,filename):
        compounds = pd.read_csv(filename)
        if compounds.shape[0] > compounds.drop_duplicates().shape[0]:
            sys.exit('USER ERROR: duplicate rows in input matrix')
        self.compounds = compounds.set_index('compound')

if __name__=='__main__':

    ## parse arguments
    args = parser.parse_args()

    ## estimate properties
    simp = SIMPOL1()
    simp.read_compounds(args.inputfile)
    out = simp.calc_properties(args.temp)
    out.reset_index().to_csv(args.outfile,index=False)

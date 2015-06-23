
################################################################################
##
## util.py
## S. Takahama (satoshi.takahama@epfl.ch)
## Nov. 2014
##
################################################################################

import pybel
import re
import sys
import pandas as pd
import numpy as np
from collections import OrderedDict
from operator import add
from functional import compose
import os
import functional
import functools

compose = functools.partial(functools.reduce, functional.compose)

userfile = os.path.join(os.path.dirname(__file__),'userfn.py')
if os.path.exists(userfile):
    execfile(userfile,globals())

class searchgroups:

    def __init__(self,groups, include=None):
        self.groups = groups
        self.include = include
        self.evalkw = re.compile("^eval[ ]([^{}]+)")     # added 17.06.2015
        self.brackets = re.compile("(?<!')\{([^{}]*)\}") # negative lookahead added 17.06.2015
        self.quotes = re.compile("'\{([^{}]*)\}")        # added 17.06.2015
        self.unquote = lambda x: x.replace("}","}'")
        
    def commonattr(self):
        return (self.groups,
                self.include,
                self.evalkw,
                self.brackets,
                self.quotes,
                self.unquote)

    def matchedpatt(self,groups):
        return (groups.map(compose((bool,self.evalkw.search,str))),
                groups.map(compose((bool,self.brackets.search,str))),
                groups.map(compose((bool,self.quotes.search,str))))

    def count(self,smilesstr):
        ##
        groups, include, evalkw, brackets, quotes, unquote = self.commonattr()
        haskw, hasbracket, hasquote = self.matchedpatt(groups)
        if include is None:
            include = [True]*len(groups)
        ##
        mol = pybel.readstring('smi',smilesstr)
        mol.addh()
        abundances = pd.Series([np.nan]*len(groups),index=groups.index)
        ## SMARTS search
        for key in groups.index[~haskw & ~hasbracket & ~hasquote]:
            abundances[key] = len(pybel.Smarts(groups[key]).findall(mol))
        ## evaluate eval keyword
        for key in groups.index[haskw]: # untested
            abundances[key] = round(eval(evalkw.search(groups[key]).group(1)))
        ## evaluated quoted expressions
        for key in groups.index[hasquote]: # untested
            abundances[key] = round(eval(unquote(groups[key]).format(**groups)))
        ## evaluate expressions
        orderedexpr = self.__orderexpr(groups,hasbracket,brackets)            
        for key in orderedexpr: #groups.index[hasbracket]:
            abundances[key] = round(eval(groups[key].format(**abundances)))
        ##
        return abundances[include].astype(int)

    def matchatoms(self,smilesstr):
        ##
        groups, include, evalkw, brackets, quotes, unquote = self.commonattr()
        haskw, hasbracket, hasquote = self.matchedpatt(groups)        
        ##
        mol = pybel.readstring('smi',smilesstr)
        mol.addh()
        tups = OrderedDict(zip(groups.index,[None]*len(groups)))
        ## SMARTS search
        for key in groups.index[~haskw & ~hasbracket & ~hasquote]:
            tups[key] = set(pybel.Smarts(groups[key]).findall(mol))
        ## evaluate eval keyword
        for key in groups.index[haskw]:
            tups[key] = eval(evalkw.search(groups[key]).group(1))
        ## evaluate quoted expressions
        for key in groups.index[hasquote]: # untested
            tups[key] = self.__substitute(mol,groups[key],quotes,groups)
        ## evaluate expressions
        orderedexpr = self.__orderexpr(groups,hasbracket,brackets)
        for key in orderedexpr: #groups.index[hasbracket]
            tups[key] = self.__substitute(mol,groups[key],brackets,tups)
        usetups = OrderedDict([(k,v) for (k,v) in tups.items() if include.ix[k]])
        alltups = reduce(set.union,usetups.values())
        allatoms = reduce(add,map(list,alltups),[])
        atomicmass = set([(atom.type,atom.atomicmass) for atom in mol.atoms
                          if atom.idx in allatoms])
        ##
        idxlabel = 'atom'
        atomtype = pd.DataFrame([(atom.idx,atom.type) for atom in mol.atoms],
                                columns=[idxlabel,'type']).set_index(idxlabel)
        matched = self.__atomtable(atomtype,usetups)
        ##
        return (matched, atomicmass)
        # self.atomicmass += atomicmass
        # allatoms = map(itemgetter(0),atomicmass)
        # atomtypes = set(allatoms)
        # return pd.Series(map(allatoms.count,atomtypes),index=atomtypes)
        # mass = sum([atom.atomicmass for atom in mol.atoms if atom.idx in allatomidx])
        # return mass

    @staticmethod
    def __substitute(mol,eq,pattern,env):
        variables = pattern.findall(eq)
        br = {}
        for idx,var in enumerate(variables):
            newvar = '_{}_'.format(idx)
            br[newvar] = env[var]
            eq = eq.replace('{{{}}}'.format(var),newvar)
        return eval(eq,br)

    @staticmethod
    def __orderexpr(groups,hasbracket,brackets):
        ##
        prepended = lambda x,y: [x]+y
        ##
        computable = set(prepended('',groups.index.tolist()))
        computed = set(prepended('',groups.index[~hasbracket].tolist()))
        remaining = groups.index[hasbracket].tolist()
        tokensdict = {grp:set(brackets.findall(groups[grp])) for grp in remaining}
        ##
        maxiter = len(groups)*10 # to break out if stuck for some reason
        ordered = []
        i = 0
        while len(remaining) > 0:
            grp = remaining.pop(0)
            tokens = tokensdict[grp]
            if not tokens.issubset(computable):
                undefined = ','.join(list(tokens-computable))
                sys.exit('"{}" uncomputable: "{}" undefined'.format(grp, undefined))
            ##
            if tokens.issubset(computed):
                computed = computed.union([grp])
                ordered.append(grp)
            else:
                remaining.append(grp)
            ##
            i += 1
            if i > maxiter:
                print 'remaining:', ','.join(remaining)
                sys.exit('exceeded maximum number of iterations {:d}'.format(maxiter))
        return ordered

    @staticmethod
    def __atomtable(atomtype,tuplist):
        ## create a table from atomtypes and matched items
        idxlabel = atomtype.index.name
        atomtype.reset_index(inplace=True)
        columns = atomtype.columns.tolist()+['match','group']
        dflist = []
        i = 1
        for group, tup in tuplist.items():
            if len(tup)==0:
                continue
            allatoms = reduce(add,map(list,tup),[])
            groupdf = pd.DataFrame(allatoms,columns=[idxlabel])
            groupdf['match'] = i
            groupdf['group'] = group
            dflist.append(groupdf)
            i += 1
        if len(dflist)==0:
            out = pd.DataFrame(columns=columns)
        else:
            allgroupdf = pd.concat(dflist)
            out = atomtype.merge(allgroupdf,on=idxlabel,how='outer')[columns]
        return out

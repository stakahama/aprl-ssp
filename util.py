import pybel
import re
import sys
import pandas as pd
import numpy as np
from collections import OrderedDict
from operator import add

class searchgroups:

    def __init__(self,groups, include=None):
        self.groups = groups
        self.include = include
        self.brackets = re.compile('\{([^{}]*)\}')

    def count(self,smilesstr):
        ##
        groups = self.groups
        include = self.include
        brackets = self.brackets
        hasbracket = groups.map(lambda x: bool(brackets.search(str(x))))
        orderedexpr = self.__orderexpr(groups,brackets,hasbracket)
        if include is None:
            include = [True]*len(groups)
        ##
        mol = pybel.readstring('smi',smilesstr)
        abundances = pd.Series([np.nan]*len(groups),index=groups.index)
        ## SMARTS search
        for key in groups.index[~hasbracket]:
            abundances[key] = len(pybel.Smarts(groups[key]).findall(mol))
        ## evaluate expressions
        for key in orderedexpr: #groups.index[hasbracket]:
            abundances[key] = round(eval(groups[key].format(**abundances)))
        ##
        return abundances[include].astype(int)

    def matchatoms(self,smilesstr):
        ##
        groups = self.groups
        include = self.include
        brackets = self.brackets
        hasbracket = groups.map(lambda x: bool(brackets.search(str(x))))
        orderedexpr = self.__orderexpr(groups,brackets,hasbracket)
        ## (if not include)?
        ##
        mol = pybel.readstring('smi',smilesstr)
        pairs = OrderedDict(zip(groups.index,[None]*len(groups)))
        ## SMARTS search
        for key in groups.index[~hasbracket]:
            pairs[key] = set(pybel.Smarts(groups[key]).findall(mol))
        ## evaluate expressions
        for key in orderedexpr: #groups.index[hasbracket]
            eq = groups[key]
            variables = brackets.findall(eq)
            br = {}
            for idx,var in enumerate(variables):
                newvar = '_{}_'.format(idx)
                br[newvar] = pairs[var]
                eq = eq.replace('{{{}}}'.format(var),newvar)
            pairs[key] = eval(eq,br)
        ##
        usepairs = OrderedDict([(k,v) for (k,v) in pairs.items() if include.ix[k]])
        matched = self.__atomtable(mol,usepairs)
        allpairs = reduce(set.union,usepairs.values())
        allatoms = reduce(add,map(list,allpairs),[])
        atomicmass = set([(atom.type,atom.atomicmass) for atom in mol.atoms
                          if atom.idx in allatoms])
        return (matched, atomicmass)
        # self.atomicmass += atomicmass
        # allatoms = map(itemgetter(0),atomicmass)
        # atomtypes = set(allatoms)
        # return pd.Series(map(allatoms.count,atomtypes),index=atomtypes)
        # mass = sum([atom.atomicmass for atom in mol.atoms if atom.idx in allatomidx])
        # return mass

    @staticmethod
    def __orderexpr(groups,brackets,hasbracket):
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
    def __atomtable(mol,pairlist):
        idxlabel = 'atom'
        columns = [idxlabel,'type','group']
        atomtype = pd.DataFrame([(atom.idx,atom.type) for atom in mol.atoms],
                                columns=[idxlabel,'type'])
        dflist = []
        for group, pairs in pairlist.items():
            if len(pairs)==0:
                continue
            allatoms = reduce(add,map(list,pairs),[])
            groupdf = pd.DataFrame(allatoms,columns=[idxlabel])
            groupdf['group'] = group
            dflist.append(groupdf)
        if len(dflist)==0:
            out = pd.DataFrame(columns=columns)
        else:
            allgroupdf = pd.concat(dflist)
            out = atomtype.merge(allgroupdf,on=idxlabel,how='outer')[columns]
        return out

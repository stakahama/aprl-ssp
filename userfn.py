#!/usr/bin/env python

################################################################################
##
## userfn.py
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
## along with APRL-SSP.  If not, see <http://www.gnu.org/licenses/>.
##
################################################################################

from itertools import chain

def count_aromatic_rings(mol):
    return len(atoms_aromatic_rings(mol))

def count_nonaromatic_rings(mol):
    return len(atoms_nonaromatic_rings(mol))

def count_nitrophenols(mol, phenol, nitro):
    ## counts nitrophenol
    return len(atoms_nitrophenols(mol, phenol, nitro))

def atoms_aromatic_rings(mol):
    return set([tuple(a.idx for a in mol.atoms if ring.IsInRing(a.idx))
                for ring in mol.sssr if ring.IsAromatic()])

def atoms_nonaromatic_rings(mol):
    return set([tuple(a.idx for a in mol.atoms if ring.IsInRing(a.idx))
                for ring in mol.sssr if not ring.IsAromatic()])

def atoms_nitrophenols(mol, phenol, nitro):
    ## returns phenols for which nitro groups are found in same ring
    def is_part_of_ring(ring, group):
        ring_atoms = atom_indices_ring(r, mol) #not the most efficient
        return len([idx for idx in group if idx in ring_atoms]) > 0
    def atom_indices_ring(ring,mol):
        return [a.idx for a in mol.atoms if ring.IsInRing(a.idx)]
    def atom_indices_group(groups):
        return list(chain.from_iterable(groups))
    _phenol = pybel.Smarts(phenol).findall(mol)
    _nitro = pybel.Smarts(nitro).findall(mol)
    _rings = [ring for ring in mol.sssr if ring.IsAromatic()]
    atomlist_ring = []        # list of rings
    atomlist_nitrophenol = [] # list of (nitro)phenol groups
    for r in _rings:
        part = {'phenol':[],'nitro':[]}
        for x in _phenol:
            if not is_part_of_ring(r, x):
                continue
            part['phenol'].append(x)
            for y in _nitro:
                if not is_part_of_ring(r, y):
                    continue
                part['nitro'].append(y)
        if part['phenol'] and part['nitro']:
            atomlist_ring.append(set(atom_indices_ring(r, mol) +
                                     atom_indices_group(part['phenol']) +
                                     atom_indices_group(part['nitro'])))
            atomlist_nitrophenol += part['phenol']
    # returning of atomlist_ring is optional
    # but the phenol groups are what are really counted so this is what is returned
    return atomlist_nitrophenol

#!/usr/bin/env python

################################################################################
##
## spiderquery.py
## S. Takahama (satoshi.takahama@epfl.ch)
## Nov. 2014
##
################################################################################

import os
from argparse import ArgumentParser, RawTextHelpFormatter
from collections import OrderedDict, namedtuple

## define arguments
parser = ArgumentParser(description='''
============================================================
Perform query of ChemSpider database. Returns table of SMILES strings,
molecular weight in CSV file; optionally save to python pickle database. 

* [Query -> CSV]
* [Query -> DB -> CSV]
* [Query -> DB]
* [DB -> CSV]
============================================================

Usage examples:

--------------------
[Query to CSV] Create CSV files (example_main.csv, example_alternates.csv)
containing matched compounds in intersect(compounds.csv, database):

$ spiderquery.py -p example -i compounds.csv

With token:

$ spiderquery.py -t ./token.txt  -p example -i compounds.csv

Further examples will assume token file is specified as ~/.chemspidertoken.

--------------------
[Query to DB] Create shelve database file (.db) containing matched compounds in compounds.csv:

$ spiderquery.py -d -p example -i compounds.csv

--------------------
[Query to DB and CSV] Create pickle DB and CSV files containing matched compounds
in compounds.csv:

$ spiderquery.py -b -p example -i compounds.csv

--------------------
[DB to CSV]

Create csv files containing all compounds in pickle database:

$ spiderquery.py -D -p example -i compounds.csv

Create csv files containing all compounds in pickle database:

$ spiderquery.py -D -p example

''',formatter_class=RawTextHelpFormatter)

## Flags (on/off):
parser.add_argument('-D','--from-db',action='store_true')
parser.add_argument('-d','--export-db-only',action='store_true')
parser.add_argument('-b','--export-db-csv',action='store_true')

## Store value:
parser.add_argument('-p','--prefix',type=str,default='queryresults',
                    help='output will generate {name}_main.csv and {name}_alternates.csv and/or database {name}_p.db.')
parser.add_argument('-i','--inputfile',type=str,help='file of compound names (optional with -D flag)')
parser.add_argument('-t','--tokenfile',type=str,default='~/.chemspidertoken',
                    help='file of chemspider token (optional)')

## parse arguments
args = parser.parse_args()

## conditional loading of modules
if not args.from_db:
    from chemspipy import ChemSpider

if not args.export_db_only:
    import pandas as pd

if args.from_db or args.export_db_only or args.export_db_csv:
    import shelve

##==============================================================================

class spiderquery:

    def __init__(self, csp=None, dbname=None,index_label='compound'):
        self.fields = OrderedDict([
            ('common_name','common_name',),
            ('CSID','csid'),
            ('SMILES','smiles'),
            ('molecular_weight','molecular_weight'),
            ])
        self.null = OrderedDict([('first',None),('rest','')])
        self.index_label = index_label
        self.csp = csp
        self.dbname = dbname

    def search(self,cmpd):
        """
        Main search function on individual character string cmpd
        """
        matched = self.null
        try:
            results = self.csp.search(cmpd)
            if results:
                matched['first'] = results[0]
                matched['rest'] = ';'.join(['{:d}:{:s}'.format(r.csid,r.common_name.encode('utf-8'))
                                            for r in results[1:] if r])
        except KeyError:
            pass
        return matched

    def search2db(self,compounds):
        """
        Query -> DB
        """
        db = shelve.open(self.dbname)
        master = db['master'] if 'master' in db.keys() else {}
        for cmpd in compounds:
            if cmpd in master.keys():
                continue
            results = self.search(cmpd)
            results['first'] = self.localized(results['first'])
            csid = str(results['first'].csid)
            master[cmpd] = csid
            if csid not in db.keys():
                db[csid] = results
        db['master'] = master
        db.close()

    def localized(self,c):
        """
        Called by self.search2db()
        """
        if c:
            for f in self.fields.values():
                getattr(c,f)
            # c.mol_3d
        return c

    def db_getter(self,db):
        """
        Returns getter function for db
        called by self.db2table()
        """
        def db_get(cmpd):
            if cmpd in db['master'].keys():
                return db[db['master'][cmpd]]
            else:
                return OrderedDict([('first',None),('rest','')])
        return db_get

    def db2table(self,compounds=None):
        """
        DB -> CSV
        """
        db = shelve.open(self.dbname)
        if not compounds:
            compounds = list(set(db.keys())-set('master'))
        out = self.results2table(compounds,self.db_getter(db))
        db.close()
        return out

    def search2table(self,compounds):
        """
        Query -> CSV
        """
        out = self.results2table(compounds,self.search)
        return out

    def results2table(self,compounds,retrievefn):
        """
        called by self.search2table() and self.db2table()
        """
        contents = []
        alternates = []
        for cmpd in compounds:
            results = retrievefn(cmpd)
            if results['first']:
                contents.append([getattr(results['first'],str(f)) for f in self.fields.values()])
            else: # no matches
                contents.append([pd.np.nan]*len(self.fields.items()))
                continue
            ## many matches
            if results['rest']!='':
                alternates.append((cmpd,results['rest']))
        contents = pd.DataFrame(contents,index=compounds,
                                columns=self.fields.keys())
        contents.index.name = self.index_label
        if len(alternates)==0:
            alternates = None
        alternates = pd.DataFrame(alternates,columns=[self.index_label,'CSID:common_name']).set_index(self.index_label)
        Tables = namedtuple('Tables',['main','alternates'])
        return Tables(main=contents,alternates=alternates)        

    def save_tables(self,tables,prefix):
        """
        export tables.main and tables.alternatives to
          {prefix}_main.csv and {prefix}_alternatives.csv, respectively
        """
        tables.main.to_csv(prefix+'_main.csv',index_label=self.index_label)
        tables.alternates.to_csv(prefix+'_alternates.csv',index_label=self.index_label)

##==============================================================================

# for debugging
# Args = namedtuple('Args',['tokenfile','inputfile','prefix',
#                           'from_db','export_db_only','export_db_csv'])
# args = Args(tokenfile='~/.chemspidertoken',inputfile='compounds.csv',prefix='example',
#             from_db=None,export_db_only=None,export_db_csv=None)
        
if __name__ == '__main__':

    ## ==================== set up chemspider ====================

    if not args.from_db:
        with open(os.path.expanduser(args.tokenfile)) as f:
            csp = ChemSpider(f.read().strip())
    else:
        csp = None

    spq = spiderquery(csp,args.prefix+'_p')
    
    ## ==================== list of compounds ====================
    
    if args.inputfile:
        with open(args.inputfile) as f:
            compounds = []
            for line in f:
                compounds.append(line.strip())
    else:
        compounds = None

    ## ==================== query/save ====================

    if args.from_db:
        tables = spq.db2table(compounds)
    else:
        if args.export_db_only:
            spq.search2db(compounds)
        elif args.export_db_csv:
            spq.search2db(compounds)
            tables = spq.db2table(compounds)
        else:
            tables = spq.search2table(compounds)            

    ## ==================== export ====================
        
    if 'tables' in globals():
        spq.save_tables(tables,args.prefix)

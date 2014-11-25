Substructure search program

[TOC]

Introduction
===

This program is made up of two parts:

* spiderquery.py: Queries the ChemSpider database for SMILES strings and other properties of a molecule.
* substructure\_search.py: Uses the Open Babel chemoinformatics tool to find number of instances of a substructure (specified by SMARTS pattern) occurring in a molecule (specified by SMILES pattern).

Both programs require python. Anaconda distribution is recommended as it comes with the pandas library. Package installation can be carried through with the pip package manager:

```sh
$ pip install packagename
```

"substructure_search.py" additionally requires OpenBabel (link to installation guide included below).

spiderquery.py
===

### Setup

Required python packages:

* chemspipy (search from chemspider)

### Arguments

Query to CSV is the default. Flags to change this behavior:

* `-d`: output to database
* `-b`: output to both database and CSV file
* `-D`: search compounds in database rather than ChemSpider

Three main arguments:

* `-p`: value of `PREFIX`. Name of prefix to database or CSV tables. When `-D` is present, this also indicates the name of database to be read in. If unspecified, will default to "queryresults".
* `-i`: value of `INPUTFILE` (optional). Name of file containing compounds. Optional if `-D`; otherwise required.
* `-t`: value of `TOKENFILE` (optional). Name of file containing ChemSpider token if not "~/.chemspidertoken".

### Examples

Input/output included in folder "examples/".

* Query to CSV:

    ```
    $ spiderquery.py -p example -i compounds.csv
	```

    This will create files called "examples\_main.csv" and "examples\_alternates.csv".

* Query to database only:

    ```
    $ spiderquery.py -d -p example -i compounds.csv
	```

    This will create a file called "examples\_p.db".  ("_p" stands for python "pickle object created by shelve module)

* Query to both database and CSV:

    ```
    $ spiderquery.py -b -p example -i compounds.csv
	```

    This will create files called "examples\_main.csv", "examples\_alternates.csv", and "examples\_db.p".

* Database to CSV:

    ```
    $ spiderquery.py -D -p example -i compounds.csv
	```

    This will read compounds from file called "examples\_p.db" and create files called "examples\_main.csv", "examples\_alternates.csv". if `-i compounds.csv` is omitted, all compounds in the database will be included.

substructure_search.py
===

### Setup

Required external programs:

* Open Babel (http://openbabel.org/wiki/Category:Installation)

Required python packages:

* pybel (interface to Open Babel)

### Arguments

Main arguments:

* `-g`: value of `GROUPFILE`. Name of file which contains columns {substructure, pattern}. An additional column, "export", consisting of 0/1 values indicating whether this substructure should be included in `OUTPUTFILE` is allowed.
* `-i`: value of `INPUTFILE`. Name of file which contains columns {compound, SMILES}.
* `-o`: value of `OUTPUTFILE`. Name of file which contains matrix of compound x substructure.
* `-e`: value of `EXPORT` (optional). Name of file which contains list of substructures to include in `OUTPUTFILE`. Overrides "export" column if present in `GROUPFILE`.

Flags:

* `-d`: When present, indicates that `GROUPFILE` exists in the subdirectory, `SMARTSpatterns/`, distributed with "substructure_search.py".

### Examples

```
$ python substructure_search.py -d -g FTIRgroups.csv \
  -i example_main.csv -o example_out.csv
```

### Additional information

Patterns specified in GROUPFILE can be derived from a combination of SMARTS patterns. For instance, `ester_all` is defined as `"[CX3,CX3H1](=O)[OX2H0][#6]"`. `nitroester` is defined as `"[#6][OX2H0][CX3,CX3H1](=O)[C;$(C[N+](=O)[O-]),$(CC[N+](=O)[O-]),$(CCC[N+](=O)[O-]),$(CCCC[N+](=O)[O-]),$(CCCCC[N+](=O)[O-])]"`. `ester` can be defined as `{ester_all}-{nitroester}`. When present, such custom patterns are computed after all the SMARTS patterns have been matched and counted. 


APRL-SSP (Substructure Search Program)
===

[TOC]

## Introduction

APRL-SS (__APRL__ __S__ubstructure __S__earch __P__rogram) is made up of three primary units:

* "spider\_query.py": Queries the ChemSpider database for SMILES strings and other properties of a molecule.
* "substructure\_search.py": Uses the Open Babel chemoinformatics tool to find number of instances of a substructure (specified by SMARTS pattern) occurring in a molecule (specified by SMILES pattern).
* "substructure\_atoms\_fulltable.py": Uses the Open Babel chemoinformatics tool to find atoms associated with a substructure (specified by SMARTS pattern) occurring in a molecule (specified by SMILES pattern).

Its application is described by

> Ruggeri, G. and Takahama, S. Technical Note: 
> Use of chemoinformatic tools to enumerate functional groups in molecules 
> for organic aerosol characterization,  _submitted_, 2015.

The program is released under the Gnu Public License (GPLv3). Please contact the corresponding author, Satoshi Takahama (satoshi.takahama@epfl.ch), with any bug reports or questions.

### Setup

You may want to add the directory of the program path so that the two python scripts can be accessed anywhere. In Linux/Mac (include in `~/.bash_profile` to make this available across various terminal sessions):

```sh
$ export PATH=/path/to/aprl-structsearch:$PATH
```

In Windows (use setx instead of set to make available across various terminal sessions):

```dos
set PATH=%PATH%;Z:/path/to/aprl-substructsearch
```
See [this page](http://www.computerhope.com/issues/ch000549.htm) for GUI instructions (there may be better sites).

In examples below, it is assumed that the path has been added and the scripts are made executable (`chmod +x scriptname.py`). Otherwise, each command must be issued as `python /path/to/scriptname.py`. All examples are run in the "examples/" subdirectory provided assuming the program directory has been included in the PATH variable.

### Python

All scripts require python. Anaconda distribution is recommended as it comes with the pandas library (assumed installed in further instructions). Python is installed with Mac OS X but a separate installation of the scientific python stack is recommended. Package installation can be carried through with the pip package manager:

```sh
$ pip install packagename
```

"substructure\_search.py" and "substructure\_atoms\_fulltable.py" additionally require Open Babel (link to installation guide included below).

## Program scripts

### --- spider\_query.py ---

#### Setup

Required python packages:

* chemspipy (search from chemspider). For more information, see [http://chemspipy.readthedocs.org/](http://chemspipy.readthedocs.org/).

#### Arguments

Query to CSV is the default. Flags to change this behavior:

* `-D`: search compounds in database rather than ChemSpider
* `-d`: output to database
* `-b`: output to both database and CSV file

Three main arguments:

* `-p`: value of `PREFIX`. Name of prefix to database or CSV tables. When `-D` is present, this also indicates the name of database to be read in. If unspecified, will default to "queryresults".
* `-i`: value of `INPUTFILE` (optional). Name of file containing compounds. Optional if `-D`; otherwise required.
* `-t`: value of `TOKENFILE` (optional). Name of file containing ChemSpider token if not "~/.chemspidertoken". See next section.

#### Examples

Input/output included in folder "examples/".

* Query to CSV:

    ```
    $ spider_query.py -p example -i compounds.csv
	```

    This will create files called "examples\_main.csv" and "examples\_alternates.csv".

* Query to database only:

    ```
    $ spider_query.py -d -p example -i compounds.csv
	```

    This will create a file called "examples\_p.db".  ("\_p" stands for python "pickle object created by shelve module)

* Query to both database and CSV:

    ```
    $ spider_query.py -b -p example -i compounds.csv
	```

    This will create files called "examples\_main.csv", "examples\_alternates.csv", and "examples\_db.p".

* Database to CSV:

    ```
    $ spider_query.py -D -p example -i compounds.csv
	```

    This will read compounds from file called "examples\_p.db" and create files called "examples\_main.csv", "examples\_alternates.csv". if `-i compounds.csv` is omitted, all compounds in the database will be included.

#### Token

If the [ChemSpider token](http://chemspipy.readthedocs.org/en/latest/guide/intro.html#obtaining-a-security-token) is not written to "~/.chemspidertoken", specify its location with the `-t` argument. For instance,

```
$ spider_query.py  -t /path/to/token.txt -p example -i compounds.csv
```

### --- substructure\_search.py ---

#### Setup

Required external programs:

* Open Babel (http://openbabel.org/wiki/Category:Installation)

Required python packages:

* pybel (interface to Open Babel)

#### Arguments

Main arguments:

* `-g`: value of `GROUPFILE`. Name of file which contains columns {substructure, pattern}. An additional column, "export", consisting of 0/1 values indicating whether this substructure should be included in `OUTPUTFILE` is allowed.
* `-i`: value of `INPUTFILE`. Name of file which contains columns {compound, SMILES}.
* `-o`: value of `OUTPUTFILE`. Name of file which contains matrix of compound x substructure.
* `-e`: value of `EXPORT` (optional). Name of file which contains list of substructures to include in `OUTPUTFILE`. Overrides "export" column if present in `GROUPFILE`.

Flags:

* `-d`: When present, indicates that `GROUPFILE` exists in the subdirectory, `SMARTSpatterns/`, distributed with "substructure_search.py".

#### Examples

```
$ substructure_search.py -d -g FTIRgroups.csv \
  -i example_main.csv -o example_out.csv
```

### --- substructure\_atoms\_fulltable.py ---

Similar to substructure\_search.py but returns atoms associated with each fragment. Technically, this information can be used to enumerate fragments as done by substructure\_search.py; but the syntax can be more complex (on account of using set operations) when specifying fragments as a function of other SMARTS patterns.

#### Setup

Required external programs:

* Open Babel (http://openbabel.org/wiki/Category:Installation)

Required python packages:

* pybel (interface to Open Babel)

#### Arguments

Main arguments:

* `-g`: value of `GROUPFILE`. Name of file which contains columns {substructure, pattern}. An additional column, "export", consisting of 0/1 values indicating whether this substructure should be included in `OUTPUTFILE` is allowed.
* `-i`: value of `INPUTFILE`. Name of file which contains columns {compound, SMILES}.
* `-o`: value of `OUTPUTFILE`. Name of file which contains matrix of compound x substructure.
* `-e`: value of `EXPORT` (optional). Name of file which contains list of substructures to include in `OUTPUTFILE`. Overrides "export" column if present in `GROUPFILE`. (*currently not implemented*)

Flags:

* `-d`: When present, indicates that `GROUPFILE` exists in the subdirectory, `SMARTSpatterns/`, distributed with "substructure_search.py".

#### Examples

```
$ substructure_atoms_fulltable.py -d -g FTIRgroups_foratoms.csv \
  -i example_main.csv -o example_out.csv
```

## Pattern file

Patterns specified in GROUPFILE can be derived from a combination of SMARTS patterns using set operations. For instance, `ester_all` is defined as `"[CX3,CX3H1](=O)[OX2H0][#6]"`. `nitroester` is defined as `"[#6][OX2H0][CX3,CX3H1](=O)[C;$(C[N+](=O)[O-]),$(CC[N+](=O)[O-]),$(CCC[N+](=O)[O-]),$(CCCC[N+](=O)[O-]),$(CCCCC[N+](=O)[O-])]"`. `ester` can be defined as `{ester_all}-{nitroester}`. When present, such custom patterns are computed after all the SMARTS patterns have been matched and counted. Additionally, functions can be provided by the user. In current implementation, functions would presumably use OpenBabel methods.

### Permissible entries:

* SMARTS pattern.
* Bracketed expression. E.g., `{ester_all}-{nitroester}`.
* Quoted expression. Substituted patterns in `'{}` are not evaluated before passing to function. E.g., `count_nitrophenol(molecule,'{phenol},'{ester})`.
* Keyword `eval` denotes another type of python expression, but can also be used to preface expressions containing brackets (redundant). E.g., `eval 0`, `eval count_aromatic_rings(molecule)`.
* Keyword `molecule` denotes the pybel Molecule class instance (useful for passing to user-defined functions). E.g., `eval count_aromatic_rings(molecule)`.

If only SMARTS patterns are used, the same pattern file can be provided to "substructure\_search.py" and "substructure\_search\_fulltable.py". When additional expressions are supplied, they must be changed such that arithmetic operations are used for matched entries in "substructure\_search.py" and set operations for "substructure\_search\_fulltable.py".

### User-supplied functions

User can define functions to be called for evaluation. These functions should be contained in a file called "userfn.py." These functions will generally accept as its first argument the pybel Molecule object upon which operations are to be performed. In the pattern file, the argument to the function should be provided as `molecule` (lowercase) as this is the object whose value will be evaluated.

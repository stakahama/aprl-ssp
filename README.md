APRL-SSP (Substructure Search Program) [![DOI](https://zenodo.org/badge/19334/stakahama/aprl-ssp.svg)](https://zenodo.org/badge/latestdoi/19334/stakahama/aprl-ssp)
===

## Introduction

APRL-SSP (**APRL** **S**ubstructure **S**earch **P**rogram) is made up of several primary units:

* "spider\_query.py": Query the ChemSpider database for SMILES strings and other properties of a molecule.
* "substructure\_search.py": Use the Open Babel chemoinformatics tool to find number of instances of a substructure (specified by SMARTS pattern) occurring in a molecule (specified by SMILES pattern).
* "substructure\_generate\_fulltable.py": Use the Open Babel chemoinformatics tool to find atoms associated with a substructure (specified by SMARTS pattern) occurring in a molecule (specified by SMILES pattern).
* "generate\_carbontypes.py": Define carbon types from associated functional groups (adapted from aprl-carbontypes).
* "validation\_triple.py": Evaluate a set of patterns for specificy and completeness (adapted from aprl-carbontypes).
* "simpol.py": Estimate compound vapor pressures using SIMPOL.1.

Its application is described by

> Ruggeri, G. and Takahama, S.: Technical Note: Development of chemoinformatic tools to enumerate functional groups in molecules for organic aerosol characterization, *Atmos. Chem. Phys.*, 16, 4401-4422, doi:10.5194/acp-16-4401-2016, 2016.

> Takahama, S. and Ruggeri, G.: Technical Note: Relating functional group measurements to carbon types for improved model–measurement comparisons of organic aerosol composition, *Atmos. Chem. Phys.*, 17, 4433–4450, doi:10.5194/acp-17-4433-2017, 2017.

The program is released under the GNU Public License (GPLv3). Please include a citation to our [manuscript](http://www.atmos-chem-phys.net/16/4401/2016/) if used. The corresponding author can be contacted with any bug reports or questions.

Tables from the manuscript are provided in the "SMARTSpatterns/" subdirectory:

* "MCMgroups.csv": Table 1, substructures 1-30.
* "FTIRextra.csv": Table 1, substructures 31-54.
* "SIMPOLgroups.csv": Table 2, SIMPOL.1 groups. Use with `-e SIMPOLexportlist.csv` (see below) to get groups exported in the right order.
* "OScBonds.csv": Table 4, bonds for calculating oxidation state for each carbon atom.
* "MCMgroupsv33.csv": Update to MCMgroups.csv to accommodate conditions for aprl-carbontypes  (compatible with MCMv3.2) and new additions to MCMv3.3 isoprene mechanism. (Carbon specificity condition is not met for anhydrides.)

Supplementary scripts:

* "substructure\_molecular\_attributes.py": Extract molecular attributes that can be retrieved from a pybel Molecule object (e.g., molecular weight).

 Scripts and input files which reproduce the validation figures in the manuscript are also described below.

## Setup

### Operating system

You may want to add the directory of the program path so that the two python scripts can be accessed anywhere. In Linux/Mac (include in `~/.bash_profile` to make this available across various terminal sessions):

```sh
$ export PATH=$PATH:/path/to/aprl-ssp
```

In Windows (use `setx` instead of `set` to make available across various terminal sessions):

```dos
set PATH=%PATH%;Z:path\to\aprl-ssp
```
See [this page](http://www.computerhope.com/issues/ch000549.htm) for GUI instructions (there may be better sites). Note that these scripts have only been tested under Ubuntu Linux. Windows users will have to take note to enter path separators as `\` rather than `/`.

In examples below, it is assumed that the path has been added and the scripts are made executable (`chmod +x scriptname.py` in Linux). Otherwise, each command must be issued as `python /path/to/aprl-ssp/scriptname.py`. All examples are run assuming this has been configured, and the working directory is the "examples/" subdirectory (included in program directory).

### Python

All scripts require python (tested with 2.7). Anaconda distribution is recommended as it comes with the pandas library (assumed to be installed in further instructions). Python is installed with Mac OS X but a separate installation of the scientific python stack is recommended. Package installation can be carried through with the pip package manager:

```sh
$ pip install packagename
```

"substructure\_search.py" and "substructure\_generate\_fulltable.py" additionally require Open Babel (link to installation guide included below).

## Program scripts

### ----- spider\_query.py -----

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
* `-i`: value of `INPUTFILE` (optional). Name of file containing compounds; a csv file with a column called "compound". Optional if `-D`; otherwise required.
* `-t`: value of `TOKEN` (optional). Name of file containing ChemSpider token if not "~/.chemspidertoken". See next section.

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
$ spider_query.py -t /path/to/token.txt -p example -i compounds.csv
```

Optionally, the argument passed can be the value of the token itself.

#### Compound names

Where possible, superfluous spaces should be removed in the name of the compound. We have found inconsistencies in returned results. For instance:

* "benzo[b]fluoranthene" returns benzo[b]fluoranthene with SMILES string of "c1ccc2c(c1)cc-3c4c2cccc4-c5c3cccc5"
* "benzo[b] fluoranthene" returns (+)-Vincamine with SMILES string of "CC[C@@]12CCCN3[C@@H]1c4c(c5ccccc5n4[C@](C2)(C(=O)OC)O)CC3"


### ----- substructure\_search.py -----

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

* `-d`: When present, indicates that `GROUPFILE` (and `EXPORT`, if provided) exists in the subdirectory, `SMARTSpatterns/`, distributed with "substructure_search.py".

#### Examples

```
$ substructure_search.py -d -g SIMPOLgroups.csv -e SIMPOLexportlist.csv \
  -i apinenemech.csv -o apinenemech_SIMPOLgroups.csv
```

In this example, SIMPOLexportlist.csv provides a list of groups in the order that they should be exported.


### ----- substructure\_generate\_fulltable.py -----

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
* `-o`: value of `OUTPUTPREFIX`. Name of file prefix to be used for generated output files: {PREFIX}\_atomcounts.csv, {PREFIX}\_groupcounts.csv, {PREFIX}\_atomicmass.csv, {PREFIX}\_atomfulltable.csv. {PREFIX}_groupcounts.csv is similar to the output of substructure\_search.py but does not contain the full set of patterns in the `INPUTFILE` -- ony the matched ones.
* `-e`: value of `EXPORT` (optional). Name of file which contains list of substructures to include in `OUTPUTFILE`. Overrides "export" column if present in `GROUPFILE`. (*currently not implemented*)

Flags:

* `-d`: When present, indicates that `GROUPFILE` (and `EXPORT`, if provided) exists in the subdirectory, `SMARTSpatterns/`, distributed with "substructure_search.py".

#### Examples

```
$ substructure_generate_fulltable.py -d -g MCMgroups.csv \
  -i apinenemech.csv -o apinenemech_MCMgroups
```

### ----- generate\_carbontypes.py -----

Generate functional group and carbon type matrices (from aprl-carbontypes; adapted to python). 

#### Arguments

Main arguments:

* `-i`: value of `INPUTFILE`. File generated by substructure\_generate\_fulltable.py; csv format.
* `-o`: value of `OUTPUTPREFIX`. Output prefix.

#### Examples

```
$ generate_carbontypes.py -i apinene_MCMgroups_atomfulltable.csv -o apinene
```

### ----- validation_triple.py -----

Adopted from aprl-carbontypes. Supercedes validation\_atoms.R, except validation\_triple.py does not currently export table of unmatched or non-unique atoms [TODO].

#### Arguments

Main arguments:

* `-f`: value of `ATOMFULLTABLE`. File generated by substructure\_generate\_fulltable.py; csv format.
* `-a`: value of `ATOMCOMMON`. File generated by substructure\_seach\.py with SMARTSpatterns/common\_atoms.csv patterns; csv format.
* `-o`: value of `OUTPUTPREFIX`. Output prefix.

#### Examples

```
$ validation_triple.py -f apinene_MCMgroups_atomfulltable.csv -a apinene_common_atoms.csv -o apinene
```

### ----- simpol.py -----

Uses SIMPOL.1 (Pankow and Asher, doi:10.5194/acp-8-2773-2008, 2008) to estimate vapor pressure p^0 (atm) and \Delta H (kJ/mole) for a given temperature.

#### Arguments

Main arguments:

* `-i`: value of `INPUTFILE`. File generated by substructure\_search\.py with SMARTSpatterns/SIMPOLgroups.csv patterns; csv format.
* `-o`: value of `OUTFILE`. Name of outputfile.
* `-t`: value of `TEMP`. Temperature at which to calculate properties.


#### Examples

```
$ simpol.py -i apinene_SIMPOLgroups.csv -o apinene_props_298.csv -t 298.15
```

## Pattern files

Patterns specified in GROUPFILE can be derived from a combination of SMARTS patterns using set operations. For instance, `ester, all` is defined as `"[CX3,CX3H1](=O)[OX2H0][#6]"`. `nitroester` is defined as `"[#6][OX2H0][CX3,CX3H1](=O)[C;$(C[N+](=O)[O-]),$(CC[N+](=O)[O-]),$(CCC[N+](=O)[O-]),$(CCCC[N+](=O)[O-]),$(CCCCC[N+](=O)[O-])]"`. `ester` can be defined as `{ester, all}-{nitroester}`. When present, such custom patterns are computed after all the SMARTS patterns have been matched and counted. Additionally, functions can be provided by the user. In current implementation, functions would presumably use OpenBabel methods.

### Permissible entries:

* SMARTS pattern.
* Bracketed expression. E.g., `{ester, all}-{nitroester}`.
* Quoted expression. Substituted patterns in `'{}` are not evaluated before passing to function. E.g., `count_nitrophenol(molecule,'{phenol},'{ester})`.
* Keyword `eval` denotes another type of python expression, but can also be used to preface expressions containing brackets (redundant). E.g., `eval 1`, `eval count_aromatic_rings(molecule)`.
* Keyword `molecule` denotes the pybel Molecule class instance (useful for passing to user-defined functions). E.g., `eval count_aromatic_rings(molecule)`.

If only SMARTS patterns are used, the same pattern file can be provided to "substructure\_search.py" and "substructure\_search\_fulltable.py". When additional expressions are supplied, they must be changed such that arithmetic operations are used for matched entries in "substructure\_search.py" and set operations for "substructure\_search\_fulltable.py".

### User-supplied functions

User can define functions to be called for evaluation. These functions should be contained in a file called "userdef.py." These functions will generally accept as its first argument the pybel Molecule object upon which operations are to be performed. In the pattern file, the argument to the function should be provided as `molecule` (lowercase) as this is the object whose value will be evaluated.

## Validation scripts (for Ruggeri and Takahama, 2016)

Three scripts are provided for generating validation figures:

* "run\_validation.py": Calls "validation\_atoms.R" and "validation\_groups.R" and applies to specific CSV files in "validation/" directory (output figures are also placed there).
* "validation\_atoms.R": Tests completeness and specificity at the atom level (2 plots generated for each set of compounds).
* "validation\_groups.R": Tests that the number of groups are matched against manually enumerated values (1 plot generated for each set of compounds).

The main script, "run\_validation.py" script can be run in the "aprl-ssp/" directory. It will read in two lists of files and apply the comparison against reference and matched groups.

"validation/filelist\_atoms.csv" for validation of atom membership criteria (completeness and specificity) includes:

* apinenemech.csv
* apinenepropenemech.csv
* tmbmech.csv
* tmbpropenemech.csv

"validation/filelist\_groups.csv" for validation of group counts against manually enumerated values:

* compounds\_SIMPOLgroups.csv
* compounds\_FTIRextra.csv

Required sofware for validation:

* python (tested with 2.7)
* R (http://cran.r-project.org), also free. Required packages: `reshape2`, `dplyr`.

    ```
     > install.packages(c("reshape2","dplyr"))
	 ```


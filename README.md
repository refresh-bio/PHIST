# PHIST

[![C/C++ CI](https://github.com/refresh-bio/PHIST/workflows/C/C++%20CI/badge.svg)](https://github.com/refresh-bio/PHIST/actions)


**Phage-Host Interaction Search Tool**

A tool to predict prokaryotic hosts for phage (meta)genomic sequences. PHIST links viruses to hosts based on the number of *k*-mers shared between their sequences.

<p align="center"><img src="phist_logo.png" width="200"></p>

## Quick start
```bash
git clone --recurse-submodules https://github.com/refresh-bio/PHIST

cd PHIST
make
mkdir ./out

python3 phist.py ./example/virus ./example/host ./out/common_kmers.csv ./out/predictions.csv

```

## Installation

PHIST uses [Kmer-db](https://github.com/refresh-bio/kmer-db) as a submodule, therefore a recursive repository clone must be performed:
```
git clone --recurse-submodules https://github.com/refresh-bio/PHIST
```
Under Linux/OS X the package can be built by running MAKE in the project directory (G++ 5.3 tested):
```
cd PHIST
make
```
Under Windows one have to build Visual Studio 2015 solutions on *kmer-db* and *utils* subdirectories (use Release 64-bit configuration, as Python script depends on the default VS output directory structure).

## Usage

PHIST takes as input two directories containing FASTA files (gzipped or not) with genomic sequences of viruses and candidate hosts (see [example](./example/)).

```
python3 phist.py [options] <virus_dir> <host_dir> <out_table> <out_predictions>
```

Positional arguments:
  * `virus_dir`         input directory with virus FASTA files (gzipped or not),
  * `host_dir`          input directory with host FASTA files (gzipped or not),
  * `out_table`         output CSV file with common k-mers table (default: *common_kmers.csv*),
  * `out_predictions`   output CSV file with hosts predictions (default: *predictions.csv*).

Options:
* `-k --k <kmer-length>`   *k*-mer length (default: 25),
* `-t, --t <num-threads>`  number of threads (default: number of cores),
* `-h, --help`             show this help message and exit,
* `--version`              show tool's version number and exit.

## Output format

PHIST outputs two CSV files. One containing a table of common *k*-mers between phages and hosts, and second file with virus-host predictions.


### Common *k*-mers table

The [common_kmers.csv](./example/common_kmers.csv) file stores numbers of common *k*-mers between phages (in columns) and hosts (in rows) in a sparse form. Specifically, zeros are omitted while non-zero *k*-mer counts are represented as pairs (*column_number* : *value*) with 1-based column indexing. Thus, rows may have different number of elements, e.g.:

| 									| 								| 					| 				|		|			|	
| :---: 							| :---: 						| :---: 			| :---:			| :---:	|  :---:	| 
| kmer-length: *k* fraction: *f* 	| phages 					| *&phi;<sub>1</sub>*					| *&phi;<sub>2</sub>* | ... 	|  *&phi;<sub>n</sub>* |
| hosts 					| total-kmers 					| &#124;*&phi;<sub>1</sub>*&#124;		| &#124;*&phi;<sub>2</sub>*&#124; 	| ... 	|  &#124;*&phi;<sub>n</sub>*&#124; |
| *h<sub>1</sub>* 					| &#124;*h<sub>1</sub>*&#124;	| *i<sub>11</sub>* : &#124;*h<sub>1</sub> &cap; &phi;<sub>i<sub>11</sub></sub>*&#124;	| *i<sub>12</sub>* : &#124;*h<sub>1</sub> &cap; &phi;<sub>i<sub>12</sub></sub>*&#124; | ||
| *h<sub>2</sub>* 					| &#124;*h<sub>2</sub>*&#124;	| *i<sub>21</sub>* : &#124;*h<sub>2</sub> &cap; &phi;<sub>i<sub>21</sub></sub>*&#124;	| *i<sub>22</sub>* : &#124;*h<sub>2</sub> &cap; &phi;<sub>i<sub>22</sub></sub>*&#124; 	| *i<sub>23</sub>* : &#124;*h<sub>2</sub> &cap; &phi;<sub>i<sub>23</sub></sub>*&#124;  	| |   
| *h<sub>2</sub>* 					| &#124;*h<sub>2</sub>*&#124;	| ||||
| ... 								| ...							| ... ||||
| *h<sub>m</sub>* 					| &#124;*h<sub>m</sub>*&#124;	| *i<sub>m1</sub>* : &#124;*h<sub>m</sub> &cap; &phi;<sub>i<sub>m1</sub></sub>*&#124;	| |||

where:
* *k* - k-mer length,
* *&phi;<sub>1</sub>*, *&phi;<sub>2</sub>*,  ...,   *&phi;<sub>n</sub>* - phage names,
* *h<sub>1</sub>*, *h<sub>2</sub>*,  ...,   *h<sub>m</sub>* - host names,
* &#124;*a*&#124; - number of k-mers in sample *a*,
* &#124;*a &cap; b*&#124; - number of k-mers common for samples *a* and *b*.


### Host predictions

The [predictions.csv](./example/predictions.csv) file assigns each phage to its most likely host (i.e., the one having most *k*-mers in common). If there are multiple potential hosts with same number of common *k*-mers, all are reported. Each virus-host interaction is followed by *p*-value and adjusted *p*-value for multiple comparisons.

| 	phage								      | 		host						| 	common *k*-mers				| 	*p*-value			|	adj. *p*-value	|				
| :---: 							       | :---: 						| :---: 			           | :---:			     | :---:	 	       | 
|  *&phi;<sub>1</sub>*   | *host*( *&phi;<sub>1</sub>*) | &#124;*&phi;<sub>1</sub>* &cap; *host*(*&phi;<sub>1</sub>*)&#124; | ... | ... |
|  *&phi;<sub>2</sub>*   | *host*( *&phi;<sub>2</sub>*) | &#124;*&phi;<sub>2</sub>* &cap; *host*(*&phi;<sub>2</sub>*)&#124; | ... | ... |
|  *&phi;<sub>3</sub>*   | *host<sub>1</sub>*( *&phi;<sub>3</sub>*) | &#124;*&phi;<sub>3</sub>* &cap; *host<sub>1</sub>*(*&phi;<sub>3</sub>*)&#124; | ... | ... |
|  *&phi;<sub>3</sub>*   | *host<sub>2</sub>*( *&phi;<sub>3</sub>*) | &#124;*&phi;<sub>3</sub>* &cap; *host<sub>2</sub>*(*&phi;<sub>3</sub>*)&#124; | ... | ... |
| ... | ... | ... | ... | ... | ... |


## Citing
[Zielezinski, A., Deorowicz, S., Gudyś, A. (2021) PHIST: fast and accurate prediction of prokaryotic hosts from metagenomic viral sequences, bioRxiv,  https://doi.org/10.1101/2021.09.06.459169](https://www.biorxiv.org/content/10.1101/2021.09.06.459169v1)

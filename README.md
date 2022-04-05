# GISTIC2
This repository contains the Matlab source code for GISTIC (**G**enomic **I**dentification of **S**ignificant **T**argets **I**n **C**ancer), version 2.0.

## Note
This repository contains a submodule of matlab functions for processing copy number data named 'snputil.'
To ensure that the snputil subdirectory is populated with files, you should clone this repository using
the *--recursive* option:
```
git clone --recursive  https://github.com/broadinstitute/gistic2.git
```
## repository directory structure
> *docs* - HTML source for links and user documentation published in the GISTIC GitHub page.
> *refgenes* - source code for creating reference genome input files (offered on an 'as is' basis).
>> *Gencode.v22.170324* - download EMBL files and build a Gencode reference genome.
>> *hg38.UCSC.add_mir.160920* - download UCSC files and build an hg38 reference genome.
> *snputil* - git submodule of utility Matlab functions for analyzing copy number data.
>> *@SegArray* - defines the SegArray class, a data compression scheme for segmented copy number.
> *source* - Matlab source code for GISTIC.
> *support* - additional files needed to build a GISTIC tarball.
> *user_docs* - HTML source for standalone documentation in the current GISTIC tarball (v2.0.23).


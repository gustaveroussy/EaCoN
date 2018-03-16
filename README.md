# EaCoN

Easy Copy Number !
A user-friendly R package to perform copy-number analysis from Affymetrix microarrays (SNP6, OncoScan, CytoScan 750K, CytoScan HD) or Whole Exome Sequencing data.	



# Authors: 
 - Bastien Job (bastien.job@inserm.fr, bastien.job@gustaveroussy.fr) : developer and maintainer
 - Thibault Dayris (thibault.dayris@gustaveroussy.fr) : beta-tester


## Core requirements

- R
- Python3-pip:
    - typing
    - snakemake
    - json
    - yaml

Please, note that each tool must be installed aside and usable.
STRonGR do not intend to test tools' -- or their dependencies' -- installation.

## Usage

You can call STRonGR with *--help* or *-h* argument to have a list of all available pipelines. If you call STRonGR with a pipeline name as a subcommand, then a dedicated help prints out.

Short hints about the command line are available without any argument.

```{bash}
$ python3 STRonGR.py Salmon_0_8_2 --help
usage: STRonGR.py Salmon_0_8_2 [-h] [-f FASTA] [-i INDEX] [--library LIBRARY]
                               [--bootstrap BOOTSTRAP]
                               [--kmer_length KMER_LENGTH] [--gencode]
                               [--perfect_hash] [--bias_correct] [-o PATH]
                               [-l PATH] [-m MAIL@PROVIEDER.COM] [-t THREADS]
                               [--json] [--dry_run] [--force_all] [--cluster]
                               [--dag] [--quiet] [--keeptmp]
                               design

    Perform lightweight quantification with Salmon (v0.8.2)!

    You may find more information about:
    Libraries:
    http://salmon.readthedocs.io/en/latest/library_type.html

    Salmon usage:
    https://salmon.readthedocs.io/en/latest/salmon.html#using-salmon

    Salmon output files:
    https://salmon.readthedocs.io/en/latest/salmon.html#using-salmon

...
```
## Join the STRonGR project !

Add things here:

- **pypelines** contains all your python classes intended to run your pipeline
- **snakefile** contains all your snakemake rules designed to run your pipeline
- **scripts** contains non pythonic scripts required for your pipelines.

Better not touch to:

- **core repository** contains deep STRonGR functions. You should not have to modify the content of this repository.

Run with:

- **STRonGR.py** automatically detects available pipelines (pair of pypeline and snakefile) and wraps it all into snakemake

Enjoy!

## Share about the project on trello:

https://trello.com/b/H1SLtjeQ/strongr

## Guidelines

Within Python: follow the [pep8](https://www.python.org/dev/peps/pep-0008/)
As potentially many people will look / use your scripts, let's code the same way!

Snakemake:

- **Use wildcards**: This is the main force of Snakemake. Thanks to them, your pipeline will be usable by anyone, with any accepted file.
- **Use local external file**: By local file, we mean files present in your working directory. External files will be *symlinked* by STRonGR for you.
- **Unique rule names**: When you write rules, always try ton make their name unique. Don't forget that snakefiles are to be included in each others! A goo rule name is: "Tool_Version_purpose".
- **No space, no dots**: Dots (.) are intended to separate extensions, underscores (\_) are here to separate fields. Use the snake_case to name things.
- **Do one thing, and do it well**: A rule should have one purpose and only one. Multiple operation are usually better handled in separate rule.

# MetaCarvel - Scaffolder for metagenomes

MetaCarvel is an updated version of previous metagenome scaffolder Bambus 2. To run MetaCarvel, you will need [Samtools](http://samtools.sourceforge.net), [Bedtools](http://bedtools.readthedocs.io/en/latest/), [Networkx](https://networkx.github.io/)(Version < 1.11), and [OGDF](http://amber-v7.cs.tu-dortmund.de/lib/exe/fetch.php/tech:ogdf-snapshot-2015-05-30.zip).

You can install Networkx as described [here](https://pypi.org/project/networkx/1.10/).
Briefly, you need to run following:
```
pip install networkx==1.10
```
If this doesn't work, you can download the source from [here](https://files.pythonhosted.org/packages/cb/fc/9b1c805b9abe249b9ce786d2ac9e6808d7776237d195365d50188a38dc30/networkx-1.10.tar.gz) and install from source.

The detailed documentation and tutorial to install and run MetaCarvel can be found [here](http://bambus3.readthedocs.io/en/latest/).

To install the software, please execute following commands:
```
cd MetaCarvel
make
```
Before running `make` command, please make sure you have modified following in the `makefile` to suit your OGDF installation location:

```
OGDF_INCL = -I /cbcb/sw/RedHat-7-x86_64/users/jayg/local/OGDF/2015.05/include/
OGDF_LINK = -L /cbcb/sw/RedHat-7-x86_64/users/jayg/local/OGDF/2015.05/_release/
```

To run MetaCarvel, run the following;

```
python run_new.py -h
usage: run_new.py [-h] -a ASSEMBLY -m MAPPING -d DIR [-f FORCE] [-r REPEATS]
              [-k KEEP] [-l LENGTH] [-b BSIZE]

MetaCarvel: A scaffolding tool for metagenomic assemblies

optional arguments:
  -h, --help            show this help message and exit
  -a ASSEMBLY, --assembly ASSEMBLY
                        assembled contigs
  -m MAPPING, --mapping MAPPING
                        mapping of read to contigs in bam format
  -d DIR, --dir DIR     output directory for results
  -f FORCE, --force FORCE
                        force re-run of pipeline, will remove any existing
                        output
  -r REPEATS, --repeats REPEATS
                        To turn repeat detection on
  -k KEEP, --keep KEEP  Set this to kepp temporary files in output directory
  -l LENGTH, --length LENGTH
                        Minimum length of contigs to consider for scaffolding
  -b BSIZE, --bsize BSIZE
                        Minium mate pair support between contigs to consider
                        for scaffolding
```

This will generate a bunch of files in the output directory. If you are interested in output of each step of the scaffolding process, these files can 
be useful. The final output files are scaffolds.fasta - which contains sequences of scaffolds  and scaffolds.agp is an agp style information for assignment of contigs to scaffolds. 


NOTE: This tool is still under active development and may produce errors while running. Please report these as github issues so that we can fix them as we develop the software. For any questions, please email jayg@cs.umd.edu. 

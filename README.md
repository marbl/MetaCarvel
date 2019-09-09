# MetaCarvel - Scaffolder for metagenomes

MetaCarvel is an updated version of previous metagenome scaffolder Bambus 2. To run MetaCarvel, you will need [Python 2.7.x](https://www.python.org/downloads/), [Samtools](http://samtools.sourceforge.net), [Bedtools](http://bedtools.readthedocs.io/en/latest/), [Networkx](https://networkx.github.io/)(Version < 1.11), [NumPy](http://www.numpy.org/),and [OGDF](http://amber-v7.cs.tu-dortmund.de/lib/exe/fetch.php/tech:ogdf-snapshot-2015-05-30.zip).

You can install Networkx as described [here](https://pypi.org/project/networkx/1.10/).
Briefly, you need to run following:
```
pip install networkx==1.10
```
If this doesn't work, you can download the source from [here](https://files.pythonhosted.org/packages/cb/fc/9b1c805b9abe249b9ce786d2ac9e6808d7776237d195365d50188a38dc30/networkx-1.10.tar.gz) and install from source.

## The detailed documentation and tutorial to install and run MetaCarvel can be found on [Wiki](https://github.com/marbl/MetaCarvel/wiki).


To run MetaCarvel, run the following;

```
python run.py -h
usage: run.py [-h] -a ASSEMBLY -m MAPPING -d DIR [-r REPEATS] [-k KEEP]
              [-l LENGTH] [-b BSIZE] [-v VISUALIZATION]

MetaCarvel: A scaffolding tool for metagenomic assemblies

optional arguments:
  -h, --help            show this help message and exit
  -a ASSEMBLY, --assembly ASSEMBLY
                        assembled contigs
  -m MAPPING, --mapping MAPPING
                        mapping of read to contigs in bam format
  -d DIR, --dir DIR     output directory for results
  -r REPEATS, --repeats REPEATS
                        To turn repeat detection on
  -k KEEP, --keep KEEP  Set this to keep temporary files in output directory
  -l LENGTH, --length LENGTH
                        Minimum length of contigs to consider for scaffolding
                        in base pairs (bp)
  -b BSIZE, --bsize BSIZE
                        Minimum mate pair support between contigs to consider
                        for scaffolding
  -v VISUALIZATION, --visualization VISUALIZATION
                        To generate .db file for AsmViz visualization program
```

This will generate a bunch of files in the output directory. If you are interested in output of each step of the scaffolding process, these files can 
be useful. The final output files are scaffolds.fasta - which contains sequences of scaffolds  and scaffolds.agp is an agp style information for assignment of contigs to scaffolds. 

Please cite MetaCarvel as follows: Ghurye, J., Treangen, T., Fedarko, M., Hervey, W. J., & Pop, M. (2019). MetaCarvel: linking assembly graph motifs to biological variants. Genome biology, 20(1), 1-14.

NOTE: This tool is still under active development and may produce errors while running. Please report these as github issues so that we can fix them as we develop the software. For any questions, please email jayg@cs.umd.edu. 

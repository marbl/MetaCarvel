# Bambus 3 - Scaffolder for metagenomes

Bambus3 is an updated version of previous metagenome scaffolder Bambus 2. To run Bambus 3, you will need [Samtools](http://samtools.sourceforge.net), [Bedtools](http://bedtools.readthedocs.io/en/latest/), and [Networkx](https://networkx.github.io/).

To install the software, please execute following commands:
```
cd bambus3
make
```

To run bambus3, run the following;

```
python run.py -h
usage: run.py [-h] -a ASSEMBLY -m MAPPING -d DIR [-f FORCE]

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
```

This will generate a bunch of files in the output directory. If you are interested in output of each step of the scaffolding process, these files can 
be useful. The final output files are scaffolds.fasta - which contains sequences of scaffolds  and scaffolds.agp is an agp style information for assignment of contigs to scaffolds. 


NOTE: This tool is still under active development and may produce errors while running. Please report these as github issues so that we can fix them as we develop the software. For any questions, please email jayg@cs.umd.edu. 

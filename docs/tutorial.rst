MetaCarvel Pipeline
===================

MetaCarvel is a scaffolding program designed specifically to address the issues concerning the metagenomic data. 
The scaffolding algorithm consists of several steps. MetaCarvel requires contig assembly fasta file, mapping of reads to contigs in BAM file and output directory as an input. It first analyzes the read to contig mapping and estimates the library insert size. Using these mappings, it generates a scaffold graph where nodes are contigs and edges are mate pairs mapped to two different contigs. This is done in :code:`libcorrect` program. Multiple links between two contigs are collapsed into a single link using :code:`bundler` program. After this, repeats in the graph are removed. High betweenness centrality nodes are identified using :code:`centrality` program and removed from the graph using :code:`repeat_filter.py` program. Contigs remaining after filtering repetitive contigs are assigned a unique orientation using :code:`orientcontigs` program. The bubbles caused due to polymorphic regions in the metagenome are identified using :code:`spqr` program. Once these bubbles are identified, they are collapsed and contig layout is performed in :code:`layout.py` program.

Installing the software
===================
MetaCarvel is open source and available on github. You need to run following steps to get it installed

.. code::

	git clone https://github.com/marbl/bambus3
	cd bambus3
	make

This would generate the executables required to run bambus3.

Preparing the data
===================
In order to run MetaCarvel, you will need two files. First one is the contig assembly of the raw reads using your favorite metagenome assembly program. 
The second one is the alignment of raw paired end reads to assembled contigs using your favorite alignment program. We strongly recommend to map paired end reads as single end reads since aligner enforces the reads to map in the distances given by the library size. If you are using Bowtie2 for read mapping, we recommend aligning reads this way:

.. code::

	bowtie2 -x idx -U forward_reads.fq | samtools view -bS - > alignment_1.bam
	bowtie2 -x idx -U reverse_reads.fq | samtools view -bS - > alignment_2.bam
	samtools alignment_total.bam alignment_1.bam alignment_2.bam
	samtools sort -n alignment_total.bam -o alignment_total_sorted.bam
We need this file sorted by the read name and not the coordinates. You will also need to provide the output directory where you want your scaffolds to be written. 


Running MetaCarvel
===================
Before you run MetaCarvel, please make sure you have following dependencies: `samtools <http://bowtie-bio.sourceforge.net/manual.shtml>`_, `bedtools <http://bedtools.readthedocs.io/en/latest/>`_ and `NetworkX <https://networkx.github.io/>`_. After you have these dependencies, running MetaCarvel is pretty straightforward. You would need to execute :code:`run.py` file. 

.. code::

	python run.py -h
	usage: run.py [-h] -a ASSEMBLY -m MAPPING -d DIR [-f FORCE] [-r REPEATS]
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

The mandatory inputs for this are the FASTA file for contig assembly, the alignment file and output directory. There are other options as well. 
:code:`-r` option lets you toggle the repeat detection process in the scaffolding pipeline. By default it is turned off, but you can provide :code:`-r true` parameter to include repeat detection. :code: `-l` option lets you set the minimum length of the contigs to consider for scaffolding. By default this is set to :code: 500. :code:`-b` option lets you set the minimum mate pair support between two contigs to consider for scaffolding. By default this is set to :code:`3`. Each step in the pipeline generates temporary files as intermediate outputs. If you want to retain these files, you need to set :code:`-k` option to :code:`true`. By default, these files will be removed. 

Interpreting the output
===================
Since MetaCarvel consists of multiple steps, it generates some intermediate files. These files can be useful if you want to dig more into the analysis. 
The :code:`libcorrect` generates two output files: :code:`contig_links`, having a raw links between the contigs based on mate pairs and :code:`contig_coverage`, having a coverage information for each contig. The :code:`bundler` program uses :code:`contig_links` file, bundles the links and outputs them in the :code:`bundled_links` file in a tsv format and in the :code:`bundled_graph.gml` in the gml format. If repeat detection is chosen then that would generate a file with name :code:`bundled_links_filtered` which will have the bundled links corresponding to all the non-repetitive contigs. The :code:`orientcontigs` program takes the bundled links as an input and outputs a graph with only one orientation corresponding to each contig. This graph is written in two formats: :code:`oriented_links` in tsv format and :code:`oriented.gml` in a gml format. The :code:`spqr` program takes orientated graph as an input and produces a potential bubbles as an output in a file called :code:`seppairs`. The :code:`layout.py` produces three main files. One is :code:`scaffolds.fa`, represnting the scaffold sequences, :code:`scaffolds.agp`, representing scaffolds in an `AGP <https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/>`_ format and :code:`scaffolds.gfa` representing scaffold graph in `GFA <https://github.com/GFA-spec/GFA-spec>`_ format. If you want to visualize the graphs used to generate scaffolds, you can either use either `oriented.gml` or `scaffolds.gfa` file and load it into any standard graph visualization software. We recommend using `MetagenomeScope `https://marbl.github.io/MetagenomeScope/>`_ as you can visualize variants and scaffold paths simply by loading the agp output into the viewer. 

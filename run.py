#!/usr/bin/env python3

import os
import argparse
import sys
import time
import subprocess


def cmd_exists(cmd):
    return subprocess.call("type " + cmd, shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

def print_msg(msg):
    ts = time.strftime("%c")
    print(f"{ts}: {msg}")

def print_err(msg):
    ts = time.strftime("%c")
    print(f"{ts}: {msg}", file=sys.stderr)


def main():
    cwd = os.path.dirname(os.path.abspath(__file__))
    print("CWD:", cwd)

    parser = argparse.ArgumentParser(description="MetaCarvel: A scaffolding tool for metagenomic assemblies")
    parser.add_argument("-a", "--assembly", help="assembled contigs", required=True)
    parser.add_argument("-m", "--mapping", help="mapping of read to contigs in bam format", required=True)
    parser.add_argument("-d", "--dir", help="output directory for results", default='out', required=True)
    parser.add_argument("-r", '--repeats',help="To turn repeat detection on", default="true")
    parser.add_argument("-k", "--keep", help="Set this to keep temporary files in output directory", default=False)
    parser.add_argument("-l", "--length", help="Minimum length of contigs to consider for scaffolding in base pairs (bp)", default=500)
    parser.add_argument("-b", "--bsize", help="Minimum mate pair support between contigs to consider for scaffolding", default=3)
    parser.add_argument("-v", '--visualization', help="Generate a .db file for the MetagenomeScope visualization tool", default=False)

    args = parser.parse_args()
    try:
        import networkx
    except ImportError:
        print_err('Looks like you do not have networkx. Please rerun with networkx module installed.')
        sys.exit(1)
    
    version = networkx.__version__
    print_msg(f"Networkx {version} found")
    
    if not cmd_exists('samtools'):
        print_err('Samtools does not exist in PATH. Terminating...')
        sys.exit(1)

    if not cmd_exists('bamToBed'):
        print_err('Bedtools does not exist in PATH. Terminating...')
        sys.exit(1)

    if not os.path.exists(args.dir):
        os.makedirs(args.dir)
    
    print_msg('Starting scaffolding...')
    if not os.path.exists(args.dir + '/alignment.bed'):
        print_msg("Converting bam file to bed file")
        try:
            subprocess.check_output('bamToBed -i ' + args.mapping + " > " + args.dir+'/alignment.bed', shell=True)
            print_msg('Finished conversion')
        except subprocess.CalledProcessError as err:
            os.system(f"rm {args.dir}/alignment.bed")
            print_err('Failed in coverting bam file to bed format, terminating scaffolding...' + str(err.output))
            sys.exit(1)
    try:
        subprocess.check_output(f"samtools faidx {args.assembly}", shell=True)
    except subprocess.CalledProcessError as err:
        print_err(str(err.output))
        sys.exit(1)

    os.system(f"cut -f 1,2 {args.assembly}.fai > {args.dir}/contig_length")
    print_msg('Finished conversion')

    print_msg('Started generating links between contigs')
    if not os.path.exists(f"{args.dir}/contig_links"):
        try:
            subprocess.check_output(f"{cwd}/libcorrect -a {args.dir}/alignment.bed \
                                    -d {args.dir}/contig_length -o {args.dir}/contig_links \
                                    -x {args.dir}/contig_coverage -c {str(args.length)}", shell=True)
            print_msg('Finished generating links between contigs')
        except subprocess.CalledProcessError as err:
            os.system(f"rm {args.dir}/contig_links")
            print_err('Failed in generate links from bed file, terminating scaffolding...' + str(err.output))
            sys.exit(1)

    print_msg('Started bulding links between contigs')
    if not os.path.exists(f"{args.dir}/bundled_links"):
        try:
            subprocess.check_output(f"{cwd}/bundler -l {args.dir}/contig_links \
                                    -o {args.dir}/bundled_links + -b {args.dir}/bundled_graph.gml \
                                    -c {str(args.bsize)}", shell=True)
            print_msg('Finished bundling of links between contigs')
        except subprocess.CalledProcessError as err:
            os.system(f"rm {args.dir}/bundled_links")
            os.system(f"rm {args.dir}/bundled_graph.gml")
            print_err('Failed to bundle links, terminating scaffolding...' + str(err.output))
            sys.exit(1)

    if args.repeats == "true":
        print_msg('Started finding and removing repeats')

        try:
            subprocess.check_output(f"{cwd}/orientcontigs -l {args.dir}/bundled_links \
                                    -c {args.dir}/contig_length --bsize -o {args.dir}/oriented.gml \
                                    -p {args.dir}/oriented_links -i {args.dir}/invalidated_counts", shell=True)
        except subprocess.CalledProcessError as err:
            print_err('Failed to find repeats, terminating scaffolding...' + str(err.output))

        try:
            subprocess.check_output(f"python {cwd}/centrality.py -g {args.dir}/bundled_links \
                                    -l {args.dir}/contig_length -o {args.dir}/high_centrality.txt", shell=True)
        except subprocess.CalledProcessError as err:
            print_err('Failed to find repeats, terminating scaffolding...' + str(err.output))
            sys.exit(1)

        try:
            subprocess.check_output(f"python {cwd}/repeat_filter.py {args.dir}/contig_coverage \
                                    {args.dir}/bundled_links {args.dir}/invalidated_counts \
                                    {args.dir}/high_centrality.txt {args.dir}/contig_length {args.dir}/repeats \
                                    > {args.dir}/bundled_links_filtered", shell=True)
        except subprocess.CalledProcessError as err:
            print_err('Failed to find repeats, terminating scaffolding...' + str(err.output))
            sys.exit(1)

        print_msg('Finished repeat finding and removal')
    else:
        os.system(f"mv {args.dir}/bundled_links {args.dir}/bundled_links_filtered")

    print_msg('Started orienting the contigs')
    try:
        subprocess.check_output(f"{cwd}/orientcontigs -l {args.dir}/bundled_links_filtered -c {args.dir}/contig_length \
                                --bsize -o {args.dir}/oriented.gml -p {args.dir}/oriented_links \
                                -i {args.dir}/invalidated_counts", shell=True)
        print_msg('Finished orienting the contigs')
    except subprocess.CalledProcessError:
        print_err('Failed to orient the contigs, terminating scaffolding...')

    print_msg('Started finding separation pairs')
    try:
        subprocess.check_output(f"{cwd}/spqr -l {args.dir}/oriented_links -o {args.dir}/seppairs", shell=True)
        print_err('Finished finding spearation pairs')
    except subprocess.CalledProcessError as err:
        print_err('Failed to decompose graph, terminating scaffolding...' + str(err.output))
        sys.exit(1)

    print_msg('Finding the layout of contigs')
    if not os.path.exists(f"{args.dir}/scaffolds.fasta"):
        try:
            subprocess.check_output(f"python {cwd}/layout.py -a {args.assembly} -b {args.dir}/bubbles.txt \
                                    -g {args.dir}/oriented.gml -s {args.dir}/seppairs -o {args.dir}/scaffolds.fa \
                                    -f {args.dir}/scaffolds.agp -e {args.dir}/scaffold_graph.gfa", shell=True)
            print_msg('Final scaffolds written, Done!')
        except subprocess.CalledProcessError as err:
            print_err('Failed to generate scaffold sequences, terminating scaffolding...' + str(err.output))

    if args.visualization == "true":
        graphpath = os.path.abspath(f"{args.dir}/oriented.gml")
        bubblepath = os.path.abspath(f"{args.dir}/bubbles.txt")
        # Output the MetagenomeScope .db file directly to args.dir. The only file
        # created by collate.py here is the mgsc.db file.
        os.system(f"python {cwd}/MetagenomeScope/graph_collator/collate.py -i {graphpath} -w -ub \
                  {bubblepath} -ubl -d {args.dir} -o mgsc")

    if not args.keep == "true":
        for fname in ['contig_length', 'contig_links', 'contig_coverage', 'bundled_links',
                      'bundled_links_filtered', 'bundled_graph.gml', 'invalidated_counts',
                      'repeats', 'oriented_links', 'oriented.gml', 'seppairs', 'alignment.bed']:
            path = os.path.join(args.dir, fname)
            if os.path.exists(path):
                os.system(f"rm {path}")

if __name__ == '__main__':
    main()

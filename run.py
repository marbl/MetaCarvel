from pathlib import Path
import argparse
import csv
import os
import shutil
import subprocess
import sys
import time


def main():
    cwd=Path(__file__).resolve().parent.absolute()

    parser = argparse.ArgumentParser(description="MetaCarvel: A scaffolding tool for metagenomic assemblies")
    parser.add_argument("-a","--assembly",help="assembled contigs",required=True)
    parser.add_argument("-m","--mapping", help="mapping of read to contigs in bam format",required=True)
    parser.add_argument("-d","--dir",help="output directory for results",default='out',required=True)
    parser.add_argument("-r",'--repeats',help="To turn repeat detection on",default="true")
    parser.add_argument("-k","--keep", help="Set this to keep temporary files in output directory",default=False)
    parser.add_argument("-l","--length",help="Minimum length of contigs to consider for scaffolding in base pairs (bp)",default=500)
    parser.add_argument("-b","--bsize",help="Minimum mate pair support between contigs to consider for scaffolding",default=3)
    parser.add_argument("-v",'--visualization',help="Generate a .db file for the MetagenomeScope visualization tool",default=False)

    args = parser.parse_args()
    try:
      import networkx
    except ImportError:
      raise ImportError('Looks like you do not have networkx. Please rerun with networkx module installed.')
      sys.exit(1)
    version = networkx.__version__
    version_id = int(version.split('.')[1])
    first = int(version.split('.')[0])
    print('Networkx ' + version + ' found', file=sys.stderr)
   # if first != 1:
   #    print(time.strftime("%c")+': Networkx should be 1.11 or earlier.. Terminating...\n', file=sys.stderr)
   #    sys.exit(1)
    if not shutil.which('samtools'):
      print(time.strftime("%c")+': Samtools does not exist in PATH. Terminating....\n', file=sys.stderr)
      sys.exit(1)

    if not shutil.which('bamToBed'):
      print(time.strftime("%c")+': Bedtools does not exist in PATH. Terminating....\n', file=sys.stderr)
      sys.exit(1)

    if not os.path.exists(args.dir):
        os.makedirs(args.dir)
    print(time.strftime("%c")+':Starting scaffolding..', file=sys.stderr)


    alignment_bed = Path(args.dir) / "alignment.bed"
    if not alignment_bed.exists():
        print("converting bam file to bed file", file=sys.stderr)
        try:
          with open(alignment_bed, "w") as alignment_bed_f:
            subprocess.run([
              "bamToBed",
              "-i",
              args.mapping,
            ], stdout=alignment_bed_f, check=True)
          print('finished conversion', file=sys.stderr)
        except subprocess.CalledProcessError as err:
          alignment_bed.unlink()
          print(time.strftime("%c")+': Failed in coverting bam file to bed format, terminating scaffolding....\n' + str(err.output), file=sys.stderr)
          sys.exit(1)
    try:
      subprocess.run([
        "samtools",
        "faidx",
        args.assembly,
      ], check=True)
    except subprocess.CalledProcessError as err:
      print(str(err.output), file=sys.stderr)
      sys.exit(1)

    with open(args.assembly+'.fai', "r", newline='') as f_in, open(args.dir+'/contig_length', "w", newline='') as f_out:
      reader = csv.DictReader(f_in, fieldnames=["name", "length", "offset", "linebases", "linewidth"], delimiter="\t")
      writer = csv.DictWriter(f_out, fieldnames=["name", "length"], delimiter="\t")
      for row in reader:
          writer.writerow({"name": row["name"], "length": row["length"]})

    print(time.strftime("%c")+':Finished conversion', file=sys.stderr)

    final_assembly = args.assembly
    final_mapping = args.mapping

    print(time.strftime("%c") + ':Started generating links between contigs', file=sys.stderr)
    contig_links = Path(args.dir) / "contig_links"
    contig_length = Path(args.dir) / "contig_length"
    contig_coverage = Path(args.dir) / "contig_coverage"
    if os.path.exists(args.dir+'/contig_links') == False:
        #print './libcorrect -l' + args.lib + ' -a' + args.dir+'/alignment.bed -d ' +args.dir+'/contig_length -o '+ args.dir+'/contig_links'
        try:
          subprocess.run([
            str(cwd / "libcorrect"),
            "-a",
            args.dir+"/alignment.bed",
            "-d",
            str(contig_length),
            "-o",
            str(contig_links),
            "-x",
            str(contig_coverage),
            "-c",
            str(args.length),
          ], check=True)
          print(time.strftime("%c") +':Finished generating links between contigs', file=sys.stderr)
        except subprocess.CalledProcessError as err:
          contig_links.unlink()
          print(time.strftime("%c")+': Failed in generate links from bed file, terminating scaffolding....\n' + str(err.output), file=sys.stderr)
          sys.exit(1)

    print(time.strftime("%c")+':Started bulding of links between contigs', file=sys.stderr)
    bundled_links = Path(args.dir) / "bundled_links"
    bundled_graph = Path(args.dir) / "bundled_graph.gml"
    if os.path.exists(args.dir+'/bundled_links') == False:
        try:
          subprocess.run([
            str(cwd / "bundler"),
            "-l",
            str(contig_links),
            "-o",
            str(bundled_links),
            "-b",
            str(bundled_graph),
            "-c",
            str(args.bsize),
          ], check=True)
          print(time.strftime("%c")+':Finished bundling of links between contigs', file=sys.stderr)
        except subprocess.CalledProcessError as err:
          bundled_links.unlink()
          bundled_graph.unlink()
          print(time.strftime("%c")+': Failed to bundle links, terminating scaffolding....\n' + str(err.output), file=sys.stderr)
          sys.exit(1)

    invalidated_counts = Path(args.dir) / "invalidated_counts"
    high_centrality = Path(args.dir) / "high_centrality.txt"
    bundled_links_filtered = Path(args.dir) / "bundled_links_filtered"
    oriented_gml = Path(args.dir) / "oriented.gml"
    oriented_links = Path(args.dir) / "oriented_links"
    repeats = Path(args.dir) / "repeats"
    if args.repeats == "true":
        print(time.strftime("%c")+':Started finding and removing repeats', file=sys.stderr)
        try:
          subprocess.run([
            str(cwd / "orientcontigs"),
            "-l",
            str(bundled_links),
            "-c",
            str(contig_length),
            "--bsize",
            "-o",
            str(oriented_gml),
            "-p",
            str(oriented_links),
            "-i",
            str(invalidated_counts),
          ], check=True)
        except subprocess.CalledProcessError as err:
            print(time.strftime("%c") + ': Failed to find repeats, terminating scaffolding...\n' + str(err.output), file=sys.stderr)

        try:
          subprocess.run([
            "python",
            str(cwd / "centrality.py"),
            "-g",
            str(bundled_links),
            "-l",
            str(contig_length),
            "-o",
            str(high_centrality),
          ], check=True)
        except subprocess.CalledProcessError as err:
          print(time.strftime("%c")+': Failed to find repeats, terminating scaffolding....\n' + str(err.output), file=sys.stderr)
          sys.exit(1)

        try:
          with open(bundled_links_filtered, "w") as bundled_links_filtered_f:
            subprocess.run([
              "python",
              str(cwd / "repeat_filter.py"),
              str(contig_coverage),
              str(bundled_links),
              str(invalidated_counts),
              str(high_centrality),
              str(contig_length),
              str(repeats),
            ], stdout=bundled_links_filtered_f, check=True)
        except subprocess.CalledProcessError as err:
            print(time.strftime("%c")+': Failed to find repeats, terminating scaffolding....\n' + str(err.output), file=sys.stderr)
            sys.exit(1)
        print(time.strftime("%c")+':Finished repeat finding and removal', file=sys.stderr)
    else:
      bundled_links.rename(bundled_links_filtered)
    print(time.strftime("%c")+':Started orienting the contigs', file=sys.stderr)
    try:
      subprocess.run([
        str(cwd / "orientcontigs"),
        "-l",
        str(bundled_links_filtered),
        "-c",
        str(contig_length),
        "--bsize",
        "-o",
        str(oriented_gml),
        "-p",
        str(oriented_links),
        "-i",
        str(invalidated_counts),
      ], check=True)
      print(time.strftime("%c")+':Finished orienting the contigs', file=sys.stderr)
    except subprocess.CalledProcessError:
      print(time.strftime("%c")+': Failed to Orient contigs, terminating scaffolding....', file=sys.stderr)

    print(time.strftime("%c")+':Started finding separation pairs', file=sys.stderr)
    seppairs = Path(args.dir) / "seppairs"
    try:
      subprocess.run([
        str(cwd / "spqr"),
        "-l",
        str(oriented_links),
        "-o",
        str(seppairs),
      ], check=True)
      print(time.strftime("%c")+':Finished finding spearation pairs', file=sys.stderr)
    except subprocess.CalledProcessError as err:
      print(time.strftime("%c")+': Failed to decompose graph, terminating scaffolding....\n' + str(err.output), file=sys.stderr)
      sys.exit(1)

    print(time.strftime("%c")+':Finding the layout of contigs', file=sys.stderr)
    bubbles = Path(args.dir) / "bubbles.txt"
    if os.path.exists(args.dir+'/scaffolds.fasta') == False:
      try:
        subprocess.run([
          "python",
          str(cwd / "layout.py"),
          "-a",
          str(args.assembly),
          "-b",
          str(bubbles),
          "-g",
          str(oriented_gml),
          "-s",
          str(seppairs),
          "-o",
          args.dir+"/scaffolds.fa",
          "-f",
          args.dir+"/scaffolds.agp",
          "-e",
          args.dir+"/scaffold_graph.gfa",
        ], check=True)
        print(time.strftime("%c")+':Final scaffolds written, Done!', file=sys.stderr)
      except subprocess.CalledProcessError as err:
        print(time.strftime("%c")+': Failed to generate scaffold sequences, terminating scaffolding....\n' + str(err.output), file=sys.stderr)

    if args.visualization == "true":
      # Output the MetagenomeScope .db file directly to args.dir. The only file
      # created by collate.py here is the mgsc.db file.
      subprocess.run([
        "python",
        str(cwd / "MetagenomeScope/graph_collator/collate.py"),
        "-i",
        str(oriented_gml),
        "-w",
        "-ub",
        str(bubbles),
        "-ubl",
        "-d",
        args.dir,
        "-o",
        "mgsc",
      ], check=True)

    if not args.keep == "true":
      for gen_f in [
        contig_length,
        contig_links,
        contig_coverage,
        bundled_links,
        bundled_links_filtered,
        bundled_graph,
        invalidated_counts,
        repeats,
        oriented_links,
        oriented_gml,
        seppairs,
        alignment_bed,
      ]:
        gen_f.unlink(missing_ok=True)

if __name__ == '__main__':
    main()

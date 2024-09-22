from pathlib import Path
import argparse
import csv
import shutil
import subprocess
import sys
import time


def main():
    cwd = Path(__file__).resolve().parent.absolute()

    parser = argparse.ArgumentParser(
        description="MetaCarvel: A scaffolding tool for metagenomic assemblies"
    )
    parser.add_argument(
        "-a", "--assembly", type=Path, help="assembled contigs", required=True
    )
    parser.add_argument(
        "-m",
        "--mapping",
        type=Path,
        help="mapping of read to contigs in bam format",
        required=True,
    )
    parser.add_argument(
        "-d",
        "--dir",
        type=Path,
        help="output directory for results",
        default="out",
        required=True,
    )
    parser.add_argument(
        "-r",
        "--repeats",
        action=argparse.BooleanOptionalAction,
        help="To turn repeat detection on",
        default=True,
    )
    parser.add_argument(
        "-k",
        "--keep",
        action=argparse.BooleanOptionalAction,
        help="Set this to keep temporary files in output directory",
        default=False,
    )
    parser.add_argument(
        "-l",
        "--length",
        type=int,
        help="Minimum length of contigs to consider for scaffolding in base pairs (bp)",
        default=500,
    )
    parser.add_argument(
        "-b",
        "--bsize",
        type=int,
        help="Minimum mate pair support between contigs to consider for scaffolding",
        default=3,
    )
    parser.add_argument(
        "-v",
        "--visualization",
        action=argparse.BooleanOptionalAction,
        help="Generate a .db file for the MetagenomeScope visualization tool",
        default=False,
    )

    args = parser.parse_args()
    try:
        import networkx
    except ImportError:
        raise ImportError(
            "Looks like you do not have networkx. Please rerun with networkx module installed."
        )
        sys.exit(1)
    version = networkx.__version__
    print("Networkx " + version + " found", file=sys.stderr)
    if not shutil.which("samtools"):
        print(
            time.strftime("%c")
            + ": Samtools does not exist in PATH. Terminating....\n",
            file=sys.stderr,
        )
        sys.exit(1)

    if not shutil.which("bamToBed"):
        print(
            time.strftime("%c")
            + ": Bedtools does not exist in PATH. Terminating....\n",
            file=sys.stderr,
        )
        sys.exit(1)

    args.dir.mkdir(parents=True, exist_ok=True)
    print(time.strftime("%c") + ":Starting scaffolding..", file=sys.stderr)

    alignment_bed = args.dir / "alignment.bed"
    if not alignment_bed.exists():
        print("converting bam file to bed file", file=sys.stderr)
        try:
            with open(alignment_bed, "w") as alignment_bed_f:
                subprocess.run(
                    [
                        "bamToBed",
                        "-i",
                        str(args.mapping),
                    ],
                    stdout=alignment_bed_f,
                    check=True,
                )
            print("finished conversion", file=sys.stderr)
        except subprocess.CalledProcessError as err:
            alignment_bed.unlink()
            print(
                time.strftime("%c")
                + ": Failed in coverting bam file to bed format, terminating scaffolding....\n"
                + str(err.output),
                file=sys.stderr,
            )
            sys.exit(1)

    assembly_idx = args.assembly.parent / (args.assembly.name + ".fai")
    try:
        subprocess.run(
            [
                "samtools",
                "faidx",
                str(args.assembly),
                "--fai-idx",
                str(assembly_idx),
            ],
            check=True,
        )
    except subprocess.CalledProcessError as err:
        print(str(err.output), file=sys.stderr)
        sys.exit(1)

    contig_length = args.dir / "contig_length"
    with open(assembly_idx, "r", newline="") as f_in, open(
        contig_length, "w", newline=""
    ) as f_out:
        reader = csv.DictReader(
            f_in,
            fieldnames=["name", "length", "offset", "linebases", "linewidth"],
            delimiter="\t",
        )
        writer = csv.DictWriter(f_out, fieldnames=["name", "length"], delimiter="\t")
        for row in reader:
            writer.writerow({"name": row["name"], "length": row["length"]})

    print(time.strftime("%c") + ":Finished conversion", file=sys.stderr)

    print(
        time.strftime("%c") + ":Started generating links between contigs",
        file=sys.stderr,
    )
    contig_links = args.dir / "contig_links"
    contig_coverage = args.dir / "contig_coverage"
    if not contig_links.exists():
        try:
            subprocess.run(
                [
                    str(cwd / "libcorrect"),
                    "-a",
                    str(alignment_bed),
                    "-d",
                    str(contig_length),
                    "-o",
                    str(contig_links),
                    "-x",
                    str(contig_coverage),
                    "-c",
                    str(args.length),
                ],
                check=True,
            )
            print(
                time.strftime("%c") + ":Finished generating links between contigs",
                file=sys.stderr,
            )
        except subprocess.CalledProcessError as err:
            contig_links.unlink()
            print(
                time.strftime("%c")
                + ": Failed in generate links from bed file, terminating scaffolding....\n"
                + str(err.output),
                file=sys.stderr,
            )
            sys.exit(1)

    print(
        time.strftime("%c") + ":Started bulding of links between contigs",
        file=sys.stderr,
    )
    bundled_links = args.dir / "bundled_links"
    bundled_graph = args.dir / "bundled_graph.gml"
    if not bundled_links.exists():
        try:
            subprocess.run(
                [
                    str(cwd / "bundler"),
                    "-l",
                    str(contig_links),
                    "-o",
                    str(bundled_links),
                    "-b",
                    str(bundled_graph),
                    "-c",
                    str(args.bsize),
                ],
                check=True,
            )
            print(
                time.strftime("%c") + ":Finished bundling of links between contigs",
                file=sys.stderr,
            )
        except subprocess.CalledProcessError as err:
            bundled_links.unlink()
            bundled_graph.unlink()
            print(
                time.strftime("%c")
                + ": Failed to bundle links, terminating scaffolding....\n"
                + str(err.output),
                file=sys.stderr,
            )
            sys.exit(1)

    invalidated_counts = args.dir / "invalidated_counts"
    high_centrality = args.dir / "high_centrality.txt"
    bundled_links_filtered = args.dir / "bundled_links_filtered"
    oriented_gml = args.dir / "oriented.gml"
    oriented_links = args.dir / "oriented_links"
    repeats = args.dir / "repeats"
    if args.repeats:
        print(
            time.strftime("%c") + ":Started finding and removing repeats",
            file=sys.stderr,
        )
        try:
            subprocess.run(
                [
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
                ],
                check=True,
            )
        except subprocess.CalledProcessError as err:
            print(
                time.strftime("%c")
                + ": Failed to find repeats, terminating scaffolding...\n"
                + str(err.output),
                file=sys.stderr,
            )

        try:
            subprocess.run(
                [
                    "python",
                    str(cwd / "centrality.py"),
                    "-g",
                    str(bundled_links),
                    "-l",
                    str(contig_length),
                    "-o",
                    str(high_centrality),
                ],
                check=True,
            )
        except subprocess.CalledProcessError as err:
            print(
                time.strftime("%c")
                + ": Failed to find repeats, terminating scaffolding....\n"
                + str(err.output),
                file=sys.stderr,
            )
            sys.exit(1)

        try:
            with open(bundled_links_filtered, "w") as bundled_links_filtered_f:
                subprocess.run(
                    [
                        "python",
                        str(cwd / "repeat_filter.py"),
                        str(contig_coverage),
                        str(bundled_links),
                        str(invalidated_counts),
                        str(high_centrality),
                        str(contig_length),
                        str(repeats),
                    ],
                    stdout=bundled_links_filtered_f,
                    check=True,
                )
        except subprocess.CalledProcessError as err:
            print(
                time.strftime("%c")
                + ": Failed to find repeats, terminating scaffolding....\n"
                + str(err.output),
                file=sys.stderr,
            )
            sys.exit(1)
        print(
            time.strftime("%c") + ":Finished repeat finding and removal",
            file=sys.stderr,
        )
    else:
        bundled_links.rename(bundled_links_filtered)
    print(time.strftime("%c") + ":Started orienting the contigs", file=sys.stderr)
    try:
        subprocess.run(
            [
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
            ],
            check=True,
        )
        print(time.strftime("%c") + ":Finished orienting the contigs", file=sys.stderr)
    except subprocess.CalledProcessError:
        print(
            time.strftime("%c")
            + ": Failed to Orient contigs, terminating scaffolding....",
            file=sys.stderr,
        )

    print(time.strftime("%c") + ":Started finding separation pairs", file=sys.stderr)
    seppairs = args.dir / "seppairs"
    try:
        subprocess.run(
            [
                str(cwd / "spqr"),
                "-l",
                str(oriented_links),
                "-o",
                str(seppairs),
            ],
            check=True,
        )
        print(
            time.strftime("%c") + ":Finished finding spearation pairs", file=sys.stderr
        )
    except subprocess.CalledProcessError as err:
        print(
            time.strftime("%c")
            + ": Failed to decompose graph, terminating scaffolding....\n"
            + str(err.output),
            file=sys.stderr,
        )
        sys.exit(1)

    print(time.strftime("%c") + ":Finding the layout of contigs", file=sys.stderr)
    bubbles = args.dir / "bubbles.txt"
    scaffolds_fasta = args.dir / "scaffolds.fa"
    if not scaffolds_fasta.exists():
        try:
            subprocess.run(
                [
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
                    str(scaffolds_fasta),
                    "-f",
                    str(args.dir / "scaffolds.agp"),
                    "-e",
                    str(args.dir / "scaffold_graph.gfa"),
                ],
                check=True,
            )
            print(
                time.strftime("%c") + ":Final scaffolds written, Done!", file=sys.stderr
            )
        except subprocess.CalledProcessError as err:
            print(
                time.strftime("%c")
                + ": Failed to generate scaffold sequences, terminating scaffolding....\n"
                + str(err.output),
                file=sys.stderr,
            )

    if args.visualization:
        # Output the MetagenomeScope .db file directly to args.dir. The only file
        # created by collate.py here is the mgsc.db file.
        subprocess.run(
            [
                "python",
                str(cwd / "MetagenomeScope/graph_collator/collate.py"),
                "-i",
                str(oriented_gml),
                "-w",
                "-ub",
                str(bubbles),
                "-ubl",
                "-d",
                str(args.dir),
                "-o",
                "mgsc",
            ],
            check=True,
        )

    if not args.keep:
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


if __name__ == "__main__":
    main()

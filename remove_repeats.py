import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l","--links",help="bundled links",required=True)
    parser.add_argument("-r","--repeats",help="repeats file",required=True)
    parser.add_argument("-c","--coverage",help="coverage file",required=True)
    parser.add_argument("-o","--output",help="bundled links with repeativie nodes removed",required=True)

    args = parser.parse_args()

    repeats = {}
    coverage = {}
    with open(args.coverage,'r') as f:
        for line in f:
            attrs = line.split()
            coverage[attrs[0]] = float(attrs[1])

    with open(args.repeats,'r') as f:
        for line in f:
            attrs = line.split()
            contig = attrs[0].replace("\"","")
            repeats[contig] = 1

    ofile = open(args.output,'w')

    with open(args.links,'r') as f:
        for line in f:
            attrs = line.split()
            if attrs[0] in repeats or attrs[2] in repeats or coverage[attrs[0]] >= 10*coverage[attrs[1]] or coverage[attrs[1]] >= 10*coverage[attrs[0]]:
                continue
            ofile.write(line)
            
if __name__ == '__main__':
    main()
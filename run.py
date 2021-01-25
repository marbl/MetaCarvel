import os
import argparse
import sys
import time
import subprocess
from subprocess import Popen, PIPE


def cmd_exists(cmd):
    return subprocess.call("type " + cmd, shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

def main():
    cwd=os.path.dirname(os.path.abspath(__file__))

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
    if not cmd_exists('samtools'):
      print(time.strftime("%c")+': Samtools does not exist in PATH. Terminating....\n', file=sys.stderr)
      sys.exit(1)

    if not cmd_exists('bamToBed'):
      print(time.strftime("%c")+': Bedtools does not exist in PATH. Terminating....\n', file=sys.stderr)
      sys.exit(1)

    if not os.path.exists(args.dir):
        os.makedirs(args.dir)
    print(time.strftime("%c")+':Starting scaffolding..', file=sys.stderr)


    if os.path.exists(args.dir+'/alignment.bed') == False:
        print("converting bam file to bed file", file=sys.stderr)
        #os.system('bamToBed -i ' + args.mapping + " > " + args.dir+'/alignment.bed')
        try:
          p = subprocess.check_output('bamToBed -i ' + args.mapping + " > " + args.dir+'/alignment.bed', shell=True)
          print('finished conversion', file=sys.stderr)
        except subprocess.CalledProcessError as err:
          os.system("rm " + args.dir+'/alignment.bed')
          print(time.strftime("%c")+': Failed in coverting bam file to bed format, terminating scaffolding....\n' + str(err.output), file=sys.stderr)
          sys.exit(1)
    try:
      #os.system('samtools faidx '+args.assembly)
      p = subprocess.check_output('samtools faidx '+args.assembly,shell=True)
    except subprocess.CalledProcessError as err:
      print(str(err.output), file=sys.stderr)
      sys.exit()
    os.system('cut -f 1,2 '+ args.assembly+'.fai > '+args.dir+'/contig_length')

    print(time.strftime("%c")+':Finished conversion', file=sys.stderr)

    final_assembly = args.assembly
    final_mapping = args.mapping

    print(time.strftime("%c") + ':Started generating links between contigs', file=sys.stderr)
    if os.path.exists(args.dir+'/contig_links') == False:
        #print './libcorrect -l' + args.lib + ' -a' + args.dir+'/alignment.bed -d ' +args.dir+'/contig_length -o '+ args.dir+'/contig_links'
        try:
          #os.system('./libcorrect -l ' + args.lib + ' -a ' + args.dir+'/alignment.bed -d ' +args.dir+'/contig_length -o '+ args.dir+'/contig_links -x '+args.dir+'/contig_coverage')
           p = subprocess.check_output(cwd+'/libcorrect -a ' + args.dir+'/alignment.bed -d ' +args.dir+'/contig_length -o '+ args.dir+'/contig_links -x '+args.dir+'/contig_coverage -c '+str(args.length),shell=True)
           print(time.strftime("%c") +':Finished generating links between contigs', file=sys.stderr)
        except subprocess.CalledProcessError as err:
            os.system('rm '+args.dir+'/contig_links')
            print(time.strftime("%c")+': Failed in generate links from bed file, terminating scaffolding....\n' + str(err.output), file=sys.stderr)
            sys.exit(1)

    print(time.strftime("%c")+':Started bulding of links between contigs', file=sys.stderr)
    if os.path.exists(args.dir+'/bundled_links') == False:
        try:
          #os.system('./bundler -l '+ args.dir+'/contig_links -o ' + args.dir+'/bundled_links + -b '+args.dir+'/bundled_graph.gml')
          p = subprocess.check_output(cwd+'/bundler -l '+ args.dir+'/contig_links -o ' + args.dir+'/bundled_links + -b '+args.dir+'/bundled_graph.gml -c '+str(args.bsize), shell=True)
          print(time.strftime("%c")+':Finished bundling of links between contigs', file=sys.stderr)
        except subprocess.CalledProcessError as err:
          os.system('rm '+args.dir+'/bundled_links')
          os.system('rm '+args.dir+'/bundled_graph.gml')
          print(time.strftime("%c")+': Failed to bundle links, terminating scaffolding....\n' + str(err.output), file=sys.stderr)
          sys.exit(1)

    if args.repeats == "true":
        print(time.strftime("%c")+':Started finding and removing repeats', file=sys.stderr)
        try:
            p = subprocess.check_output(cwd+'/orientcontigs -l '+args.dir+'/bundled_links -c '+ args.dir+'/contig_length --bsize -o ' +args.dir+'/oriented.gml -p ' + args.dir+'/oriented_links -i '+args.dir+'/invalidated_counts',shell=True)

        except subprocess.CalledProcessError as err:
            print(time.strftime("%c") + ': Failed to find repeats, terminating scaffolding...\n' + str(err.output), file=sys.stderr)

        try:
            p = subprocess.check_output('python '+cwd+'/centrality.py  -g '+args.dir+'/bundled_links -l ' + args.dir+ '/contig_length -o  '+args.dir+'/high_centrality.txt' ,shell=True)
        except subprocess.CalledProcessError as err:
                print(time.strftime("%c")+': Failed to find repeats, terminating scaffolding....\n' + str(err.output), file=sys.stderr)
                sys.exit(1)

        try:
            p = subprocess.check_output('python '+cwd+'/repeat_filter.py  '+args.dir+'/contig_coverage ' + args.dir+ '/bundled_links ' + args.dir+'//invalidated_counts ' + args.dir+'/high_centrality.txt ' + args.dir+ '/contig_length '+ args.dir+'/repeats > ' + args.dir+'//bundled_links_filtered',shell=True)
        except subprocess.CalledProcessError as err:
            print(time.strftime("%c")+': Failed to find repeats, terminating scaffolding....\n' + str(err.output), file=sys.stderr)
            sys.exit(1)
        print(time.strftime("%c")+':Finished repeat finding and removal', file=sys.stderr)
    else:
        os.system('mv '+args.dir+'/bundled_links ' + args.dir+'/bundled_links_filtered')
    print(time.strftime("%c")+':Started orienting the contigs', file=sys.stderr)
    # if os.path.exists(args.dir+'/oriented_links') == False:
      #os.system('./orientcontigs -l '+args.dir+'/bundled_links_filtered -c '+ args.dir+'/contig_length --bsize -o ' +args.dir+'/oriented.gml -p ' + args.dir+'/oriented_links' )
    try:
        p = subprocess.check_output(cwd+'/orientcontigs -l '+args.dir+'/bundled_links_filtered -c '+ args.dir+'/contig_length --bsize -o ' +args.dir+'/oriented.gml -p ' + args.dir+'/oriented_links -i '+args.dir+'/invalidated_counts',shell=True)
        print(time.strftime("%c")+':Finished orienting the contigs', file=sys.stderr)
    except subprocess.CalledProcessError:
        print(time.strftime("%c")+': Failed to Orient contigs, terminating scaffolding....', file=sys.stderr)

    print(time.strftime("%c")+':Started finding separation pairs', file=sys.stderr)
    #if os.path.exists(args.dir+'/seppairs') == False:
    #os.system('./spqr -l ' + args.dir+'/oriented_links -o ' + args.dir+'/seppairs')
    try:
        p = subprocess.check_output(cwd+'/spqr -l ' + args.dir+'/oriented_links -o ' + args.dir+'/seppairs',shell=True)
        print(time.strftime("%c")+':Finished finding spearation pairs', file=sys.stderr)
    except subprocess.CalledProcessError as err:
        print(time.strftime("%c")+': Failed to decompose graph, terminating scaffolding....\n' + str(err.output), file=sys.stderr)
        sys.exit(1)

    print(time.strftime("%c")+':Finding the layout of contigs', file=sys.stderr)
    if os.path.exists(args.dir+'/scaffolds.fasta') == False:
        try:
            p = subprocess.check_output('python '+cwd+'/layout.py -a '+ args.assembly +' -b '+args.dir+'/bubbles.txt' +' -g ' + args.dir+'/oriented.gml -s '+args.dir+'/seppairs -o '+args.dir+'/scaffolds.fa -f '+args.dir+'/scaffolds.agp -e '+args.dir+'/scaffold_graph.gfa',shell=True)
            print(time.strftime("%c")+':Final scaffolds written, Done!', file=sys.stderr)
        except subprocess.CalledProcessError as err:
            print(time.strftime("%c")+': Failed to generate scaffold sequences, terminating scaffolding....\n' + str(err.output), file=sys.stderr)

    if args.visualization == "true":
        #try:
      graphpath = os.path.abspath(args.dir+'/oriented.gml')
      bubblepath = os.path.abspath(args.dir+'/bubbles.txt')
      # Output the MetagenomeScope .db file directly to args.dir. The only file
      # created by collate.py here is the mgsc.db file.
      os.system('python '+cwd+'/MetagenomeScope/graph_collator/collate.py -i '
              + graphpath + ' -w -ub ' + bubblepath + ' -ubl -d ' + args.dir
              + ' -o mgsc')
          #p = subprocess.check_output('python '+cwd+'/MetagenomeScope/graph_collator/collate.py -i ' + graphpath + ' -w -ub ' + bubblepath + ' -ubl -d ' + args.dir + ' -o mgsc')
        #except subprocess.CalledProcessError as err:
            #print >> sys.stderr, time.strftime("%c")+": Failed to run MetagenomeScope \n" + str(err.output)

    if not args.keep == "true":
      if os.path.exists(args.dir+'/contig_length'):
       os.system("rm "+args.dir+'/contig_length')
      if os.path.exists(args.dir+'/contig_links'):
       os.system("rm "+args.dir+'/contig_links')
      if os.path.exists(args.dir+'/contig_coverage'):
        os.system("rm "+args.dir+'/contig_coverage')
      if os.path.exists(args.dir+'/bundled_links'):
        os.system("rm "+args.dir+'/bundled_links')
      if os.path.exists(args.dir+'/bundled_links_filtered'):
        os.system("rm "+args.dir+'/bundled_links_filtered')
      if os.path.exists(args.dir+'/bundled_graph.gml'):
        os.system("rm "+args.dir+'/bundled_graph.gml')
      if os.path.exists(args.dir+'/invalidated_counts'):
        os.system("rm "+args.dir+'/invalidated_counts')
      if os.path.exists(args.dir+'/repeats'):
        os.system("rm "+args.dir+'/repeats')
      if os.path.exists(args.dir+'/oriented_links'):
        os.system("rm "+args.dir+'/oriented_links')
      if os.path.exists(args.dir+'/oriented.gml'):
        os.system("rm "+args.dir+'/oriented.gml')
      if os.path.exists(args.dir+'/seppairs'):
        os.system("rm "+args.dir+'/seppairs')
      if os.path.exists(args.dir+'/alignment.bed'):
        os.system("rm "+args.dir+'/alignment.bed')
if __name__ == '__main__':
    main()

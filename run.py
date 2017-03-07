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
    bin=os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser(description="Bambus3: A scaffolding tool for metagenomic assemblies")
    parser.add_argument("-a","--assembly",help="assembled contigs",required=True)
    parser.add_argument("-m","--mapping", help="mapping of read to contigs in bam format",required=True)
    parser.add_argument("-d","--dir",help="output directory for results",default='out',required=True)
    parser.add_argument("-f",'--force',help="force re-run of pipeline, will remove any existing output",default=False)
    parser.add_argument("-r",'--repeats',help="To turn repeat detection on",default=False)
    parser.add_argument("-k","--keep", help="Set this to keep temporary files in output directory",default=False)
    parser.add_argument("-l","--length",help="Minimum length of contigs to consider for scaffolding",default=500)
    parser.add_argument("-b","--bsize",help="Minimum mate pair support between contigs to consider for scaffolding",default=3)

    args = parser.parse_args()
    try:
      import networkx
    except ImportError:
      raise ImportError('Looks like you do not have networkx. Please rerun with networkx module installed.')
      sys.exit(1)
    if not cmd_exists('samtools'):
      print >> sys.stderr, time.strftime("%c")+': Samtools does not exist in PATH. Terminating....\n'
      sys.exit(1)

    if not cmd_exists('bamToBed'):
      print >> sys.stderr, time.strftime("%c")+': Bedtools does not exist in PATH. Terminating....\n'
      sys.exit(1)

    if not os.path.exists(args.dir):
        os.makedirs(args.dir)
    print >> sys.stderr, time.strftime("%c")+':Starting scaffolding..'

    if args.force:
       if os.path.exists(args.dir+'/alignment.bed'): os.unlink(args.dir+'/alignment.bed')
       if os.path.exists(args.dir+'/contig_length'): os.unlink(args.dir+'/contig_length')
       if os.path.exists(args.dir+'/contig_links'): os.unlink(args.dir+'/contig_links')
       if os.path.exists(args.dir+"/bundled_links"): os.unlink(args.dir+"/bundled_links")
       if os.path.exists(args.dir+'/oriented_links'): os.unlink(args.dir+'/oriented_links')
       if os.path.exists(args.dir+'/seppairs'): os.unlink(args.dir+'/seppairs')
       if os.path.exists(args.dir+'/scaffolds.fasta'): os.unlink(args.dir+'/scaffolds.fasta')

    if os.path.exists(args.dir+'/alignment.bed') == False:
        print >> sys.stderr, "converting bam file to bed file"
        #os.system('bamToBed -i ' + args.mapping + " > " + args.dir+'/alignment.bed')
        try:
          p = subprocess.check_output('bamToBed -i ' + args.mapping + " > " + args.dir+'/alignment.bed', shell=True)
          print >> sys.stderr, 'finished conversion'
        except subprocess.CalledProcessError as err:
          os.system("rm " + args.dir+'/alignment.bed')
          print >> sys.stderr, time.strftime("%c")+': Failed in coverting bam file to bed format, terminating scaffolding....\n' + str(err.output)
          sys.exit(1)
    try:
      #os.system('samtools faidx '+args.assembly)
      p = subprocess.check_output('samtools faidx '+args.assembly,shell=True)
    except subprocess.CalledProcessError as err:
      print >> sys.stderr, str(err.output)
      sys.exit()
    os.system('cut -f 1,2 '+ args.assembly+'.fai > '+args.dir+'/contig_length')

    print >> sys.stderr, time.strftime("%c")+':Finished conversion'

    final_assembly = args.assembly
    final_mapping = args.mapping

    print >> sys.stderr, time.strftime("%c") + ':Started generating links between contigs'
    if os.path.exists(args.dir+'/contig_links') == False:
        #print './libcorrect -l' + args.lib + ' -a' + args.dir+'/alignment.bed -d ' +args.dir+'/contig_length -o '+ args.dir+'/contig_links'
        try:
          #os.system('./libcorrect -l ' + args.lib + ' -a ' + args.dir+'/alignment.bed -d ' +args.dir+'/contig_length -o '+ args.dir+'/contig_links -x '+args.dir+'/contig_coverage')
           p = subprocess.check_output('./libcorrect -a ' + args.dir+'/alignment.bed -d ' +args.dir+'/contig_length -o '+ args.dir+'/contig_links -x '+args.dir+'/contig_coverage -c '+str(args.length),shell=True)
           print >> sys.stderr, time.strftime("%c") +':Finished generating links between contigs'
        except subprocess.CalledProcessError as err:
            os.system('rm '+args.dir+'/contig_links')
            print >> sys.stderr, time.strftime("%c")+': Failed in generate links from bed file, terminating scaffolding....\n' + str(err.output)
            sys.exit(1)

    print >> sys.stderr, time.strftime("%c")+':Started bulding of links between contigs'
    if os.path.exists(args.dir+'/bundled_links') == False:
        try:
          #os.system('./bundler -l '+ args.dir+'/contig_links -o ' + args.dir+'/bundled_links + -b '+args.dir+'/bundled_graph.gml')
          p = subprocess.check_output('./bundler -l '+ args.dir+'/contig_links -o ' + args.dir+'/bundled_links + -b '+args.dir+'/bundled_graph.gml -c '+str(args.bsize), shell=True)
          print >> sys.stderr, time.strftime("%c")+':Finished bundling of links between contigs'
        except subprocess.CalledProcessError as err:
          os.system('rm '+args.dir+'/bundled_links')
          os.system('rm '+args.dir+'/bundled_graph.gml')
          print >> sys.stderr, time.strftime("%c")+': Failed to bundle links, terminating scaffolding....\n' + str(err.output)
          sys.exit(1)

    if args.repeats:
      print >> sys.stderr, time.strftime("%c")+':Started finding and removing repeats'
      try:
       #os.system('./vc_algo -g '+ args.dir+'/bundled_graph.gml -r ' + args.dir+'/repeats -f 0.02')
        p = subprocess.check_output('./centrality -g '+ args.dir+'/bundled_graph.gml -r ' + args.dir+'/repeats -f 0.02',shell=True)

      except subprocess.CalledProcessError as err:
        print >> sys.stderr, time.strftime("%c")+': Failed to find repeats, terminating scaffolding....\n' + str(err.output)
        sys.exit(1)
       #print './vc_algo -g '+ args.dir+'/bundled_graph.gml -r ' + args.dir+'/repeats -f 0.02'
      try:
        p = subprocess.check_output('./orientcontigs -l '+args.dir+'/bundled_links -c '+ args.dir+'/contig_length --bsize -o ' +args.dir+'/oriented.gml -p ' + args.dir+'/oriented_links -i '+args.dir+'/invalidated_counts',shell=True)
      
      except subprocess.CalledProcessError as err:
        print >> sys.stderr, time.strftime("%c") + ': Failed to find repeats, terminating scaffolding...\n' + str(err.output)

      try:
        #os.system('python remove_repeats.py -l '+ args.dir+'/bundled_links -r ' + args.dir+'/repeats -o '+ args.dir+'/bundled_links_filtered -c '+args.dir+'/contig_coverage -s '+args.dir+'/contig_length')
        p = subprocess.check_output('python repeat_filter.py '+args.dir+'/contig_coverage ' + args.dir+ '/bundled_links ' +args.dir+'/invalidated_counts ' + args.dir+'/repeats ' + args.dir+'/contig_length > '+args.dir+'/bundled_links_filtered' ,shell=True)
      except subprocess.CalledProcessError as err:
        print >> sys.stderr, time.strftime("%c")+': Failed to find repeats, terminating scaffolding....\n' + str(err.output)
        sys.exit(1)
      print >> sys.stderr, time.strftime("%c")+':Finished repeat finding and removal'
    else:
      os.system('mv '+args.dir+'/bundled_links ' + args.dir+'/bundled_links_filtered')

    print >> sys.stderr, time.strftime("%c")+':Started orienting the contigs'
    # if os.path.exists(args.dir+'/oriented_links') == False:
      #os.system('./orientcontigs -l '+args.dir+'/bundled_links_filtered -c '+ args.dir+'/contig_length --bsize -o ' +args.dir+'/oriented.gml -p ' + args.dir+'/oriented_links' ) 
    try:
      p = subprocess.check_output('./orientcontigs -l '+args.dir+'/bundled_links_filtered -c '+ args.dir+'/contig_length --bsize -o ' +args.dir+'/oriented.gml -p ' + args.dir+'/oriented_links -i '+args.dir+'/invalidated_counts',shell=True)
      print >> sys.stderr, time.strftime("%c")+':Finished orienting the contigs'
    except subprocess.CalledProcessError:
      print >> sys.stderr, time.strftime("%c")+': Failed to Orient contigs, terminating scaffolding....'
      
    print >> sys.stderr, time.strftime("%c")+':Started finding separation pairs'
    #if os.path.exists(args.dir+'/seppairs') == False:
    #os.system('./spqr -l ' + args.dir+'/oriented_links -o ' + args.dir+'/seppairs')
    try:
      p = subprocess.check_output('./spqr -l ' + args.dir+'/oriented_links -o ' + args.dir+'/seppairs',shell=True)
      print >> sys.stderr, time.strftime("%c")+':Finished finding spearation pairs'
    except subprocess.CalledProcessError as err:
      print >> sys.stderr, time.strftime("%c")+': Failed to decompose graph, terminating scaffolding....\n' + str(err.output)
      sys.exit(1)
    
    print >> sys.stderr, time.strftime("%c")+':Finding the layout of contigs'
    if os.path.exists(args.dir+'/scaffolds.fasta') == False:
      #os.system('python layout.py -a '+ args.assembly + ' -g ' + args.dir+'/oriented.gml -s '+args.dir+'/seppairs -o '+args.dir+'/scaffolds.fa -f '+args.dir+'/scaffolds.agp -b '+args.dir+'/bubbles')
      try:
        p = subprocess.check_output('python layout.py -a '+ args.assembly + ' -g ' + args.dir+'/oriented.gml -s '+args.dir+'/seppairs -o '+args.dir+'/scaffolds.fa -f '+args.dir+'/scaffolds.agp -e '+args.dir+'/scaffold_graph.gfa',shell=True)
        print >> sys.stderr, time.strftime("%c")+':Final scaffolds written, Done!'
      except subprocess.CalledProcessError as err:
        print >> sys.stderr, time.strftime("%c")+': Failed to generate scaffold sequenes , terminating scaffolding....\n' + str(err.output)

    if not args.keep:
      os.system("rm "+args.dir+'/contig_length')
      os.system("rm "+args.dir+'/contig_links')
      os.system("rm "+args.dir+'/contig_coverage')
      os.system("rm "+args.dir+'/bundled_links')
      os.system("rm "+args.dir+'/bundled_links_filtered')
      os.system("rm "+args.dir+'/bundled_graph.gml')
      os.system("rm "+args.dir+'/invalidated_counts')
      os.system("rm "+args.dir+'/repeats')
      os.system("rm "+args.dir+'/oriented_links')
      os.system("rm "+args.dir+'/oriented.gml')
      os.system("rm "+args.dir+'/seppairs')
      os.system("rm "+args.dir+'/alignment.bed')
if __name__ == '__main__':
    main()

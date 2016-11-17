import os
import argparse
import sys
import time

def main():
    bin=os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser()
    parser.add_argument("-a","--assembly",help="assembled contigs",required=True)
    parser.add_argument("-m","--mapping", help="mapping of read to contigs in bam format",required=True)
    parser.add_argument("-d","--dir",help="output directory for results",default='out',required=True)
    parser.add_argument("-l","--lib",help="library description for mate pairs",default='out',required=True)
    parser.add_argument("-f",'--force',help="force re-run of pipeline, will remove any existing output",default=False)

    args = parser.parse_args()

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
        print >> sys.stderr, "convertng bam file to bed file"
        os.system('bamToBed -i ' + args.mapping + " > " + args.dir+'/alignment.bed')
        print >> sys.stderr, 'finished conversion'
    os.system('samtools faidx '+args.assembly)
    os.system('cut -f 1,2 '+ args.assembly+'.fai > '+args.dir+'/contig_length')

    print >> sys.stderr, time.strftime("%c")+':Finished conversion'

    final_assembly = args.assembly
    final_mapping = args.mapping

    print >> sys.stderr, time.strftime("%c") + ':Started generating links between contigs'
    if os.path.exists(args.dir+'/contig_links') == False:
       	#print './libcorrect -l' + args.lib + ' -a' + args.dir+'/alignment.bed -d ' +args.dir+'/contig_length -o '+ args.dir+'/contig_links'
	os.system('./libcorrect -l ' + args.lib + ' -a ' + args.dir+'/alignment.bed -d ' +args.dir+'/contig_length -o '+ args.dir+'/contig_links -x '+args.dir+'/contig_coverage')
    print >> sys.stderr, time.strftime("%c") +':Finished generating links between contigs'

    print >> sys.stderr, time.strftime("%c")+':Started bulding of links between contigs'
    if os.path.exists(args.dir+'/bundled_links') == False:
       os.system('./bundler -l '+ args.dir+'/contig_links -o ' + args.dir+'/bundled_links + -b '+args.dir+'/bundled_graph.gml')
    print >> sys.stderr, time.strftime("%c")+':Finished bunlding of links between contigs'

    print >> sys.stderr, time.strftime("%c")+':Started finding and removing repeats'
    if os.path.exists(args.dir+'/repeats') == False:
       os.system('./vc_algo -g '+ args.dir+'/bundled_graph.gml -r ' + args.dir+'/repeats -f 0.02')
       #print './vc_algo -g '+ args.dir+'/bundled_graph.gml -r ' + args.dir+'/repeats -f 0.02'
       os.system('python remove_repeats.py -l '+ args.dir+'/bundled_links -r ' + args.dir+'/repeats -o '+ args.dir+'/bundled_links_filtered')
    print >> sys.stderr, time.strftime("%c")+':Finished repeat finding and removal'

    print >> sys.stderr, time.strftime("%c")+':Started orienting the contigs'
    if os.path.exists(args.dir+'/oriented_links') == False:
       os.system('./orientcontigs -l '+args.dir+'/bundled_links_filtered -c '+ args.dir+'/contig_length --bsize -o ' +args.dir+'/oriented.gml -p ' + args.dir+'/oriented_links' ) 
    print >> sys.stderr, time.strftime("%c")+':Finished orienting the contigs'
    
    print >> sys.stderr, time.strftime("%c")+':Started finding separation pairs'
    if os.path.exists(args.dir+'/seppairs') == False:
      os.system('./spqr -l ' + args.dir+'/oriented_links -o ' + args.dir+'/seppairs')
    print >> sys.stderr, time.strftime("%c")+':Finished finding spearation pairs'

    print >> sys.stderr, time.strftime("%c")+':Finding the layout of contigs'
    if os.path.exists(args.dir+'/scaffolds.fasta') == False:
      os.system('python layout.py -a '+ args.assembly + ' -g ' + args.dir+'/oriented.gml -s '+args.dir+'/seppairs -o '+args.dir+'/scaffolds.fa')
    print >> sys.stderr, time.strftime("%c")+':Final scaffolds written, Done!'
if __name__ == '__main__':
    main()

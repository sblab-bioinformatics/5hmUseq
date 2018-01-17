#!/usr/bin/env python

import re
import sys
import argparse
import os
parser = argparse.ArgumentParser(description= """
  DESCRIPTION
  Extract barcodes and inserts from fastq file based on given tags marking the position of barcodes and 
  inxsert.
""", formatter_class= argparse.RawTextHelpFormatter, prog= os.path.basename(__file__))

parser.add_argument('input', 
help='''Input fastq, might be from stdin
''')

parser.add_argument('--left_tag', '-l', 
default= 'ATCGAG',
help='''Sequence of the left tag. Barcode is taken before it.
''')

parser.add_argument('--right_tag', '-r', 
default= 'CGTGTC',
help='''Sequence of the left tag. Barcode is taken before it.
''')

parser.add_argument('--barcode_len', '-bl', 
default= 10,
type= int,
help='''Length of barcodes.
''')

parser.add_argument('--version', action='version', version='%(prog)s 0.1.0')

def getInsertAndBarcodes(xseq, left_tag= 'ATCGAG', right_tag= 'CGTGTC', blen= 10):
    """ Find the position of the left_tag and extract the barcode that comes before it. The same for 
    the right_tag: Find it and extract the barcode on the left. Finally get the insert between left_tag
    and barcode on the right. 
    Return the indexes of insert, and barcodes.
    Test:
    xseq= 'NNaxxxxxxxxzATCGAGAATCCCGGTGCCGATACCHACTCTTGHAGAZayyyyyyyyzCGTGTCAGATATATACATCCGAT'
    x= getInsertAndBarcodes(xseq)
    xseq[x['ins'][0] : x['ins'][1]] == 'ATCGAGAATCCCGGTGCCGATACCHACTCTTGHAGAZ'
    xseq[x['l_bar'][0] : x['l_bar'][1]] == 'axxxxxxxxz'
    xseq[x['r_bar'][0] : x['r_bar'][1]] == 'ayyyyyyyyz'

    ## Missing tags:
    getInsertAndBarcodes('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx') == None
    ## Insert too short (right barcode overlaps left barcode)
    getInsertAndBarcodes('NNaxxxxxxxxzATCGAGyyyzCGTGTCAGATATATACATCCGAT') == None
    """
    xseq= xseq.upper(); left_tag= left_tag.upper(); right_tag= right_tag.upper()
    if left_tag not in xseq or right_tag not in xseq:
        return None
    l_idx= re.search(left_tag, xseq).start()
    left= l_idx-blen
    if left < 0:
        left= 0
    l_bar= (left, l_idx)

    r_idx= re.search(right_tag, xseq).start()
    r_bar= (r_idx-blen, r_idx)

    if l_bar[1] >= r_bar[0]:
        return None
    ins= (l_bar[1], r_bar[0])
    return {'ins': ins, 'l_bar': l_bar, 'r_bar': r_bar}

args= parser.parse_args()

if args.input == '-':
    fin= sys.stdin
else:
    fin= open(args.input)

fq= []
for line in fin:
    line= line.strip()
    fq.append(line)
    if len(fq) == 4:
        xseq= fq[1]
        qual= fq[3]
        fq= []
        idx= getInsertAndBarcodes(xseq, args.left_tag, args.right_tag, args.barcode_len)
        if idx is None:
            continue
        rname= '@' + xseq[idx['l_bar'][0] : idx['l_bar'][1]] + ':'  + xseq[idx['r_bar'][0] : idx['r_bar'][1]]
        ins= xseq[idx['ins'][0] : idx['ins'][1]]
        insq= qual[idx['ins'][0] : idx['ins'][1]]
        print '\n'.join([rname, ins, '+', insq])
        
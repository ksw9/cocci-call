#!/usr/bin/env python3

"""
This script takes in a genome fasta file and a seqid2taxid.map file and converts sequences ids to be compatible with Kraken2
"""

### ---------------------------------------- ###

def parseArgs():
    
    # Load fasta file
    fasta_path = argv[argv.index("--fasta") + 1]
    fasta = open(fasta_path).read().split('\n')

    # Load seqid2taxid
    seqid2taxid_path = argv[argv.index("--seqid2taxid") + 1]
    seqid2taxid = open(seqid2taxid_path).read().split('\n')
    seqid2taxid = {st.split('\t')[0] : st.split('\t')[1] for st in seqid2taxid if len(st)}
    
    return fasta, seqid2taxid

### ---------------------------------------- ###

def parseFasta(fasta, seqid2taxid):
    
    parsed_fasta = []
    for f in fasta:
        
        if not len(f):
            
            continue
        
        elif f[0] == ">" and "kraken:taxid" not in f:
            
            seqid = f[1:].replace(' ', '').replace('|', '')
            taxid = seqid2taxid[seqid]
            new_name = f'>{seqid}|kraken:taxid|{taxid} {seqid} |'
            parsed_fasta.append(new_name)
        
        elif f[0] == ">" and "kraken:taxid" in f:
            
            seqid = f[1:].replace(' ', '')
            seqid = seqid if seqid[-1] != '|' else seqid[:-1]
            taxid = seqid2taxid[seqid]
            simple_seqid = seqid.replace('kraken:taxid', '').replace(taxid, '').replace('|', '')
            new_name = f'>{simple_seqid}|kraken:taxid|{taxid} {simple_seqid} |'
            parsed_fasta.append(new_name)
        
        else:
            
            parsed_fasta.append(f)
    
    parsed_fasta = '\n'.join(parsed_fasta)
    
    with open('corrected_fasta.fna', 'w') as fasta_out:
        
        fasta_out.write(parsed_fasta)

### ------------------MAIN------------------ ###

from sys import argv

fasta, seqid2taxid = parseArgs()

parseFasta(fasta, seqid2taxid)

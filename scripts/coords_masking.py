#!/usr/bin/env python3

"""
This script parses a coords file from Nucmer and a fasta file to mask the latter
"""

### ---------------------------------------- ###

def parse_args():

    print(f'args = {argv}')
    
    # Fasta
    fasta_file = argv[argv.index("--fasta_file") + 1]
    
    # Coords
    coords_file = argv[argv.index("--coords_file") + 1]

    return fasta_file, coords_file

### ---------------------------------------- ###

def load_fasta(path):
    
    contig_sequences = {}
    
    with open(path) as fasta:
        
        for f in fasta:
            
            f = f.replace('\n', '')
            
            if f.startswith('>'):
                
                try:
                    
                    contig_sequences[contig_name] = contig_seq
                
                except:
                    
                    pass
                
                contig_name = f.split(' ')[0].replace('>', '')
                contig_seq = ''
                
            else:
                
                contig_seq += f
        
        if contig_name not in contig_sequences.keys():
            
            contig_sequences[contig_name] = contig_seq
                
    return contig_sequences

### ---------------------------------------- ###

def load_coords(path):
    
    masked_regions = []
    
    with open(path) as coords:
        
        for c in coords:
            
            c = c.replace('\n', '').split('|')
            
            if len(c) != 5:
                
                continue
                
            ref, query, _, _, tag = c
            
            s1, e1 = [r for r in ref.split(' ') if len(r)]
            s2, e2 = [q for q in query.split(' ') if len(q)]
            contig = [t for t in tag.replace('\t', ' ').split(' ') if len(t)][0]
                
            if contig != '[TAGS]' and s1 != s2 and e1 != e2:
                
                masked_regions.append([contig, int(s1) - 1, int(e1) - 1])
        
    masked_regions.sort()
        
    return masked_regions

### ---------------------------------------- ###

def resolve_overlaps(raw_data):
    
    clean_data = []
    
    for contig in set([rd[0] for rd in raw_data]):
        
        # Paint contig regions from raw_data
        contig_good_regions = [0] * max([rd[2] for rd in raw_data if rd[0] == contig])
        
        for (c,start,stop) in [rd for rd in raw_data if rd[0] == contig]:
            
            for i in range(start, stop):
                
                contig_good_regions[i] = 1
        
        # Get bed intervals
        end = -1
        for n,s in enumerate(contig_good_regions):
            
            if s == 1 and end == -1:
                
                start = n
                end = n
                
            elif s == 1 and end != -1:
                
                end = n
                
            elif s != 1 and end != -1:
                
                end += 1
                clean_data.append([contig, start, end])
                end = -1
                
            else:
                
                continue
    
    return clean_data

### ---------------------------------------- ###

def mask_seqs(con_seq, mas_reg):
    
    masked_con_seq = {}
    for con, seq in con_seq.items():
        
        mask = [mr[1:] for mr in mas_reg if mr[0] == con]
        
        masked_seq = ''
        for m in mask:
            
            masked_seq += seq[len(masked_seq) : m[0]]
            masked_seq += 'N' * (m[1] - m[0])
            
            if len(masked_seq) != m[1]: break ### TESTING
        
        masked_seq += seq[len(masked_seq):]
        
        masked_con_seq[con] = masked_seq
    
    return masked_con_seq    

### ---------------------------------------- ###

def save_fasta(seq_dict, output_name, chars_per_line=70):
    
    with open(output_name, 'w') as fasta_out:
        
        for contig, seq in seq_dict.items():
            
            fasta_out.write(f'>{contig}\n')
            
            for i in range(0, len(seq), chars_per_line):
                
                fasta_out.write(f'{seq[i : i + chars_per_line]}\n')
            
            if i + chars_per_line < len(seq):
                
                i += 1
                fasta_out.write(f'{seq[i : i + chars_per_line]}\n')

### ------------------MAIN------------------ ###

from sys import argv

### IMPORT DATA ---------------------------- ###

# Parse cli
fasta_file, coords_file = parse_args()

# Load fasta
contig_sequences = load_fasta(fasta_file)

# Load coords file and find regions that need masking
masked_regions = load_coords(coords_file)

# Resolve overlaps of masked regions
masked_regions = resolve_overlaps(masked_regions)

### MASK FASTA ----------------------------- ###

# Mask fasta
masked_contig_sequeces = mask_seqs(contig_sequences, masked_regions)

# Save fasta to file
output_name = 'masked.fasta'
save_fasta(masked_contig_sequeces, output_name, 70)

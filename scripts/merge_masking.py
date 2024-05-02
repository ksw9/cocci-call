#!/usr/bin/env python3

"""
This script parses a coords file from Nucmer and a fasta file to mask the latter
"""

### ---------------------------------------- ###

def parse_args():

    print(f'args = {argv}')
    
    # Fasta
    fasta_files = argv[argv.index("--fasta_files") + 1].split(';')
    
    # Load fasta files
    fastas = [load_fasta(ff) for ff in fasta_files]

    return fastas

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

def merge_fasta_masking(fasta1, fasta2):
    
    merged = {contig : ''.join([f1 if (f1 != 'N') & (f2 != 'N')
                                else 'N'
                                for f1,f2 in zip(fasta1[contig], fasta2[contig])])
              for contig in fasta1.keys()}
    
    return merged

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

# Parse cli and load fasta files
fasta_1, fasta_2 = parse_args()

### MERGE FASTA MASKING -------------------- ###

# Mask fasta
final_fasta = merge_fasta_masking(fasta_1, fasta_2)

# Save fasta to file
output_name = 'merged.fasta'
save_fasta(final_fasta, output_name, 70)

import argparse
import itertools
import subprocess
import regex as re
import pandas as pd
from pyfaidx import Fasta
from Bio.Seq import Seq
import pickle

#PAM, start_editwindow, end_editwindow
BE_library = {}
BE_library['VQR'] = ['.GA', 4, 8]
BE_library['R33A'] = ['.GG', 5, 7]
BE_library['R33A_K34A'] = ['.GG', 5, 6]

#find the sgRNA given the region/row in gtf file 
def find_sgRNA(row, BE, BE_lib, length_5end, total_trans_len):
    chr = row[0]
    start = row[2]
    end = row[3]
    genome = Fasta('/storage/cynthiawu/iSeq/hg38.fa')
    seq = genome[chr][int(start):int(end)].seq
    revcomp_seq = str(Seq(seq).reverse_complement())
    BE_info = BE_lib[BE]
    sgRNA = []
    match = list(re.finditer(BE_info[0], seq[20:], overlapped=True))
    match_rev = list(re.finditer(BE_info[0], revcomp_seq[20:], overlapped=True))
    total_len = len(seq)
    for m in match:
        index = m.start()
        window_start = BE_info[1]
        window_end = BE_info[2]
        window = seq[window_start+index-1: window_end+index]
        if 'C' in window:
            row_copy = row.copy()
            row_copy.append(int(start)+index+1)
            row_copy.append(seq[index: index+23])
            row_copy.append(seq[index: index+23])
            position_edits = [(index+window_start) for index, value in enumerate(window) if value == 'C']
            row_copy.append(len(position_edits))
            row_copy.append(position_edits)
            row_copy.append('+')
            distance_5end = index + length_5end +1
            distance_3end = total_trans_len - distance_5end
            row_copy.append(distance_5end)
            row_copy.append(distance_3end)
            sgRNA.append(row_copy)
    for m in match_rev:
        index = m.start()
        window_start = BE_info[1]
        window_end = BE_info[2]
        window = revcomp_seq[window_start+index-1: window_end+index]
        if 'C' in window:
            row_copy = row.copy()
            row_copy.append(int(end)-index-23+1)
            guide = revcomp_seq[index: index+23]
            row_copy.append(guide)
            row_copy.append(str(Seq(guide).reverse_complement()))
            position_edits = [(index+window_start) for index, value in enumerate(window) if value == 'C']
            row_copy.append(len(position_edits))
            row_copy.append(position_edits)
            row_copy.append('-')
            distance_5end = total_len - index + length_5end +1
            distance_3end = total_trans_len - distance_5end
            row_copy.append(distance_5end)
            row_copy.append(distance_3end)
            sgRNA.append(row_copy)
    return sgRNA
    
#get the region (exon) of the gene from genes_loc pickle and use it to call tabix to get from the indexed gtf file
def get_regions(gene, genes_loc):
    filename = '/storage/cynthiawu/RNAseq/GRCh38/gencode.v31.primary_assembly.annotation.sorted.gtf.gz'
    regions = []
    coordinate = genes_loc[gene]
    location = f'{coordinate[0]}:{coordinate[1]}-{coordinate[2]}'
    tabix = subprocess.Popen(['tabix', filename, location], stdout=subprocess.PIPE, universal_newlines=True)
    for line in tabix.stdout.readlines():
       regions.append(((re.split(';|\t', line.rstrip()))))
    
    regions = pd.DataFrame(regions)
    regions = regions[[0, 2, 3, 4, 9, 14, 15, 17]]
#     types = ['exon', 'UTR']
    types = ['CDS']
    regions = regions.loc[regions[2].isin(types)]
    # transcript_id, exon number, exon_id, protein_id
    regions[9] = regions[9].map(lambda name: name.split('"')[1])
    regions[14] = regions[14].map(lambda name: int(name.split(' ')[2]))
    regions[15] = regions[15].map(lambda name: name.split('"')[1])
    regions[17] = regions[17].map(lambda name: name.split('"')[1])

    return regions

#get the sgRNAs for each region in the gtf file    
def get_sgRNAs_for_region(BE, regions):
    positions = []
    info = []
    transcripts_len = calculate_exon_len(regions)
    for index, row in regions.iterrows():
        transcript = row[9]
        exon = row[14]
        length_5end =transcripts_len[transcript][exon]
        total_trans_len = transcripts_len[transcript]['total']
        info = info + find_sgRNA(list(row), BE , BE_library, length_5end, total_trans_len)
    info = pd.DataFrame(info)
    info = info[[0,1,4,5,6,7,8,9,10,11,12,13,14,15]]
    info = info.rename(columns={0:'chr', 1:'type', 4:'transcript_id', 5:'exon_number', 6:'exon_id', 7:'protein_id', 8:'position', 9:'sgRNA+PAM', 10:'seq', 11:'#edits', 12:'edit_position', 13:'strand', 14:'5end_distance', 15:'3end_distance'})
    return info
	
#dictionary (transcript_id) of dictionary (exon #) containing the length of the transcript up to the exon, also contain total transcript length
def calculate_exon_len(regions):
    regions = regions.loc[regions[2] == 'CDS'].sort_values(by=[9, 14])
    transcripts = {}
    for index, row in regions.iterrows():
        transcript_id = row[9]
        exon = int(row[14])
        length = int(row[4]) - int(row[3])
        if transcript_id in transcripts:
            transcripts[transcript_id][exon] = transcripts[transcript_id]['total']
            transcripts[transcript_id]['total'] += length
        else:
            exons = {}
            exons[exon] = 0
            transcripts[transcript_id] = exons
            transcripts[transcript_id]['total'] = length
    return transcripts

def run_BE_tool(gene, BE, output):
    with open('/storage/cynthiawu/BE_tool/genes_loc.pkl', 'rb') as f:
        genes_loc = pickle.load(f)
    regions = get_regions(gene, genes_loc)
    regions = get_sgRNAs_for_region(BE, regions)
    outfile = f'{output}{gene}_{BE}.csv'
    regions.to_csv(outfile)
    
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-g', '--gene', type=str, dest='gene', help='Gene Name')
    parser.add_argument('-b', '--BE', type=str, dest='BE', help='Base Editor')
    parser.add_argument('-o', '--output', type=str, dest='output', help='Destination folder for output')
    input_values = parser.parse_args()

    run_BE_tool(input_values.gene, input_values.BE, input_values.output)

if __name__ == '__main__':
    main()

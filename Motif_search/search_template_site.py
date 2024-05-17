
import RNA
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import sys
import re
from optparse import OptionParser
from copy import deepcopy


columns = ['seqname', 'start', 'end', 'score', 'strand']

colors_dict = {
    'RED' : '\033[91m',
    'GREEN' : '\033[92m',
    'BLUE' : '\033[94m',
    'CYAN' : '\033[96m',
    'END' : '\033[0m',
}
strand_dict = {
        'forward': '+',
        'reverse': '-',
        '+': '+',
        '-': '-',
        }
tata_distance = 50
tata_motif = 'TATA[T/A]A[T/A]'
min_dist = 10


def load_candidate_anchors(bedfile):
    bed = pd.read_csv(bedfile, sep='\t', header=None, names=columns, index_col=False)
    return bed


def load_genome_handle(filename):
    handle_ = SeqIO.parse(open(filename, 'r'), format='fasta')
    return handle_


def load_repeat_region(filename):
    handle_ = SeqIO.parse(open(filename, 'r'), format='fasta')
    repeat_seqs = [str(rx.seq.reverse_complement().replace('-', '')) for rx in handle_]
    return repeat_seqs


def generate_cyclic_permutations(seq):
    n = len(seq)
    cycles = [str().join([seq[i - j] for i in range(n)]) for j in range(n)]
    return cycles


def write_bed(records, filename):
    names, starts, ends, IDs, scores, strands = [], [], [], [], [], []
    with open(filename, 'w') as f:
        for record in records:
            id_ = record.id.split('_')
            names.append('%s_%s' % (id_[0], id_[1]))
            starts.append(id_[2])
            ends.append(id_[3])
            strands.append(strand_dict.get(id_[4], '.'))
            scores.append('.')
            # f.write('%s\t%s\t%s\t%s\t%s\t%s \n' % (name, start, end, ID, score, strand_dict.get(strand, '.')))
    
    names = ['.'.join(x.split('.')[1:]) for x in names]
    data = {
            'sequence': names,
            'start': map(int, starts),
            'end': map(int, ends),
            'score': scores,
            'strand': strands,
            }
    df = pd.DataFrame(data=data)
    df.sort_values(by=['sequence', 'start'], inplace=True)
    df.drop_duplicates(subset=['sequence', 'start', 'end'], inplace=True)

    IDs = ['TR_%s' % n for n in range(1, len(df) +1)]
    df['ID'] = IDs

    df[['sequence', 'start', 'end', 'ID', 'score', 'strand']].to_csv(filename, index=False, header=None, sep='\t')


def extract_sequences(handle_, candidates_bed, overhang=0, prefix='', search_tata=False):
    candidate_records = []
    relevant_sequences = candidates_bed['seqname'].unique()
    upstream_overhang = overhang
    downstream_overhang = overhang

    for sequence in handle_:
        if prefix == '':
            name = sequence.id
        else:
            name = '%s.%s' % (prefix, sequence.id)
        if name not in relevant_sequences:
            continue
        bed_ = candidates_bed[candidates_bed['seqname'] == name]
        
        for i in range(0, len(bed_)):
            start = bed_['start'].iloc[i] -upstream_overhang
            end = bed_['end'].iloc[i] +downstream_overhang
            strand = strand_dict.get(bed_['strand'].iloc[i], '+')
            this_seq = sequence.seq[start:end]

            record = SeqRecord(
                    this_seq.upper(),
                    id='%s_%s_%s_forward' % (name, start, end),
                    name='%s_%s_%s_forward' % (name, start, end),
                    description='',
                )
            candidate_records.append(record)
            record = SeqRecord(
                    this_seq.reverse_complement().upper(),
                    id='%s_%s_%s_reverse' % (name, start, end),
                    name='%s_%s_%s_reverse' % (name, start, end),
                    description='',
                )
            candidate_records.append(record)
    
    print('Total candidate sequences: %s' % len(candidate_records))
    return candidate_records
    

def find_hairpin_bounds(seq):
    i = 0
    j = len(seq) -1
    struc_start = -1
    struc_end = -1

    while (j > i):
        if (seq[i] == '('):
            struc_start = i
        else:
            i += 1
        if seq[j] == ')':
            struc_end = j
        else:
            j -= 1
        if (struc_start >= 0) and (struc_end > 0):
            break

    return struc_start, struc_end


def search_tata_box(seq):
    match = re.search(tata_motif, str(seq))
    
    if match:
        start, end = match.span()
        return seq[start:end], start, end
    else:
        return '', 0, 0


def assemble_hits(motif_sites, bed):
    records = []
    print(motif_sites)
    data = {
            'sequence' : [x[0].strip('Caenorhabditis_elegans.').split('_')[0] for x in motif_sites],
            'start': [x[1] for x in motif_sites],
            'end': [x[2] for x in motif_sites],
            'empty': ['.' for n in range(len(motif_sites))],
            'strand' : [strand_dict.get(x[0].split('_')[-1]) for x in motif_sites],
        }
    motif_df = pd.DataFrame(data)
    motif_df.drop_duplicates(subset=['sequence', 'start', 'end'], inplace=True)
    motif_df['ID'] = ['TR_{}'.format(n+1) for n in range(len(motif_df))]
    motif_df[['sequence', 'start', 'end', 'ID', 'empty', 'strand']].to_csv(bed, sep='\t', index=False, header=None)
    print(motif_df)


def search_template(record, permutations, length, min_length=0, search_tata=False):
    if (search_tata) == True:
        i = tata_distance
    else:
        i = 0
    j = i+length
    seq = record.seq 
    #
    # TODO: Fix/Replace this try-catch block!!! 
    # 
    try:
        start, end = int(str(record.id).split('_')[2]), int(str(record.id).split('_')[3])
        seq_name = str(record.id).split('_')[1]
    except Exception:
        start, end = int(str(record.id).split('_')[1]), int(str(record.id).split('_')[2])
        seq_name = str(record.id).split('_')[0]
    template_site = ''
    seq_start, seq_end = 0, 0
    motif_hit = False

    while j < len(seq):
        subseq = str(seq[i:j])
        if subseq in permutations:
            if (len(template_site) == 0):
                template_site = subseq
                seq_start = i
                seq_end = j
            else:
                template_site += subseq[-1]
                seq_end = j
        elif (len(template_site) > 0):
            if (len(template_site) < min_length) or (len(template_site) > 14):
                template_site = ''
                continue

            elif (search_tata == True):
                promoter_candidate_region = str(seq[:i])
                box_, tata_start, tata_end = search_tata_box(promoter_candidate_region)
                if (tata_end > 0) and (tata_end +min_dist < seq_start):
                    hairpin_seq = str(seq[tata_end+min_dist+3:seq_start])
                    if len(hairpin_seq) < 5:
                        struc, mfe = '', 0
                    else:
                        struc, mfe = RNA.fold_compound(hairpin_seq).mfe()
                    struc_start, struc_end = find_hairpin_bounds(struc)
                    
                    if (struc_end == -1):
                        i += 1
                        j += 1
                        continue

                    print("\nFound in record: %s " % str().join(record.id.split('_')[:2]))
                    print("TATA box-like motif ", box_, " at position %s - %s." % (start+tata_start, start+tata_end))
                    print("Suspected start of transcription around ~%s." % (start+tata_end+min_dist))
                    print("Predicted structure %s at %s - %s." % (struc[struc_start:struc_end], start+tata_end+min_dist+3+struc_start, start+tata_end+min_dist+3+struc_end))
                    print("Plausible template motif ", template_site, " at position %s - %s. \n" % (start+seq_start, start+seq_end))
                    print('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' %
                        (seq[:tata_start],
                        colors_dict.get('GREEN', '\033[92m'),
                        seq[tata_start:tata_end],
                        colors_dict.get('END', '\033[0m'),
                        seq[tata_end:tata_end+min_dist],
                        colors_dict.get('BLUE'),
                        seq[tata_end+min_dist:tata_end+min_dist+3],
                        colors_dict.get('END'),
                        seq[tata_end+min_dist+3:tata_end+min_dist+3+struc_start],
                        colors_dict.get('CYAN'),
                        seq[tata_end+min_dist+3+struc_start:tata_end+min_dist+3+struc_end],
                        colors_dict.get('END'),
                        seq[tata_end+min_dist+3+struc_end:seq_start],
                        colors_dict.get('RED', '\033[91m'),
                        seq[seq_start:seq_end],
                        colors_dict.get('END', '\033[0m'),
                        seq[seq_end:]
                        )
                        )
                    motif_hit = True
                    template_site = ''
                    gen_start_index = start + tata_end + min_dist
                    gen_end_index   = start + tata_end + min_dist + 400
                    return (str(record.id), gen_start_index, gen_end_index)
            else:
                print("\nFound in record: %s " % record.id)
                gen_start_index = start + seq_start - 200 
                gen_end_index   = start + seq_end + 200
                return (str(record.id), gen_start_index, gen_end_index)

        i += 1
        j += 1

    return []

    '''
                print(template_site, " at position %s - %s." % (start+seq_start, start+seq_end))
                print('%s%s%s%s%s' % 
                    (seq[:seq_start],
                     colors_dict.get('RED', '\033[91m'),
                     seq[seq_start:seq_end],
                     colors_dict.get('END', '\033[0m'),
                     seq[seq_end:]
                    )
                    )
                motif_hit = True
                template_site = ''
    
    
        i += 1
        j += 1

    
    if (motif_hit == True):
        return_record = deepcopy(record)
        if search_tata == False:
            tata_end = 0
            min_dist = 0
            
        else:
            return_record.seq = record.seq[tata_end+min_dist:]
        id_split = str(record.id).split('_')
        if (len(id_split) < 5):
            return_record.id = '_'.join([id_split[0], id_split[1], str(start+tata_end+min_dist), id_split[2], id_split[3]])
        else:
            return_record.id = '_'.join([id_split[0], id_split[1], str(start+tata_end+min_dist), id_split[3], id_split[4]])
        return [return_record]

    return []
    '''


def main():
    parser = OptionParser()
    parser.add_option('-b', '--bed', dest='bed', default='', help='BED file with candidate loci.')
    parser.add_option('-r', '--range', dest='range', type='int', default=0, help='Search range around candidate loci (both directions, Default: 0).')
    parser.add_option('-t', '--tata-box', dest='tata', action='store_true', default=False, help='If this flag is set, will also try to find a TATA-box upstream.')
    parser.add_option('-g', '--genome', dest='genome', type='string', default='', help='Genome file.')
    parser.add_option('-m', '--motif', dest='motif', type='string', default='', help='Fasta file with Telomer repeat motif.')
    parser.add_option('-p', '--prefix', dest='prefix', type='string', default='', help='Name prefix for sequence search.')
    parser.add_option('-l', '--min-length', dest='length', type='int', default=0, help='Report only sites that are at least this many nucleotides long (Default: Length of motif).')
    # parser.add_option('-F', '--fasta-output', dest='fasta', type='string', default='', help='Write sequences producing matches to this FASTA file.')
    parser.add_option('-B', '--bed-output', dest='bed_out', type='string', default='', help='Write coordinates of sequences producing matches to this BED annotation file.')
    options, args = parser.parse_args()

    bed = load_candidate_anchors(options.bed)
    telomer_motif = load_repeat_region(options.motif)[0]
    permutations = generate_cyclic_permutations(telomer_motif)
    genome_handle = load_genome_handle(options.genome)
    candidate_records = extract_sequences(genome_handle, bed, overhang=options.range, prefix=options.prefix, search_tata=options.tata)

    print('\nSet of cyclic permutations: ', permutations)
    length = len(telomer_motif)
    motif_sites = []

    for record in candidate_records:
        search_result = search_template(record, permutations, length, min_length=options.length, search_tata=options.tata)
        if len(search_result):
            motif_sites.append(search_result)

    print('\n\nTotal sites discovered: %s' % len(motif_sites))

    motif_records = assemble_hits(motif_sites, options.bed_out)
    

    '''
    if options.fasta != '':
        SeqIO.write(motif_sites, open(options.fasta, 'w'), format='fasta')
    
    if options.bed_out != '':
        write_bed(options.bed_out)
    '''

main()


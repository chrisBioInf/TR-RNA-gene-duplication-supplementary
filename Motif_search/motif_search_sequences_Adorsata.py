

import sys
import pandas as pd
from Bio import SeqIO
from optparse import OptionParser


colors_dict = {
    'RED' : '\033[91m',
    'GREEN' : '\033[92m',
    'BLUE' : '\033[94m',
    'CYAN' : '\033[96m',
    'END' : '\033[0m',
}

chr_dict = {"BX284601.5" : "chrI",
            "BX284602.5" : "chrII",
            "BX284603.4" : "chrIII",
            "BX284604.4" : "chrIV",
            "BX284605.5" : "chrV",
            "BX284606.5" : "chrX"
            }

columns = ['name', 'source', 'type', 'start', 'end', '1000', 'strand', 'empty', 'attribute']


def parse_record_id(string):
    s = string.split(':')
    name_ = s[0]
    cords = s[1].split('(')
    strand_ = cords[1].strip(')')
    start_, end_ = cords[0].split('-')

    return name_, int(start_), int(end_), strand_


def get_coverage(attribute):
    cov = int(float(attribute.strip('"').split(';')[2].split()[1].strip('"')))
    return cov


def load_repeat_region(filename):
    handle_ = SeqIO.parse(open(filename, 'r'), format='fasta')
    repeat_seqs = [str(rx.seq.reverse_complement().replace('-', '')) for rx in handle_]
    return repeat_seqs


def generate_cyclic_permutations(seq):
    n = len(seq)
    cycles = [str().join([seq[i - j] for i in range(n)]) for j in range(n)]
    return cycles


def main():
    parser = OptionParser()
    parser.add_option('-m', '--motif', dest='motif', type='string', help='Fasta file with Telomer repeat motif.')
    parser.add_option('-f', '--file', dest='input', type='string', help='Fasta file with candidate sequences.')
    parser.add_option('-r', '--reference', dest='reference', type='string', default='', help='Reference GTF to pull coverage and coordinates from.')
    parser.add_option('-l', '--min-length', dest='length', type='int', default=0, help='Report only sites that are at least this many nucleotides long (Default: 1,0 x Length of motif).')
    # parser.add_option('-B', '--bed-output', dest='bed_out', type='string', default='', help='Write coordinates of sequences producing matches to this BED annotation file.')
    options, args = parser.parse_args()

    if not options.input:
        parser.error('No file with candidate sequences given (--file).')
    if not options.motif:
        parser.error('No file with telomere repeat motifs given (--motif).')

    motifs = load_repeat_region(options.motif)
    df_ = pd.read_csv(options.reference, sep='\t', names=columns)
    print(df_)
    
    TR_chrs = []
    TR_starts = []
    TR_ends = []
    TR_names = []
    TR_scores = []
    TR_strands = []
    TR_templates = []
    
    TR_count = 0
    Motif_count = 0

    for motif in motifs:
        cycles = generate_cyclic_permutations(motif)
        print('Candidate motifs: %s \n' % ' '.join(cycles))
        handle = SeqIO.parse(open(options.input, 'r'), format='fasta')
        
        if (options.length == 0):
            min_length = 1.0 * len(motif)
        else:
            min_length = options.length

        for record in handle:
            seq = str(record.seq).upper()
            i, j = 0, len(motif)
            template_site = ''
            seq_start, seq_end = 0, 0
            motif_hit = False

            while j < len(seq):
                subseq = str(seq[i:j])
                if subseq in cycles:
                    if (len(template_site) == 0):
                        template_site = subseq
                        seq_start = i
                        seq_end = j
                    else:
                        template_site += subseq[-1]
                        seq_end = j
                elif (len(template_site) > 0):
                    if (len(template_site) >= min_length) and (len(template_site) <= 20):
                        Motif_count += 1
                        print(seq_start, seq_end, template_site)
                        print(record.id)
                        print('%s%s%s%s%s' % (seq[:seq_start], colors_dict.get('RED'), seq[seq_start:seq_end], colors_dict.get('END'), seq[seq_end:]))
                        chr_, start_, end_, strand_ = parse_record_id(record.id)
                        ref_data = df_[(df_['name'] == chr_ ) & (df_['start'] == int(start_)+1) & (df_['end'] == int(end_))]
                        cov_ = get_coverage(ref_data['attribute'].iloc[0])
                        chr_ = chr_dict.get(chr_, chr_)

                        if (start_ not in TR_starts) and (end_ not in TR_ends): 
                            TR_count += 1
                            TR_chrs.append(chr_)
                            TR_starts.append(start_)
                            TR_ends.append(end_)
                            TR_names.append('TR_%s' % TR_count)
                            TR_scores.append(cov_)
                            TR_strands.append(strand_)
                            TR_templates.append(template_site)

                        TR_chrs.append(chr_)
                        TR_starts.append(start_ + seq_start)
                        TR_ends.append(start_ + seq_start + len(template_site))
                        TR_names.append('Motif_%s' % Motif_count)
                        TR_scores.append(cov_)
                        TR_strands.append(strand_)
                        TR_templates.append(template_site)

                    template_site = ''
                else:
                    pass

                i += 1
                j += 1

    df = pd.DataFrame(data=
                        {'sequence': TR_chrs,
                         'start': TR_starts,
                         'end': TR_ends,
                         'id': TR_names,
                         'score': TR_scores,
                         'strand': TR_strands,
                         'template': TR_templates,
                         }
                      )
    df.sort_values(by=['sequence', 'start'], ascending=True, inplace=True)
    
    df = df[df['id'].str.startswith('TR')]
    
    print(df)
    df.to_csv('C_elegans_Candidates.bed', sep='\t', index=False, header=None)

main()


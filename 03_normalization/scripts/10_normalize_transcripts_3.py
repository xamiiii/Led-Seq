#   purpose of this scipt: take the raw counts (transcriptome wide) from negative and positive strand and convert to
#   normalized counts ONLY in non-overlapping regions of transcripts. positions within those transcripts that are
#   without signal, are assigned zero before normalization.
# ______________________________________________________________________________________________________________________
#   usage: python script raw_counts_posStrand raw_counts_negStrand in:gff_like  out:output_table
#
#   input is: 2  3-column-raw_count files: 1) signal == raw count 2) sequence name (e.g. chromosome id) 3) position in sequence (0-BASED!)
#   gff-(like)-file with NON-OVERLAPPING transcript regions/genes: 1) sequence name 4) start 5) end  7) strand 9) name/id/etc  (1-BASED!)
#
#   output is: 5-column-file: 1) sequence name 2)position 3)normed count (n.p. if not possible to normalize 4)transcript_name
#   5)strand ("+"/"-")                                                                                        (1-BASED!)
# ______________________________________________________________________________________________________________________

import numpy as np
from Bio import SeqIO

rc_pos_file = snakemake.input['pos_strand']
rc_neg_file = snakemake.input['neg_strand']
input = snakemake.input['no_overlaps']
genome = snakemake.input['genome']
output = snakemake.output['normalize_counts']

infile_rc_pos = open(rc_pos_file, 'r')
infile_rc_neg = open(rc_neg_file, 'r')
outfile = open(output, 'w')
infile_transcripts = open(input, 'r')

"""open up sequence dict in memory"""
genome_dict = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))

def find_sequence_in_genome(seqid, strand, pos):
    base = genome_dict[seqid].seq[pos-1:pos]   # make it 0-based
    if strand == "-":
        base = base.complement()
    base = str(base.transcribe())
    return base

rc_pos = {}
rc_neg = {}

for line in infile_rc_pos:
    columns = line.rstrip().split()
    rc = int(columns[0])
    seq_id = columns[1]
    position = int(columns[2])+1    # make it 1-based

    if seq_id in rc_pos:
        rc_pos[seq_id][position] = rc
    else:
        rc_pos[seq_id] = dict()
        rc_pos[seq_id][position] = rc

for line in infile_rc_neg:
    columns = line.rstrip().split()
    rc = int(columns[0])
    seq_id = columns[1]
    position = int(columns[2])+1    # make it 1-based

    if seq_id in rc_neg:
        rc_neg[seq_id][position] = rc
    else:
        rc_neg[seq_id] = dict()
        rc_neg[seq_id][position] = rc

def normalize(signal_array):
    normalized_counts = []

    q98 = np.percentile(signal_array, 98)
    q90 = np.percentile(signal_array, 90)
    # print(q90, q98)

    high_counts = []
    for e in signal_array:
        if q90 <= e <= q98:
            high_counts.append(e)
    # print(high_counts)
    if len(high_counts) > 0:
        norm_factor = np.mean(high_counts)

        if norm_factor == 0.0:
            return False
        else:
            for e in signal_array:
                normed_count = round((float(e) / norm_factor), 2)

                # capping
                if normed_count > 7:
                    normed_count = 7.00

                normalized_counts.append(normed_count)

            # print(normalized_counts)
            return normalized_counts
    else:
        return False

# transcript_lengths = []
# counter = 0
for line in infile_transcripts:
    # counter += 1
    # if counter >= 2:
    #     break
    columns = line.rstrip().split()
    seq_id = columns[0]
    start = int(columns[3])
    end = int(columns[4])
    strand = columns[6]
    name = columns[8]
    # generate array with positions for this transcript
    positions = []
    # generate array with raw counts, if no signal, default is 0
    raw_counts = []

    bases = []
    neighboring_bases = []

    if strand == "+":
        for i in range(start, end+1):
            base = find_sequence_in_genome(seq_id, strand, i)
            neighboring_base = ""
            if i == end:
                neighboring_base = "N"
            else:
                neighboring_base = find_sequence_in_genome(seq_id, strand, i + 1)
            bases.append(base)
            neighboring_bases.append(neighboring_base)
            positions.append(i)
            if seq_id in rc_pos:
                if i in rc_pos[seq_id]:
                    raw_counts.append(rc_pos[seq_id][i])
                else:
                    raw_counts.append(0)
            else:
                raw_counts.append(0)
    elif strand == "-":
        for i in range(start, end+1):
            base = find_sequence_in_genome(seq_id, strand, i)
            neighboring_base = ""
            if i == start:
                neighboring_base = "N"
            else:
                neighboring_base = find_sequence_in_genome(seq_id, strand, i - 1)
            bases.append(base)
            neighboring_bases.append(neighboring_base)
            positions.append(i)
            if seq_id in rc_neg:
                if i in rc_neg[seq_id]:
                    raw_counts.append(rc_neg[seq_id][i])
                else:
                    raw_counts.append(0)
            else:
                raw_counts.append(0)
    # print(raw_counts)

    normalized_counts = normalize(raw_counts)
    if normalized_counts != False:          #normalization was possible
        for i in range(len(positions)):
            outfile.write(seq_id + "\t" + strand + "\t" + str(positions[i]) + "\t"+ bases[i]+ "\t"+ neighboring_bases[i] + "\t" + name + "\t" + str(raw_counts[i]) + "\t" + str(normalized_counts[i]) + "\n")
    else:                                   # normalization not possible, not enough data
        for i in range(len(positions)):
            outfile.write(seq_id + "\t" + strand + "\t" + str(positions[i]) + "\t"+ bases[i]+ "\t"+ neighboring_bases[i] + "\t" + name + "\t" + str(raw_counts[i]) + "\t" + "n.p." + "\n")

outfile.close()
infile_transcripts.close()
infile_rc_pos.close()
infile_rc_neg.close()
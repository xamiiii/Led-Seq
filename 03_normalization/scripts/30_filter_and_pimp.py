import os

ref_structures_file = snakemake.config['ref_structures']
covered_transcripts_file = snakemake.input['covered_transcripts_gff']
library_file = snakemake.input['normalize_counts']
output = snakemake.output['filtered_annotated_data']

library_name = os.path.basename(library_file).split(".")[0]
if "OH" in library_name:
    library_type = "OH"
else:
    library_type = "cP"


"""store the reference structures and sequences in dict"""
ref_structures = {}
line_counter=0
current_id = ""

with open(ref_structures_file, 'r') as file:
    for line in file:
        line_counter += 1
        line=line.rstrip()
        if line_counter == 1:           #transcript_id
            if line.startswith(">"):
                current_id=line[1:]
            else: print ("error")
        if line_counter == 2:           # sequence
           ref_structures[current_id]=[line]
        if line_counter == 3:           # structure
           ref_structures[current_id].append(line)
           line_counter = 0
           current_id = ""

"""read in the relevant transcripts (that are sufficiently covered)"""
covered_transcripts = {}

with open(covered_transcripts_file, 'r') as file:
    for line in file:
        columns = line.rstrip().split()
        seqid, start, end, strand, transcript_id = columns[0], int(columns[3]), int(columns[4]), columns[6], columns[8]
        covered_transcripts[transcript_id]=[seqid, start, end, strand]


def transform_postion(position, start, end, strand):
    if strand == "+":
        i = position - start + 1  # convert pos to relative position in sequence
    elif strand == "-":
        i = end - position + 1  # convert pos to relative position in sequence
    return i


def extract_structure(transcript_id, position):
    char = ref_structures[transcript_id][1][position-1]
    if char == "(" or char == ")":
        return 1
    elif char == ".":
        return 0


"""iterate over library data and write out filtered and pimped information"""
with open(output, 'w') as output_file:
    with open(library_file, 'r') as file:
        for line in file:
            columns = line.rstrip().split()
            transcript_id = columns[5]
            if transcript_id in covered_transcripts:
                seqid, strand, position, base, neighboring_base, raw_count, normalized_count = columns[0], columns[1], int(columns[2]), columns[3], columns[4], columns[6], columns[7]
                transcript_start = covered_transcripts[transcript_id][1]
                transcript_end = covered_transcripts[transcript_id][2]
                i=transform_postion(position, transcript_start, transcript_end, strand)
                if library_type == "cP":
                    if i <= 11:
                        normalized_count = "NA"
                elif library_type == "OH":
                    transcript_length = transcript_end - transcript_start + 1
                    if i > (transcript_length - 12):
                        normalized_count = "NA"

                if transcript_id in ref_structures:
                    structure = extract_structure(transcript_id, i)
                    # sequence = ref_structures[transcript_id][0][i-1]
                else:
                    structure = "-"
                    # sequence = "-"

                output_file.write(seqid + "\t" + strand + "\t" + str(position)+ "\t"+ base + "\t"+ neighboring_base + "\t" + transcript_id + "\t" + raw_count + "\t" + normalized_count + "\t" + str(i)+ "\t" + str(structure) + "\t" + library_name + "\t" + library_type + "\n")

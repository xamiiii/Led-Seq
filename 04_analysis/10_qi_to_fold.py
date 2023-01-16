import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
import sys
import RNA


signal_file = sys.argv[1]
ref_structures_file = sys.argv[2]
temperature = int(sys.argv[3])
output_folder = sys.argv[4]

"""create and open output files"""
mfe_structures_no_sc= open(output_folder + "/benchmark_mfe_structures_no_sc.fa", 'w')

mfe_structures_with_sc_pb= open(output_folder + "/benchmark_mfe_structures_with_sc_pb.fa", 'w')
mfe_structures_with_sc_h2o= open(output_folder + "/benchmark_mfe_structures_with_sc_h2o.fa", 'w')
mfe_structures_with_sc_both= open(output_folder + "/benchmark_mfe_structures_with_sc_both.fa", 'w')

cP_mfe_structures_with_sc_pb= open(output_folder + "/cP_benchmark_mfe_structures_with_sc_pb.fa", 'w')
cP_mfe_structures_with_sc_h2o= open(output_folder + "/cP_benchmark_mfe_structures_with_sc_h2o.fa", 'w')
cP_mfe_structures_with_sc_both= open(output_folder + "/cP_benchmark_mfe_structures_with_sc_both.fa", 'w')

OH_mfe_structures_with_sc_pb= open(output_folder + "/OH_benchmark_mfe_structures_with_sc_pb.fa", 'w')
OH_mfe_structures_with_sc_h2o= open(output_folder + "/OH_benchmark_mfe_structures_with_sc_h2o.fa", 'w')
OH_mfe_structures_with_sc_both= open(output_folder + "/OH_benchmark_mfe_structures_with_sc_both.fa", 'w')

centroid_structures_no_sc= open(output_folder + "/benchmark_centroid_structures_no_sc.fa", 'w')

centroid_structures_with_sc_pb= open(output_folder + "/benchmark_centroid_structures_with_sc_pb.fa", 'w')
centroid_structures_with_sc_h2o= open(output_folder + "/benchmark_centroid_structures_with_sc_h2o.fa", 'w')
centroid_structures_with_sc_both= open(output_folder + "/benchmark_centroid_structures_with_sc_both.fa", 'w')

cP_centroid_structures_with_sc_pb= open(output_folder + "/cP_benchmark_centroid_structures_with_sc_pb.fa", 'w')
cP_centroid_structures_with_sc_h2o= open(output_folder + "/cP_benchmark_centroid_structures_with_sc_h2o.fa", 'w')
cP_centroid_structures_with_sc_both= open(output_folder + "/cP_benchmark_centroid_structures_with_sc_both.fa", 'w')

OH_centroid_structures_with_sc_pb= open(output_folder + "/OH_benchmark_centroid_structures_with_sc_pb.fa", 'w')
OH_centroid_structures_with_sc_h2o= open(output_folder + "/OH_benchmark_centroid_structures_with_sc_h2o.fa", 'w')
OH_centroid_structures_with_sc_both= open(output_folder + "/OH_benchmark_centroid_structures_with_sc_both.fa", 'w')


soft_constraints_output_pb= open(output_folder + "/benchmark_soft_constraints_pb.txt", 'w')
soft_constraints_output_h2o= open(output_folder + "/benchmark_soft_constraints_h2o.txt", 'w')
soft_constraints_output_both= open(output_folder + "/benchmark_soft_constraints_both.txt", 'w')

cP_soft_constraints_output_pb= open(output_folder + "/cP_benchmark_soft_constraints_pb.txt", 'w')
cP_soft_constraints_output_h2o= open(output_folder + "/cP_benchmark_soft_constraints_h2o.txt", 'w')
cP_soft_constraints_output_both= open(output_folder + "/cP_benchmark_soft_constraints_both.txt", 'w')

OH_soft_constraints_output_pb= open(output_folder + "/OH_benchmark_soft_constraints_pb.txt", 'w')
OH_soft_constraints_output_h2o= open(output_folder + "/OH_benchmark_soft_constraints_h2o.txt", 'w')
OH_soft_constraints_output_both= open(output_folder + "/OH_benchmark_soft_constraints_both.txt", 'w')


"""store the reference structures in dict"""
data = {}
line_counter=0
current_id = ""
pos_unpaired_counter = 0
pos_counter = 0

with open(ref_structures_file, 'r') as file:
    for line in file:
        line_counter += 1
        line=line.rstrip()
        if line.startswith(">"):        #transcript_id
            line_counter = 1
            current_id=line[1:]
        if line_counter == 2:           # sequence [0]
            data[current_id] = {}
            for i, char in enumerate(line):
                data[current_id][i+1]=[char]   #1-based position in sequence

        if line_counter == 3:           # structure [1]
            for i, char in enumerate(line):
                pos_counter += 1
                if char == ".":
                    pos_unpaired_counter += 1
                data[current_id][i + 1].append(char)

    line_counter = 0
    current_id = ""

p0 = pos_unpaired_counter/pos_counter


with open(signal_file, 'r') as file:
    for line in file:
        line=line.rstrip()
        columns = line.split("\t")
        transcript_id, position, qi_cP, qi_OH, qi_both, qi_h2o_cP, qi_h2o_OH, qi_h2o_both = columns[0], int(columns[1]), columns[8], columns[9], columns[10], columns[11], columns[12], columns[13]
        if qi_cP != "-": qi_cP=float(qi_cP)
        if qi_OH != "-": qi_OH=float(qi_OH)
        if qi_both != "-": qi_both=float(qi_both)
        if qi_h2o_cP != "-": qi_h2o_cP=float(qi_h2o_cP)
        if qi_h2o_OH != "-": qi_h2o_OH=float(qi_h2o_OH)
        if qi_h2o_both != "-": qi_h2o_both=float(qi_h2o_both)

        data[transcript_id][position].append(qi_cP)         # qi_cP [2]
        data[transcript_id][position].append(qi_OH)         # qi_OH [3]
        data[transcript_id][position].append(qi_both)       # qi_both [4]
        data[transcript_id][position].append(qi_h2o_cP)     # qi_h2o_cP [5]
        data[transcript_id][position].append(qi_h2o_OH)         # qi_h2o_OH [6]
        data[transcript_id][position].append(qi_h2o_both)       # qi_h2o_both [7]

print("finished reading data")

# def calc_pi(seq, temp):
#     md = RNA.md()
#     md.temperature = temp

#     fc = RNA.fold_compound(seq, md)
#     (mfe_struct, mfe) = fc.mfe()

#     # rescale Boltzmann factors for partition function computation
#     fc.exp_params_rescale(mfe)

#     # compute partition function
#     (propensity, ensemble_energy) = fc.pf()
#     basepair_probs = fc.bpp()

#     up_probabilities = []
#     for i in range(1, len(seq)+1):
#         p_paired = 0
#         for j in range(1, len(seq) + 1):
#             if i<j:
#                 p_paired += basepair_probs[i][j]
#             elif i>j:
#                 p_paired += basepair_probs[j][i]

#         p_unpaired = 1- p_paired
#         up_probabilities.append(p_unpaired)

#     return up_probabilities

# print("calculating pi...")
# for transcript_id in data:
#     if transcript_id in ["16S_ribosomal_RNA", "23S_ribosomal_RNA"]:
#         continue
#     sequence=""
#     for pos in sorted(data[transcript_id]):
#         sequence+=data[transcript_id][pos][0]
#     pi_predicted = calc_pi(sequence, temperature)
#     for i, pi in enumerate(pi_predicted):
#         data[transcript_id][i+1].append(pi)         # pi [8]

#     print(transcript_id + " ready")

# print("finished calculating pi")


def calcSoftConstraints(kT, qi):#, pi):
    alpha = 1.2
    vi = -kT * alpha * max((np.log(qi / (1 - qi)) - np.log(p0 / (1 - p0))),0)   #max(ln.., 0) forbids destabilization/pushing pairdness
    return vi

def calcSoftConstraints_with_h2o(kT, qi, qi_h2o):
    alpha = 0.5
    beta =  0.5

    vi = (-kT * alpha * max((np.log(qi / (1 - qi)) - np.log(p0 / (1 - p0))),0)) +  (-kT * beta * max((np.log(qi_h2o / (1 - qi_h2o)) - np.log(p0 / (1 - p0))),0))  #max(ln.., 0) forbids destabilization/pushing pairdness

    return vi

def setSoftConstraints_with_h2o(transcript_id):
    m_d = RNA.md()
    m_d.temperature = temperature
    params = RNA.exp_param(m_d)
    kT = params.kT / 1000.

    sc = []
    sc_cP = []
    sc_OH = []

    for pos in sorted(data[transcript_id]):
        qi_cP = data[transcript_id][pos][2]
        qi_OH = data[transcript_id][pos][3]
        qi_both = data[transcript_id][pos][4]
        qi_h2o_cP = data[transcript_id][pos][5]
        qi_h2o_OH = data[transcript_id][pos][6]
        qi_h2o_both = data[transcript_id][pos][7]
        #pi = data[transcript_id][pos][8]

        vi = None
        vi_cP = None
        vi_OH = None

        if qi_cP != "-":
            if qi_h2o_cP != "-":
                vi_cP = calcSoftConstraints_with_h2o(kT, qi_cP, qi_h2o_cP)
            else:
                print("h2o samples does not match")

        if qi_OH != "-":
            if qi_h2o_OH != "-":
                vi_OH = calcSoftConstraints_with_h2o(kT, qi_OH, qi_h2o_OH)
            else:
                print["h2o samples does not match"]

        if qi_both != "-":
            if qi_h2o_both != "-":
                vi = calcSoftConstraints_with_h2o(kT, qi_both, qi_h2o_both)
            else:
                print["h2o samples does not match"]

        if vi is not None:
            sc.append(vi)
        else:
            if vi_cP is not None:
                sc.append(vi_cP)
            elif vi_OH is not None:
                sc.append(vi_OH)
            else:
                print("error")
                quit()


        if vi_cP is not None:
            sc_cP.append(vi_cP)
        else:
            sc_cP.append(-0.0)

        if vi_OH is not None:
            sc_OH.append(vi_OH)
        else:
            sc_OH.append(-0.0)

    return sc, sc_cP, sc_OH

def setSoftConstraints(transcript_id):
    m_d = RNA.md()
    m_d.temperature = temperature
    params = RNA.exp_param(m_d)
    kT = params.kT / 1000.

    sc_pb = []
    sc_cP_pb = []
    sc_OH_pb = []
    sc_h2o = []
    sc_cP_h2o = []
    sc_OH_h2o = []

    for pos in sorted(data[transcript_id]):
        qi_cP = data[transcript_id][pos][2]
        qi_OH = data[transcript_id][pos][3]
        qi_both = data[transcript_id][pos][4]
        qi_h2o_cP = data[transcript_id][pos][5]
        qi_h2o_OH = data[transcript_id][pos][6]
        qi_h2o_both = data[transcript_id][pos][7]
        #pi = data[transcript_id][pos][8]

        vi_pb = None
        vi_cP_pb = None
        vi_OH_pb = None

        if qi_cP != "-":
            vi_cP_pb = calcSoftConstraints(kT, qi_cP)#, pi)

        if qi_OH != "-":
            vi_OH_pb = calcSoftConstraints(kT, qi_OH)#, pi)

        if qi_both != "-":
            vi_pb = calcSoftConstraints(kT, qi_both)#, pi)

        if vi_pb is not None:
            sc_pb.append(vi_pb)
        else:
            if vi_cP_pb is not None:
                sc_pb.append(vi_cP_pb)
            elif vi_OH_pb is not None:
                sc_pb.append(vi_OH_pb)
            else:
                print("error")
                quit()

        if vi_cP_pb is not None:
            sc_cP_pb.append(vi_cP_pb)
        else:
            sc_cP_pb.append(-0.0)

        if vi_OH_pb is not None:
            sc_OH_pb.append(vi_OH_pb)
        else:
            sc_OH_pb.append(-0.0)

        vi_h2o = None
        vi_cP_h2o = None
        vi_OH_h2o = None

        if qi_h2o_cP != "-":
            vi_cP_h2o = calcSoftConstraints(kT, qi_h2o_cP)#, pi)

        if qi_h2o_OH != "-":
            vi_OH_h2o = calcSoftConstraints(kT, qi_h2o_OH)#, pi)

        if qi_h2o_both != "-":
            vi_h2o = calcSoftConstraints(kT, qi_h2o_both)#, pi)

        if vi_h2o is not None:
            sc_h2o.append(vi_h2o)
        else:
            if vi_cP_h2o is not None:
                sc_h2o.append(vi_cP_h2o)
            elif vi_OH_h2o is not None:
                sc_h2o.append(vi_OH_h2o)
            else:
                print("error")
                quit()

        if vi_cP_h2o is not None:
            sc_cP_h2o.append(vi_cP_h2o)
        else:
            sc_cP_h2o.append(-0.0)

        if vi_OH_h2o is not None:
            sc_OH_h2o.append(vi_OH_h2o)
        else:
            sc_OH_h2o.append(-0.0)

    return sc_pb, sc_cP_pb, sc_OH_pb, sc_h2o, sc_cP_h2o, sc_OH_h2o

def addSoftConstraintsToFc(fc, soft_constraints):

    for i, v in enumerate(soft_constraints):
        fc.sc_add_up(i+1, v)        #position in sequence must be 1-based


def foldMfe(transcript_id, file_structures_no_sc,
            file_structures_with_sc_pb, file_structures_with_sc_h2o, file_structures_with_sc_both,
            file_sc_pb, file_sc_h2o, file_sc_both,
            cP_file_structures_with_sc_pb, cP_file_structures_with_sc_h2o, cP_file_structures_with_sc_both,
            cP_file_sc_pb, cP_file_sc_h2o, cP_file_sc_both,
            OH_file_structures_with_sc_pb, OH_file_structures_with_sc_h2o, OH_file_structures_with_sc_both,
            OH_file_sc_pb, OH_file_sc_h2o, OH_file_sc_both):

    """get sequence"""
    seq=""
    for pos in sorted(data[transcript_id]):
        seq+=data[transcript_id][pos][0]

    """calculate soft_constraints and save in file"""
    soft_constraints_pb, cP_soft_constraints_pb, OH_soft_constraints_pb, soft_constraints_h2o, cP_soft_constraints_h2o, OH_soft_constraints_h2o = setSoftConstraints(transcript_id)
    soft_constraints_both, cP_soft_constraints_both, OH_soft_constraints_both = setSoftConstraints_with_h2o(transcript_id)

    for pos, sc in enumerate(soft_constraints_pb):
        file_sc_pb.write(transcript_id + "\t" +  str(pos + 1) + "\t" + str(sc) + "\n")
    for pos, sc in enumerate(soft_constraints_h2o):
        file_sc_h2o.write(transcript_id + "\t" +  str(pos + 1) + "\t" + str(sc) + "\n")
    for pos, sc in enumerate(soft_constraints_both):
        file_sc_both.write(transcript_id + "\t" +  str(pos + 1) + "\t" + str(sc) + "\n")

    for pos, sc in enumerate(cP_soft_constraints_pb):
        cP_file_sc_pb.write(transcript_id + "\t" +  str(pos + 1) + "\t" + str(sc) + "\n")
    for pos, sc in enumerate(cP_soft_constraints_h2o):
        cP_file_sc_h2o.write(transcript_id + "\t" +  str(pos + 1) + "\t" + str(sc) + "\n")
    for pos, sc in enumerate(cP_soft_constraints_both):
        cP_file_sc_both.write(transcript_id + "\t" +  str(pos + 1) + "\t" + str(sc) + "\n")

    for pos, sc in enumerate(OH_soft_constraints_pb):
        OH_file_sc_pb.write(transcript_id + "\t" +  str(pos + 1) + "\t" + str(sc) + "\n")
    for pos, sc in enumerate(OH_soft_constraints_h2o):
        OH_file_sc_h2o.write(transcript_id + "\t" +  str(pos + 1) + "\t" + str(sc) + "\n")
    for pos, sc in enumerate(OH_soft_constraints_both):
        OH_file_sc_both.write(transcript_id + "\t" +  str(pos + 1) + "\t" + str(sc) + "\n")


    md = RNA.md()
    md.temperature = temperature

    """compute WITHOUT soft constraints"""
    fc = RNA.fold_compound(seq, md)
    (structure_no_sc, mfe_no_sc) = fc.mfe()
    file_structures_no_sc.write(">" + transcript_id + "," + str(mfe_no_sc) + "\n" + seq + "\n" + structure_no_sc + "\n")

    """compute WITH soft constraints (Pb)"""
    fc2 = RNA.fold_compound(seq, md)
    addSoftConstraintsToFc(fc2, soft_constraints_pb)
    (structure_with_sc_pb, mfe_with_sc_pb) = fc2.mfe()
    file_structures_with_sc_pb.write(">" + transcript_id + "," + str(mfe_with_sc_pb) + "\n" + seq + "\n" + structure_with_sc_pb + "\n")

    """compute WITH soft constraints based on cP library only (Pb)"""
    fc3 = RNA.fold_compound(seq, md)
    addSoftConstraintsToFc(fc3, cP_soft_constraints_pb)
    (cP_structure_with_sc_pb, cP_mfe_with_sc_pb) = fc3.mfe()
    cP_file_structures_with_sc_pb.write(">" + transcript_id + "," + str(cP_mfe_with_sc_pb) + "\n" + seq + "\n" + cP_structure_with_sc_pb + "\n")

    """compute WITH soft constraints based on OH library only (Pb)"""

    fc4 = RNA.fold_compound(seq, md)
    addSoftConstraintsToFc(fc4, OH_soft_constraints_pb)
    (OH_structure_with_sc_pb, OH_mfe_with_sc_pb) = fc4.mfe()
    OH_file_structures_with_sc_pb.write(">" + transcript_id + "," + str(OH_mfe_with_sc_pb) + "\n" + seq + "\n" + OH_structure_with_sc_pb + "\n")

    """compute WITH soft constraints (H2O)"""
    fc5 = RNA.fold_compound(seq, md)
    addSoftConstraintsToFc(fc5, soft_constraints_h2o)
    (structure_with_sc_h2o, mfe_with_sc_h2o) = fc5.mfe()
    file_structures_with_sc_h2o.write(
        ">" + transcript_id + "," + str(mfe_with_sc_h2o) + "\n" + seq + "\n" + structure_with_sc_h2o + "\n")

    """compute WITH soft constraints based on cP library only (H2O)"""
    fc6 = RNA.fold_compound(seq, md)
    addSoftConstraintsToFc(fc6, cP_soft_constraints_h2o)
    (cP_structure_with_sc_h2o, cP_mfe_with_sc_h2o) = fc6.mfe()
    cP_file_structures_with_sc_h2o.write(
        ">" + transcript_id + "," + str(cP_mfe_with_sc_h2o) + "\n" + seq + "\n" + cP_structure_with_sc_h2o + "\n")

    """compute WITH soft constraints based on OH library only (H2O)"""

    fc7 = RNA.fold_compound(seq, md)
    addSoftConstraintsToFc(fc7, OH_soft_constraints_h2o)
    (OH_structure_with_sc_h2o, OH_mfe_with_sc_h2o) = fc7.mfe()
    OH_file_structures_with_sc_h2o.write(
        ">" + transcript_id + "," + str(OH_mfe_with_sc_h2o) + "\n" + seq + "\n" + OH_structure_with_sc_h2o + "\n")

    """compute WITH soft constraints (both)"""
    fc8 = RNA.fold_compound(seq, md)
    addSoftConstraintsToFc(fc8, soft_constraints_both)
    (structure_with_sc_both, mfe_with_sc_both) = fc8.mfe()
    file_structures_with_sc_both.write(
        ">" + transcript_id + "," + str(mfe_with_sc_both) + "\n" + seq + "\n" + structure_with_sc_both + "\n")

    """compute WITH soft constraints based on cP library only (both)"""
    fc9 = RNA.fold_compound(seq, md)
    addSoftConstraintsToFc(fc9, cP_soft_constraints_both)
    (cP_structure_with_sc_both, cP_mfe_with_sc_both) = fc9.mfe()
    cP_file_structures_with_sc_both.write(
        ">" + transcript_id + "," + str(cP_mfe_with_sc_both) + "\n" + seq + "\n" + cP_structure_with_sc_both + "\n")

    """compute WITH soft constraints based on OH library only (both)"""

    fc10 = RNA.fold_compound(seq, md)
    addSoftConstraintsToFc(fc10, OH_soft_constraints_both)
    (OH_structure_with_sc_both, OH_mfe_with_sc_both) = fc10.mfe()
    OH_file_structures_with_sc_both.write(
        ">" + transcript_id + "," + str(OH_mfe_with_sc_both) + "\n" + seq + "\n" + OH_structure_with_sc_both + "\n")




    return None


def foldCentroid(transcript_id, file_structures_no_sc,
            file_structures_with_sc_pb, file_structures_with_sc_h2o, file_structures_with_sc_both,
            file_sc_pb, file_sc_h2o, file_sc_both,
            cP_file_structures_with_sc_pb, cP_file_structures_with_sc_h2o, cP_file_structures_with_sc_both,
            cP_file_sc_pb, cP_file_sc_h2o, cP_file_sc_both,
            OH_file_structures_with_sc_pb, OH_file_structures_with_sc_h2o, OH_file_structures_with_sc_both,
            OH_file_sc_pb, OH_file_sc_h2o, OH_file_sc_both):
    """get sequence"""
    seq=""
    for pos in sorted(data[transcript_id]):
        seq+=data[transcript_id][pos][0]

    """calculate soft_constraints and save in file"""
    soft_constraints_pb, cP_soft_constraints_pb, OH_soft_constraints_pb, soft_constraints_h2o, cP_soft_constraints_h2o, OH_soft_constraints_h2o = setSoftConstraints(transcript_id)
    soft_constraints_both, cP_soft_constraints_both, OH_soft_constraints_both = setSoftConstraints_with_h2o(transcript_id)

    md = RNA.md()
    md.temperature = temperature

    """compute WITHOUT soft constraints"""
    fc = RNA.fold_compound(seq, md)
    (structure_no_sc, mfe_no_sc) = fc.mfe()
    fc.exp_params_rescale(mfe_no_sc)
    (pp, pf) = fc.pf()
    (centroid_struct_no_sc, dist_no_sc) = fc.centroid()
    file_structures_no_sc.write(">" + transcript_id + "," + str(dist_no_sc) + "\n" + seq + "\n" + centroid_struct_no_sc + "\n")

    """compute WITH soft constraints (Pb)"""
    fc2 = RNA.fold_compound(seq, md)
    addSoftConstraintsToFc(fc2, soft_constraints_pb)
    (structure_with_sc_pb, mfe_with_sc_pb) = fc2.mfe()
    fc2.exp_params_rescale(mfe_with_sc_pb)
    (pp, pf) = fc2.pf()
    (centroid_structure_with_sc_pb, dist_with_sc_pb) = fc2.centroid()
    file_structures_with_sc_pb.write(
        ">" + transcript_id + "," + str(dist_with_sc_pb) + "\n" + seq + "\n" + centroid_structure_with_sc_pb + "\n")

    """compute WITH soft constraints based on cP library only (Pb)"""
    fc3 = RNA.fold_compound(seq, md)
    addSoftConstraintsToFc(fc3, cP_soft_constraints_pb)
    (cP_structure_with_sc_pb, cP_mfe_with_sc_pb) = fc3.mfe()
    fc3.exp_params_rescale(cP_mfe_with_sc_pb)
    (pp, pf) = fc3.pf()
    (cP_centroid_structure_with_sc_pb, cP_dist_with_sc_pb) = fc3.centroid()
    cP_file_structures_with_sc_pb.write(
        ">" + transcript_id + "," + str(cP_dist_with_sc_pb) + "\n" + seq + "\n" + cP_centroid_structure_with_sc_pb + "\n")

    """compute WITH soft constraints based on OH library only (Pb)"""

    fc4 = RNA.fold_compound(seq, md)
    addSoftConstraintsToFc(fc4, OH_soft_constraints_pb)
    (OH_structure_with_sc_pb, OH_mfe_with_sc_pb) = fc4.mfe()
    fc4.exp_params_rescale(OH_mfe_with_sc_pb)
    (pp, pf) = fc4.pf()
    (OH_centroid_structure_with_sc_pb, OH_dist_with_sc_pb) = fc4.centroid()
    OH_file_structures_with_sc_pb.write(
        ">" + transcript_id + "," + str(OH_dist_with_sc_pb) + "\n" + seq + "\n" + OH_centroid_structure_with_sc_pb + "\n")

    """compute WITH soft constraints (H2O)"""
    fc5 = RNA.fold_compound(seq, md)
    addSoftConstraintsToFc(fc5, soft_constraints_h2o)
    (structure_with_sc_h2o, mfe_with_sc_h2o) = fc5.mfe()
    fc5.exp_params_rescale(mfe_with_sc_h2o)
    (pp, pf) = fc5.pf()
    (centroid_structure_with_sc_h2o, dist_with_sc_h2o) = fc5.centroid()
    file_structures_with_sc_h2o.write(
        ">" + transcript_id + "," + str(dist_with_sc_h2o) + "\n" + seq + "\n" + centroid_structure_with_sc_h2o + "\n")

    """compute WITH soft constraints based on cP library only (H2O)"""
    fc6 = RNA.fold_compound(seq, md)
    addSoftConstraintsToFc(fc6, cP_soft_constraints_h2o)
    (cP_structure_with_sc_h2o, cP_mfe_with_sc_h2o) = fc6.mfe()
    fc6.exp_params_rescale(cP_mfe_with_sc_h2o)
    (pp, pf) = fc6.pf()
    (cP_centroid_structure_with_sc_h2o, cP_dist_with_sc_h2o) = fc6.centroid()
    cP_file_structures_with_sc_h2o.write(
        ">" + transcript_id + "," + str(cP_dist_with_sc_h2o) + "\n" + seq + "\n" + cP_centroid_structure_with_sc_h2o + "\n")

    """compute WITH soft constraints based on OH library only (H2O)"""

    fc7 = RNA.fold_compound(seq, md)
    addSoftConstraintsToFc(fc7, OH_soft_constraints_h2o)
    (OH_structure_with_sc_h2o, OH_mfe_with_sc_h2o) = fc7.mfe()
    fc7.exp_params_rescale(OH_mfe_with_sc_h2o)
    (pp, pf) = fc7.pf()
    (OH_centroid_structure_with_sc_h2o, OH_dist_with_sc_h2o) = fc7.centroid()
    OH_file_structures_with_sc_h2o.write(
        ">" + transcript_id + "," + str(OH_dist_with_sc_h2o) + "\n" + seq + "\n" + OH_centroid_structure_with_sc_h2o + "\n")

    """compute WITH soft constraints (both)"""
    fc8 = RNA.fold_compound(seq, md)
    addSoftConstraintsToFc(fc8, soft_constraints_both)
    (structure_with_sc_both, mfe_with_sc_both) = fc8.mfe()
    fc8.exp_params_rescale(mfe_with_sc_both)
    (pp, pf) = fc8.pf()
    (centroid_structure_with_sc_both, dist_with_sc_both) = fc8.centroid()
    file_structures_with_sc_both.write(
        ">" + transcript_id + "," + str(dist_with_sc_both) + "\n" + seq + "\n" + centroid_structure_with_sc_both + "\n")

    """compute WITH soft constraints based on cP library only (both)"""
    fc9 = RNA.fold_compound(seq, md)
    addSoftConstraintsToFc(fc9, cP_soft_constraints_both)
    (cP_structure_with_sc_both, cP_mfe_with_sc_both) = fc9.mfe()
    fc9.exp_params_rescale(cP_mfe_with_sc_both)
    (pp, pf) = fc9.pf()
    (cP_centroid_structure_with_sc_both, cP_dist_with_sc_both) = fc9.centroid()
    cP_file_structures_with_sc_both.write(
        ">" + transcript_id + "," + str(cP_dist_with_sc_both) + "\n" + seq + "\n" + cP_centroid_structure_with_sc_both + "\n")

    """compute WITH soft constraints based on OH library only (both)"""

    fc10 = RNA.fold_compound(seq, md)
    addSoftConstraintsToFc(fc10, OH_soft_constraints_both)
    (OH_structure_with_sc_both, OH_mfe_with_sc_both) = fc10.mfe()
    fc10.exp_params_rescale(OH_mfe_with_sc_both)
    (pp, pf) = fc10.pf()
    (OH_centroid_structure_with_sc_both, OH_dist_with_sc_both) = fc10.centroid()
    OH_file_structures_with_sc_both.write(
        ">" + transcript_id + "," + str(OH_dist_with_sc_both) + "\n" + seq + "\n" + OH_centroid_structure_with_sc_both + "\n")





for transcript_id in data:
    if transcript_id in ["16S_ribosomal_RNA", "23S_ribosomal_RNA"]:
        continue
    foldMfe(transcript_id, mfe_structures_no_sc,
            mfe_structures_with_sc_pb, mfe_structures_with_sc_h2o,
            mfe_structures_with_sc_both, soft_constraints_output_pb, soft_constraints_output_h2o,
            soft_constraints_output_both, cP_mfe_structures_with_sc_pb, cP_mfe_structures_with_sc_h2o,
            cP_mfe_structures_with_sc_both, cP_soft_constraints_output_pb, cP_soft_constraints_output_h2o,
            cP_soft_constraints_output_both, OH_mfe_structures_with_sc_pb, OH_mfe_structures_with_sc_h2o,
            OH_mfe_structures_with_sc_both, OH_soft_constraints_output_pb, OH_soft_constraints_output_h2o,
            OH_soft_constraints_output_both)
    foldCentroid(transcript_id, centroid_structures_no_sc,
            centroid_structures_with_sc_pb, centroid_structures_with_sc_h2o,
            centroid_structures_with_sc_both, soft_constraints_output_pb, soft_constraints_output_h2o,
            soft_constraints_output_both, cP_centroid_structures_with_sc_pb, cP_centroid_structures_with_sc_h2o,
            cP_centroid_structures_with_sc_both, cP_soft_constraints_output_pb, cP_soft_constraints_output_h2o,
            cP_soft_constraints_output_both, OH_centroid_structures_with_sc_pb, OH_centroid_structures_with_sc_h2o,
            OH_centroid_structures_with_sc_both, OH_soft_constraints_output_pb, OH_soft_constraints_output_h2o,
            OH_soft_constraints_output_both)
    print('finished folding ' + transcript_id)

mfe_structures_no_sc.close()

mfe_structures_with_sc_pb.close()
mfe_structures_with_sc_h2o.close()
mfe_structures_with_sc_both.close()
cP_mfe_structures_with_sc_pb.close()
cP_mfe_structures_with_sc_h2o.close()
cP_mfe_structures_with_sc_both.close()

OH_mfe_structures_with_sc_pb.close()
OH_mfe_structures_with_sc_h2o.close()
OH_mfe_structures_with_sc_both.close()

centroid_structures_no_sc.close()

centroid_structures_with_sc_pb.close()
centroid_structures_with_sc_h2o.close()
centroid_structures_with_sc_both.close()

cP_centroid_structures_with_sc_pb.close()
cP_centroid_structures_with_sc_h2o.close()
cP_centroid_structures_with_sc_both.close()

OH_centroid_structures_with_sc_pb.close()
OH_centroid_structures_with_sc_h2o.close()
OH_centroid_structures_with_sc_both.close()


soft_constraints_output_pb.close()
soft_constraints_output_h2o.close()
soft_constraints_output_both.close()

cP_soft_constraints_output_pb.close()
cP_soft_constraints_output_h2o.close()
cP_soft_constraints_output_both.close()

OH_soft_constraints_output_pb.close()
OH_soft_constraints_output_h2o.close()
OH_soft_constraints_output_both.close()



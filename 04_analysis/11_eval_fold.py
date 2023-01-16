import sys
import numpy as np
from re import finditer
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from tabulate import tabulate
plt.rcParams.update({'font.size': 13})
plt.rcParams.update({'axes.titlesize': 13})
plt.rc('xtick', labelsize=13)
plt.rc('ytick', labelsize=13)

prediction_no_sc_file = sys.argv[1]
prediction_with_sc_file = sys.argv[2]
ref_structures_file = sys.argv[3]
alg = sys.argv[4]                   #   mfe/centroid/..
output_folder = sys.argv[5]
plotname = sys.argv[6]

"""store the structures in dicts"""

prediction_no_sc = {}
prediction_with_sc = {}
ref_structures = {}

line_counter=0
current_id = ""
with open(prediction_no_sc_file) as file:
    for line in file:
        line_counter += 1
        line=line.rstrip()
        if line.startswith(">"):        #transcript_name
            name=line[1:].split(",")[0]
            line_counter = 1
            prediction_no_sc[name]=[]
            current_id=name
        if line_counter == 2:           # RNA sequence
            prediction_no_sc[current_id].append(line)
        if line_counter == 3:           # structure
            prediction_no_sc[current_id].append(line)
    line_counter = 0
    current_id = ""


with open(prediction_with_sc_file) as file:
    for line in file:
        line_counter += 1
        line=line.rstrip()
        if line.startswith(">"):        #transcript_name
            name = line[1:].split(",")[0]
            line_counter = 1
            prediction_with_sc[name] = []
            current_id = name
        if line_counter == 2:           # RNA sequence
            prediction_with_sc[current_id].append(line)
        if line_counter == 3:           # structure
            prediction_with_sc[current_id].append(line)
    line_counter = 0
    current_id = ""

with open(ref_structures_file) as file:
    for line in file:
        line_counter += 1
        line=line.rstrip()
        if line.startswith(">"):        #transcript_name
            line_counter = 1
            ref_structures[line[1:]]=[]
            current_id=line[1:]
        if line_counter == 2:           # RNA sequence
            ref_structures[current_id].append(line)
        if line_counter == 3:           # structure
           ref_structures[current_id].append(line)
    line_counter = 0
    current_id = ""


def calc_MCC_F_val(prediction, actual):
    TP = 0
    TN = 0
    FP = 0
    FN = 0
    prediction_bp=set()
    reference_bp=set()
    
    for i, (a, b) in enumerate(zip(prediction, actual)):
            if a != -1:
                if i<a:
                    prediction_bp.add((i,a))
                else:
                    prediction_bp.add((a,i))
            if b != -1:
                if i<b:
                    reference_bp.add((i,b))
                else:
                    reference_bp.add((b,i))
    for (i,j) in reference_bp:
        if (i,j) in prediction_bp:
            TP += 1
        else:
            FN+=1
    for i in range(0,len(actual)-1):
        for j in range (i+1, len(actual)):
            if (i,j) not in reference_bp and (i,j) not in prediction_bp:
                TN += 1

    for (i,j) in prediction_bp:
        if (i,j) not in reference_bp:
            FP+=1

    PPV = TP / (TP + FP)
    Sensitivity = TP / (TP + FN)
    MCC = ((TP * TN) - (FP * FN)) / np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    F_Value = 0.5 * (PPV + Sensitivity)

    return MCC, F_Value


def parse_dot_bracket(input):
    output = np.full(len(input), -1)
    more = True
    while more:
        more = False

        #finds matched parenthesis
        for x in finditer(r"\([^()]*\)", input):
            more = True
            output[x.start()] = x.end()-1
            output[x.end()-1] = x.start()

            #its recursive...
            input=input[0:x.start()] + "." + input[x.start()+1:x.end()-1] + "." + input[x.end():]

    return output


def calc_overall_MCC_F_val_PPV_Sensitivity(predictions, references):
    PPVs=[]
    Sensitivities=[]
    F_vals=[]
    MCCs=[]   
    for transcript_name in references.keys():
        TP = 0
        TN = 0
        FP = 0
        FN = 0
        prediction_bp=set()
        reference_bp=set()
        reference = parse_dot_bracket(references[transcript_name][1])
        prediction = parse_dot_bracket(predictions[transcript_name][1]) 
    
    
        for i, (a, b) in enumerate(zip(prediction, reference)):
            if a != -1:
                if i<a:
                    prediction_bp.add((i,a))
                else:
                    prediction_bp.add((a,i))
            if b != -1:
                if i<b:
                    reference_bp.add((i,b))
                else:
                    reference_bp.add((b,i))
        
        for (i,j) in reference_bp:
            if (i,j) in prediction_bp:
                TP += 1
            else:
                FN+=1
        for i in range(0,len(reference)-1):
            for j in range (i+1, len(reference)):
                if (i,j) not in reference_bp and (i,j) not in prediction_bp:
                    TN += 1
        for (i,j) in prediction_bp:
            if (i,j) not in reference_bp:
                FP+=1

        PPV = TP / (TP + FP)
        Sensitivity = TP / (TP + FN)
        MCC = ((TP * TN) - (FP * FN)) / np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
        F_Value = 0.5 * (PPV + Sensitivity)

        PPVs.append(PPV)
        Sensitivities.append(Sensitivity)
        MCCs.append(MCC)
        F_vals.append(F_Value)

    MCC_all=sum(MCCs)/len(references.keys())
    F_Value_all=sum(F_vals)/len(references.keys())
    PPV_all=sum(PPVs)/len(references.keys())
    Sensitivity_all=sum(Sensitivities)/len(references.keys())
    
    return MCC_all, F_Value_all, PPV_all, Sensitivity_all

bad_candidates=[]
good_candidates=[]
for transcript_name in prediction_no_sc.keys():
    reference = parse_dot_bracket(ref_structures[transcript_name][1])
    prediction_no_sc_parsed = parse_dot_bracket(prediction_no_sc[transcript_name][1])
    prediction_with_sc_parsed = parse_dot_bracket(prediction_with_sc[transcript_name][1])
    MCC_prediction_no_sc, F_prediction_no_sc = calc_MCC_F_val(prediction_no_sc_parsed, reference)
    MCC_prediction_with_sc, F_prediction_with_sc = calc_MCC_F_val(prediction_with_sc_parsed, reference)

    if MCC_prediction_no_sc > MCC_prediction_with_sc:
        bad_candidates.append(transcript_name)
    if MCC_prediction_no_sc < MCC_prediction_with_sc:
        good_candidates.append(transcript_name)

    prediction_no_sc[transcript_name].append(MCC_prediction_no_sc)
    prediction_no_sc[transcript_name].append(F_prediction_no_sc)
    prediction_with_sc[transcript_name].append(MCC_prediction_with_sc)
    prediction_with_sc[transcript_name].append(F_prediction_with_sc)

#comparison scatterplot MCC
fig = plt.figure(figsize=(6,6))
ax = fig.gca()
ax.plot([0, 1], [0, 1], transform=ax.transAxes, color="#0571b0")
for transcript_name in prediction_no_sc.keys():
    ax.scatter(prediction_no_sc[transcript_name][2], prediction_with_sc[transcript_name][2], color="#ca0020", alpha=0.7)
    if transcript_name == "tRNA_Ile_GAU":
        ax.annotate(r'$\mathregular{tRNA^{Ile(GAU)}}$', (prediction_no_sc[transcript_name][2]-0.04, prediction_with_sc[transcript_name][2]+0.02), fontsize=12 ,  color = "#0571b0") 
    # if transcript_name == "tRNA_Leu_CAG_1":
    #     ax.annotate(r'$\mathregular{tRNA^{Leu(CAG)}}$', (prediction_no_sc[transcript_name][2]-0.04, prediction_with_sc[transcript_name][2]+0.02), fontsize=12 ,  color = "#0571b0") 
    elif transcript_name == "5S_ribosomal_RNA":
        ax.annotate("5S rRNA", (prediction_no_sc[transcript_name][2]-0.06, prediction_with_sc[transcript_name][2]+0.02), fontsize=12, color = "#0571b0")
    elif transcript_name == "23S_ribosomal_RNA_domain_IV":
        ax.annotate("  23S rRNA \n(domain IV)", (prediction_no_sc[transcript_name][2]-0.09, prediction_with_sc[transcript_name][2]-0.08), fontsize=12,  color = "#0571b0")

width=0.03
height=0.03
# for (x, y) in [(prediction_no_sc["tRNA_Leu_CAG_1"][2],prediction_with_sc["tRNA_Leu_CAG_1"][2]),(prediction_no_sc["tRNA_Ile_GAU"][2],prediction_with_sc["tRNA_Ile_GAU"][2]),(prediction_no_sc["5S_ribosomal_RNA"][2],prediction_with_sc["5S_ribosomal_RNA"][2]),(prediction_no_sc["23S_ribosomal_RNA_domain_IV"][2],prediction_with_sc["23S_ribosomal_RNA_domain_IV"][2])]:
for (x, y) in [(prediction_no_sc["tRNA_Ile_GAU"][2],prediction_with_sc["tRNA_Ile_GAU"][2]),(prediction_no_sc["5S_ribosomal_RNA"][2],prediction_with_sc["5S_ribosomal_RNA"][2]),(prediction_no_sc["23S_ribosomal_RNA_domain_IV"][2],prediction_with_sc["23S_ribosomal_RNA_domain_IV"][2])]:
    ax.add_patch(Rectangle(
        xy=(x-width/2, y-height/2) ,width=width, height=height,
        linewidth=1.5, color= "#0571b0", fill=False))

ax.set_xlabel("MCC (without constraints)")
ax.set_ylabel("MCC (with probing signal)")
plt.tight_layout()
plt.savefig(output_folder+"/"+plotname + "_MCC_scatter.pdf", dpi=600)
plt.close()

# print out overall MCC, F-val, PPV, Sensitivity
MCC_prediction_no_sc_all, F_val_prediction_no_sc_all, PPV_prediction_no_sc_all, Sens_prediction_no_sc_all = calc_overall_MCC_F_val_PPV_Sensitivity(prediction_no_sc, ref_structures)
MCC_prediction_with_sc_all, F_val_prediction_with_sc_all, PPV_prediction_with_sc_all, Sens_prediction_with_sc_all = calc_overall_MCC_F_val_PPV_Sensitivity(prediction_with_sc, ref_structures)

stats=[[alg, MCC_prediction_no_sc_all, F_val_prediction_no_sc_all, PPV_prediction_no_sc_all, Sens_prediction_no_sc_all], 
[alg+" with sc", MCC_prediction_with_sc_all, F_val_prediction_with_sc_all, PPV_prediction_with_sc_all, Sens_prediction_with_sc_all]]

print("\nstatistics folding accuracy")
print(tabulate(stats, headers=["prediction", "MCC", "F-val", "PPV", "Sensitivity"]))
print("bad_candidates: "+str(bad_candidates))
print("\n")
print("good_candidates: "+str(good_candidates))
print("\n")
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

plt.rcParams.update({'font.size': 19})
plt.rcParams.update({'axes.titlesize': 19})
plt.rc('legend', fontsize=19)
plt.rc('xtick', labelsize=19)
plt.rc('ytick', labelsize=19)

signal_file = sys.argv[1]
benchmark_gff_file = sys.argv[2]
output_folder = sys.argv[3]


"""read in benchmark set"""
benchmark_transcripts = []

with open(benchmark_gff_file, 'r') as file:
    for line in file:
        line = line.rstrip()
        columns = line.split("\t")
        transcript_id = columns[8]
        benchmark_transcripts.append(transcript_id)




"""Pb treated samples - analyze mean values of replicates"""

signal = {}
paired_signal_cP = []
unpaired_signal_cP = []
paired_signal_OH = []
unpaired_signal_OH = []


with open(signal_file) as file:
    for line in file:
        line=line.rstrip()
        columns = line.split("\t")
        transcript_id, nc, position, structure, sequence, library, library_type = columns[5],columns[7],int(columns[8]),columns[9],columns[3],columns[10], columns[11]
        if transcript_id not in signal:
            signal[transcript_id] = {}
        if library not in signal[transcript_id]:
            signal[transcript_id][library] = {}
        if position not in signal[transcript_id][library]:
            signal[transcript_id][library][position] = [nc, structure]

for transcript_id in benchmark_transcripts:
    for pos in signal[transcript_id]["1-cP"]:
        if signal[transcript_id]["1-cP"][pos][0] != "NA" and signal[transcript_id]["2-cP"][pos][0] != "NA":
            if signal[transcript_id]["1-cP"][pos][1]=="0":
                unpaired_signal_cP.append(np.mean([float(signal[transcript_id]["1-cP"][pos][0]),float(signal[transcript_id]["2-cP"][pos][0])]))
            elif signal[transcript_id]["1-cP"][pos][1]=="1":
                paired_signal_cP.append(np.mean([float(signal[transcript_id]["1-cP"][pos][0]),float(signal[transcript_id]["2-cP"][pos][0])]))

        if signal[transcript_id]["1-OH"][pos][0] != "NA" and signal[transcript_id]["2-OH"][pos][0] != "NA":
            if signal[transcript_id]["1-OH"][pos][1]=="0":
                unpaired_signal_OH.append(np.mean([float(signal[transcript_id]["1-OH"][pos][0]),float(signal[transcript_id]["2-OH"][pos][0])]))
            elif signal[transcript_id]["1-OH"][pos][1]=="1":
                paired_signal_OH.append(np.mean([float(signal[transcript_id]["1-OH"][pos][0]),float(signal[transcript_id]["2-OH"][pos][0])]))


fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(10,4))
"cP"

ax1.hist(paired_signal_cP, bins=50, color="#0571b0",linewidth=1.5, histtype="step", label="paired")
ax1.hist(unpaired_signal_cP, bins=50, color="#ca0020",linewidth=1.5, histtype="step", label="unpaired")
ax1.set_xlabel(r'$S$')
ax1.set_ylabel('$\mathrm{frequency}$')
ax1.set_title("2'3'-cP")
ax1.set_yscale('log')

"OH"
ax2.hist(paired_signal_OH, bins=50, color="#0571b0", linewidth=1.5, histtype="step", label="paired")
ax2.hist(unpaired_signal_OH, bins=50, color="#ca0020", linewidth=1.5, histtype="step", label="unpaired")
ax2.set_title("5'-OH")
ax2.set_xlabel(r'$S$')
ax2.set_yscale('log')
plt.legend()
plt.tight_layout()
plt.savefig(output_folder + "/"+"histograms_Pb.pdf", bbox_inches='tight')
plt.close()

print(stats.mannwhitneyu(paired_signal_cP,unpaired_signal_cP))
print(stats.mannwhitneyu(paired_signal_OH,unpaired_signal_OH))



"""H2O treated samples"""

signal = {}
paired_signal = []
unpaired_signal = []

with open(signal_file) as file:
    for line in file:
        line=line.rstrip()
        columns = line.split("\t")
        transcript_id, nc, structure, library =columns[5],  columns[7], columns[9], columns[10]
        if library not in signal:
            signal[library]=[[],[]]
        if structure != "-":
            if structure == "0" and nc != "NA":
                signal[library][0].append(float(nc))
                unpaired_signal.append(float(nc))
            elif structure == "1" and nc != "NA":
                signal[library][1].append(float(nc))
                paired_signal.append(float(nc))



fig2, (ax3, ax4) = plt.subplots(1, 2, sharey=True, figsize=(10,4))
"3-cP"

ax3.hist(signal["3-cP"][1], bins=50, color="#0571b0",linewidth=1.5, histtype="step", label="paired")
ax3.hist(signal["3-cP"][0], bins=50, color="#ca0020",linewidth=1.5, histtype="step", label="unpaired")
ax3.set_xlabel(r'$S$')
ax3.set_ylabel('$\mathrm{frequency}$')
ax3.set_title("2'3'-cP $\mathregular{Pb^{2\!+}}\!(-)$")
ax3.set_yscale('log')

"3-OH"
ax4.hist(signal["3-OH"][1], bins=50, color="#0571b0", linewidth=1.5, histtype="step", label="paired")
ax4.hist(signal["3-OH"][0], bins=50, color="#ca0020", linewidth=1.5, histtype="step", label="unpaired")
ax4.set_title("5'-OH $\mathregular{Pb^{2\!+}}\!(-)$")
ax4.set_xlabel(r'$S$')
ax4.set_yscale('log')


plt.legend()
plt.tight_layout()
plt.savefig(output_folder + "/"+"histograms_H2O.pdf", bbox_inches='tight')
plt.close()
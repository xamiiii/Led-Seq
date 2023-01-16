import sys
from scipy import stats
import matplotlib.pyplot as plt

signal_file = sys.argv[1]
signal = {}
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
            signal[transcript_id][library][position] = nc


ncs_1_cP=[]
ncs_2_cP=[]
ncs_1_OH=[]
ncs_2_OH=[]
pos_counter_cP=0
pos_counter_OH=0
for transcript_id in signal:
    if "1-cP" in signal[transcript_id] and "2-cP" in signal[transcript_id]:
        for pos in signal[transcript_id]["1-cP"]:
            if signal[transcript_id]["1-cP"][pos] != "NA" and signal[transcript_id]["2-cP"][pos] != "NA":
                pos_counter_cP+=1
                ncs_1_cP.append(float(signal[transcript_id]["1-cP"][pos]))
                ncs_3_cP.append(float(signal[transcript_id]["2-cP"][pos]))
    if "1-OH" in signal[transcript_id] and "2-OH" in signal[transcript_id]:
        for pos in signal[transcript_id]["1-OH"]:
            if signal[transcript_id]["1-OH"][pos] != "NA" and signal[transcript_id]["2-OH"][pos] != "NA":
                pos_counter_OH+=1
                ncs_1_OH.append(float(signal[transcript_id]["1-OH"][pos]))
                ncs_2_OH.append(float(signal[transcript_id]["2-OH"][pos]))

print(stats.pearsonr(ncs_1_cP, ncs_2_cP), pos_counter_cP)
plt.scatter(ncs_1_cP, ncs_2_cP)
plt.savefig("./correlation_cP_replicates.pdf")

print(stats.pearsonr(ncs_1_OH, ncs_2_OH), pos_counter_OH)
plt.scatter(ncs_1_OH, ncs_2_OH)
plt.savefig("./correlation_OH_replicates.pdf")

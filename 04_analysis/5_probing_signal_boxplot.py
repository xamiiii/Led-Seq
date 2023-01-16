import sys
import matplotlib.pyplot as plt

signal_file = sys.argv[1]
output_folder = sys.argv[2]

signal = {}
paired_signal=[]
unpaired_signal=[]

with open(signal_file) as file:
    for line in file:
        line=line.rstrip()
        columns = line.split("\t")
        nc, structure, library = columns[7], columns[9], columns[10]
        if library not in signal:
            signal[library]=[[],[]]
        if structure != "-":
            if structure == "0" and nc != "NA":
                signal[library][0].append(float(nc))
                unpaired_signal.append(float(nc))
            elif structure == "1" and nc != "NA":
                signal[library][1].append(float(nc))
                paired_signal.append(float(nc))
labels = ['all_libraries']
A=[unpaired_signal, paired_signal]
position1=1
position2=2

fig = plt.figure()
ax = fig.add_subplot()
bp = ax.boxplot(A, positions = [position1, position2], widths=0.6)
for flier in bp['fliers']:
    flier.set(marker='.')

for library in signal:
    position1+=3
    position2+=3
    bp = ax.boxplot(signal[library], positions = [position1, position2], widths=0.6)
    labels.append(library)
    for flier in bp['fliers']:
        flier.set(marker='.')

x_ticks=[]
for i, v in enumerate(labels):
    position=1.5+i*3
    x_ticks.append(position)
ax.set_xticks(x_ticks)
ax.set_xticklabels(labels, rotation='vertical')
ax.set_ylabel('normalized reactivity')
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(output_folder+"/boxplot.pdf")
plt.close()




import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import logomaker as lm
import pandas as pd
plt.rcParams.update({'font.size': 19})
plt.rcParams.update({'axes.titlesize': 19})
plt.rc('xtick', labelsize=19)
plt.rc('ytick', labelsize=19)


signal_file = sys.argv[1]
covered_transcripts=set()
signal = {}
with open(signal_file) as file:
    for line in file:
        line=line.rstrip()
        columns = line.split("\t")
        transcript_id, nc, position, structure, sequence, library, library_type, base, neighboring_base = columns[5], columns[7] ,int(columns[8]),columns[9],columns[3],columns[10], columns[11], columns[3], columns[4]
        covered_transcripts.add(transcript_id)
        if library not in signal:
            signal[library] = {}
        if transcript_id not in signal[library]:
            signal[library][transcript_id] = {}
        signal[library][transcript_id][position] = [nc, base]


"""for the first 100 Codons of CDS"""
mean_ncs_cP=[]
mean_ncs_OH=[]
codons_cP=[]
codons_OH=[]


for pos in range(1, 351):
    nc_sum_cP = 0
    counter_cP = 0
    nc_sum_OH = 0
    counter_OH = 0
    current_codons_cP=[]
    current_codons_OH = []

    for transcript_id in covered_transcripts:
        if transcript_id in signal["1-cP"] and transcript_id in signal["2-cP"]:
            # if 350 not in signal["1-cP"][transcript_id] or 350 not in signal["2-cP"][transcript_id]:
            #     continue
            if pos in signal["1-cP"][transcript_id] and pos in signal["2-cP"][transcript_id]:
                current_codons_cP.append(signal["1-cP"][transcript_id][pos][1])
                if signal["1-cP"][transcript_id][pos][0] != "NA" and signal["1-cP"][transcript_id][pos][0] != "n.p." and signal["2-cP"][transcript_id][pos][0] != "NA" and signal["2-cP"][transcript_id][pos][0] != "n.p.":
                    nc_sum_cP += (float(signal["1-cP"][transcript_id][pos][0])+float(signal["2-cP"][transcript_id][pos][0]))/2
                    counter_cP += 1
            

        if transcript_id in signal["1-OH"] and transcript_id in signal["2-OH"]:
            # if 350 not in signal["1-OH"][transcript_id] or 350 not in signal["2-OH"][transcript_id]:
            #     continue
            if pos in signal["1-OH"][transcript_id] and pos in signal["2-OH"][transcript_id]:
                current_codons_OH.append(signal["1-OH"][transcript_id][pos][1])
                if signal["1-OH"][transcript_id][pos][0] != "NA" and signal["1-OH"][transcript_id][pos][0] != "n.p." and signal["2-OH"][transcript_id][pos][0] != "NA" and signal["2-OH"][transcript_id][pos][0] != "n.p.":
                    nc_sum_OH += (float(signal["1-OH"][transcript_id][pos][0])+float(signal["2-OH"][transcript_id][pos][0]))/2
                    counter_OH += 1

    if counter_cP > 0:
        mean_nc_cP = nc_sum_cP / counter_cP
        mean_ncs_cP.append(mean_nc_cP)
    else:
        mean_ncs_cP.append(0)
    codons_cP.append(current_codons_cP)

    if counter_OH > 0:
        mean_nc_OH = nc_sum_OH / counter_OH
        mean_ncs_OH.append(mean_nc_OH)
    else:
        mean_ncs_OH.append(0)
    codons_OH.append(current_codons_OH)
    print(pos, counter_cP, counter_OH)

codonbases_cP={"A":[0,0,0], "U":[0,0,0], "C":[0,0,0], "G":[0,0,0]}
codonbases_OH={"A":[0,0,0], "U":[0,0,0], "C":[0,0,0], "G":[0,0,0]}

for i, position in enumerate(codons_cP[50:230][0::3]):
    # if i==0:
    #     print(position)
    for base in position:
        codonbases_cP[base][0]+=1
for i, position in enumerate(codons_cP[50:230][1::3]):
    # if i==0:
    #     print(position)
    for base in position:
        codonbases_cP[base][1]+=1
for i, position in enumerate(codons_cP[50:230][2::3]):
    # if i==0:
    #     print(position)
    for base in position:
        codonbases_cP[base][2]+=1

"""convert counts to frequencies"""
for pos in [0,1,2]:
    sum=codonbases_cP["A"][pos]+codonbases_cP["U"][pos]+codonbases_cP["C"][pos]+codonbases_cP["G"][pos]
    codonbases_cP["A"][pos]=codonbases_cP["A"][pos]/sum
    codonbases_cP["U"][pos] = codonbases_cP["U"][pos] / sum
    codonbases_cP["C"][pos] = codonbases_cP["C"][pos] / sum
    codonbases_cP["G"][pos] = codonbases_cP["G"][pos] / sum

"""draw logo"""
# create color scheme
color_scheme = {
    'A' : "#f4a582",
    'C' : "#67a9cf",
    'G' : "#ca0020",
    'U': "#2166ac"
}

logo_cP = lm.Logo(pd.DataFrame(data=codonbases_cP), figsize=(3,3), color_scheme=color_scheme)#'colorblind_safe')
logo_cP.ax.set_xlabel('codon position')
logo_cP.ax.set_ylabel('frequency')
plt.xticks(ticks=[0,1,2], labels=["1", "2", "3"])
plt.tight_layout()
plt.savefig("./sequence_logo_cP_100codons.pdf")


for i, position in enumerate(codons_OH[50:230][0::3]):
    # if i==0:
        # print(position)
    for base in position:
        codonbases_OH[base][0]+=1
for i, position in enumerate(codons_OH[50:230][1::3]):
    # if i==0:
        # print(position)
    for base in position:
        codonbases_OH[base][1]+=1
for i, position in enumerate(codons_OH[50:230][2::3]):
    # if i==0:
        # print(position)
    for base in position:
        codonbases_OH[base][2]+=1

"""convert counts to frequencies"""
for pos in [0,1,2]:
    sum=codonbases_OH["A"][pos]+codonbases_OH["U"][pos]+codonbases_OH["C"][pos]+codonbases_OH["G"][pos]
    codonbases_OH["A"][pos]=codonbases_OH["A"][pos]/sum
    codonbases_OH["U"][pos] = codonbases_OH["U"][pos] / sum
    codonbases_OH["C"][pos] = codonbases_OH["C"][pos] / sum
    codonbases_OH["G"][pos] = codonbases_OH["G"][pos] / sum

"""draw logo"""
logo_OH = lm.Logo(pd.DataFrame(data=codonbases_OH), figsize=(3,3),color_scheme=color_scheme)
logo_OH.ax.set_xlabel('codon position')
logo_OH.ax.set_ylabel('frequency')
plt.xticks(ticks=[0,1,2], labels=["1", "2", "3"])
plt.tight_layout()
plt.savefig("./sequence_logo_OH_100codons.pdf")


fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4,  sharey=True, figsize=(10,4))
medianprops = dict(color="#ca0020")
ax1_bplot=ax1.boxplot([mean_ncs_cP[14:50][0::3],mean_ncs_cP[14:50][1::3],mean_ncs_cP[14:50][2::3]], showfliers=False, patch_artist=True, medianprops=medianprops, widths=(0.8, 0.8, 0.8))
ax2_bplot=ax2.boxplot([mean_ncs_cP[50:230][0::3], mean_ncs_cP[50:230][1::3], mean_ncs_cP[50:230][2::3]], showfliers=False, patch_artist=True, medianprops=medianprops, widths=(0.8, 0.8, 0.8))

ax3_bplot=ax3.boxplot([mean_ncs_OH[2:50][0::3], mean_ncs_OH[2:50][1::3], mean_ncs_OH[2:50][2::3]], showfliers=False, patch_artist=True, medianprops=medianprops, widths=(0.8, 0.8, 0.8))
ax4_bplot=ax4.boxplot([mean_ncs_OH[50:230][0::3], mean_ncs_OH[50:230][1::3], mean_ncs_OH[50:230][2::3]], showfliers=False, patch_artist=True, medianprops=medianprops, widths=(0.8, 0.8, 0.8))

ax1.set_title("2'3'-cP", loc="right")
ax3.set_title("5'-OH", loc="right")
for boxplot in (ax2_bplot, ax4_bplot):
    for patch in boxplot['boxes']:
        patch.set_facecolor("#f4a582")
for boxplot in (ax1_bplot, ax3_bplot):
    for patch in boxplot['boxes']:
        patch.set_facecolor(color="#67a9cf")

ax1.set_ylabel(r'$S$')
ax1.set_xlabel("5'-UTR")
ax2.set_xlabel("CDS")
ax3.set_xlabel("5'-UTR")
ax4.set_xlabel("CDS")
ax1.set_ylim(top=0.9)
plt.tight_layout()
plt.savefig("./peridicity_both_60codons.pdf")
plt.close()



print("UTR cP 12 codons 1-2: ", stats.ttest_ind(mean_ncs_cP[14:50][0::3],mean_ncs_cP[14:50][1::3]))
print("UTR cP 12 codons 1-3: ", stats.ttest_ind(mean_ncs_cP[14:50][0::3],mean_ncs_cP[14:50][2::3]))
print("UTR cP 12 codons 2-3: ", stats.ttest_ind(mean_ncs_cP[14:50][1::3],mean_ncs_cP[14:50][2::3]))
print("CDS cP 60 codons 1-2: ", stats.ttest_ind(mean_ncs_cP[50:230][0::3],mean_ncs_cP[50:230][1::3]))
print("CDS cP 60 codons 1-3: ", stats.ttest_ind(mean_ncs_cP[50:230][0::3],mean_ncs_cP[50:230][2::3]))
print("CDS cP 60 codons 2-3: ", stats.ttest_ind(mean_ncs_cP[50:230][1::3],mean_ncs_cP[50:230][2::3]))
print("CDS cP 100 codons 1-2: ", stats.ttest_ind(mean_ncs_cP[50:350][0::3],mean_ncs_cP[50:350][1::3]))
print("CDS cP 100 codons 1-3: ", stats.ttest_ind(mean_ncs_cP[50:350][0::3],mean_ncs_cP[50:350][2::3]))
print("CDS cP 100 codons 2-3: ", stats.ttest_ind(mean_ncs_cP[50:350][1::3],mean_ncs_cP[50:350][2::3]))
print("UTR OH 16 codons 1-2: ", stats.ttest_ind(mean_ncs_OH[2:50][0::3],mean_ncs_OH[2:50][1::3]))
print("UTR OH 16 codons 1-3: ", stats.ttest_ind(mean_ncs_OH[2:50][0::3],mean_ncs_OH[2:50][2::3]))
print("UTR OH 16 codons 2-3: ", stats.ttest_ind(mean_ncs_OH[2:50][1::3],mean_ncs_OH[2:50][2::3]))
print("CDS OH 60 codons 1-2: ", stats.ttest_ind(mean_ncs_OH[50:230][0::3],mean_ncs_OH[50:230][1::3]))
print("CDS OH 60 codons 1-3: ", stats.ttest_ind(mean_ncs_OH[50:230][0::3],mean_ncs_OH[50:230][2::3]))
print("CDS OH 60 codons 2-3: ", stats.ttest_ind(mean_ncs_OH[50:230][1::3],mean_ncs_OH[50:230][2::3]))
print("CDS OH 100 codons 1-2: ", stats.ttest_ind(mean_ncs_OH[50:350][0::3],mean_ncs_OH[50:350][1::3]))
print("CDS OH 100 codons 1-3: ", stats.ttest_ind(mean_ncs_OH[50:350][0::3],mean_ncs_OH[50:350][2::3]))
print("CDS OH 100 codons 2-3: ", stats.ttest_ind(mean_ncs_OH[50:350][1::3],mean_ncs_OH[50:350][2::3]))


fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True,  figsize=(10,6))#sharey=True,
ax1.plot([47,48],mean_ncs_cP[49:51],'.-', color="lightgrey", linewidth=1.15, markersize=4)
ax1.plot(np.arange(0,11),mean_ncs_cP[2:13],'.-', color="lightgrey", linewidth=1.15, markersize=4)
ax1.plot(np.arange(48,228), mean_ncs_cP[50:230], '.-', color="#f4a582", label="CDS", linewidth=1.15, markersize=4)
ax1.plot([48,107],[np.mean(mean_ncs_cP[50:110]),np.mean(mean_ncs_cP[50:110])], color="#ca0020", ls='--', linewidth=1.25)
ax1.plot([108,228],[np.mean(mean_ncs_cP[110:230]),np.mean(mean_ncs_cP[110:230])], color="#ca0020", ls='--', linewidth=1.25)
ax1.plot(np.arange(10,48), mean_ncs_cP[12:50],'.-', color="#0571b0", label="5'-UTR", linewidth=1.15, markersize=4)
ax1.plot([10,47], [np.mean(mean_ncs_cP[12:50]),np.mean(mean_ncs_cP[12:50])], color="#0571b0", ls='--', linewidth=1.25)
ax1.set_xlim(left=-4, right=233)
ax1.set_ylabel(r'$S$')
ax2.set_xlabel("nucleotide position")
plt.xticks(ticks=[8, 28, 48, 67, 87, 107, 127, 147, 167, 187, 207, 227], labels=[-40,-20, 1, 20, 40, 60, 80, 100, 120, 140, 160, 180])#, 200, 220, 240, 260, 280, 300])
ax1.text(4, 0.96, "2'3'-cP", fontsize=19)
ax2.text(4, 0.665, "5'-OH", fontsize=19)
ax2.set_ylabel(r'$S$')
ax2.plot([47,48],mean_ncs_OH[49:51],'.-', color="lightgrey", linewidth=1.15, markersize=4)
ax2.plot(np.arange(48,228), mean_ncs_OH[50:230], '.-', color="#f4a582", label="CDS", linewidth=1.15, markersize=4)
ax2.plot([48,107],[np.mean(mean_ncs_OH[50:110]),np.mean(mean_ncs_OH[50:110])], color="#ca0020", ls='--', linewidth=1.25)
ax2.plot([108,228],[np.mean(mean_ncs_OH[110:230]),np.mean(mean_ncs_OH[110:230])], color="#ca0020", ls='--', linewidth=1.25)
ax2.plot(np.arange(48), mean_ncs_OH[2:50],'.-', color="#0571b0", label="5'-UTR", linewidth=1.15, markersize=4)
ax2.plot([0,47], [np.mean(mean_ncs_OH[2:50]),np.mean(mean_ncs_OH[2:50])], color="#0571b0", ls='--', linewidth=1.25)
plt.legend(loc="upper right")
plt.tight_layout()
plt.savefig("./nc_aug_mRNA_both_60codons.pdf")
plt.close()

print("cP mean UTR: ",np.mean(mean_ncs_cP[12:50]))
print("cP mean CDS 16 Codons: ",np.mean(mean_ncs_cP[50:98]))
print("cP mean CDS 1.-20. Codon: ",np.mean(mean_ncs_cP[50:110]))
print("cP mean CDS 21.-100 Codons: ",np.mean(mean_ncs_cP[110:350]))
print("cP mean CDS 21.-60. Codons: ",np.mean(mean_ncs_cP[110:230]))
print("cP UTR vs CDS 1_20: ", stats.ttest_ind(mean_ncs_cP[12:50],mean_ncs_cP[50:110]))
print("cP CDS 1_20 vs CDS 21_100: ", stats.ttest_ind(mean_ncs_cP[50:110],mean_ncs_cP[110:350]))
print("cP CDS 1_20 vs CDS 21_60: ", stats.ttest_ind(mean_ncs_cP[50:110],mean_ncs_cP[110:230]))

print("OH mean UTR: ",np.mean(mean_ncs_OH[2:50]))
print("OH mean CDS 16 Codons: ",np.mean(mean_ncs_OH[50:98]))
print("OH mean CDS 1.-20. Codon: ",np.mean(mean_ncs_OH[50:110]))
print("OH mean CDS 21.-100 Codons: ",np.mean(mean_ncs_OH[110:350]))
print("OH mean CDS 21.-60 Codons: ",np.mean(mean_ncs_OH[110:230]))
print("OH UTR vs CDS 1_20: ", stats.ttest_ind(mean_ncs_OH[2:50],mean_ncs_OH[50:110]))
print("OH CDS 1_20 vs CDS 21_100: ", stats.ttest_ind(mean_ncs_OH[50:110],mean_ncs_OH[110:350]))
print("OH CDS 1_20 vs CDS 21_60: ", stats.ttest_ind(mean_ncs_OH[50:110],mean_ncs_OH[110:230]))
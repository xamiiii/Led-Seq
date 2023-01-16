import sys
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'axes.titlesize': 18})
plt.rcParams['xtick.major.pad']= 5
plt.rcParams['axes.titlepad'] = 10
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)
plt.rc('legend', fontsize=18)

library_file_1cP = sys.argv[1]
library_file_2cP = sys.argv[2]
library_file_3cP = sys.argv[3]
library_file_1OH = sys.argv[4]
library_file_2OH = sys.argv[5]
library_file_3OH = sys.argv[6]

bases = {"AA":[0,0,0,0,0,0], "AU": [0,0,0,0,0,0], "AG":[0,0,0,0,0,0], "AC":[0,0,0,0,0,0],
         "UA":[0,0,0,0,0,0], "UU": [0,0,0,0,0,0], "UG":[0,0,0,0,0,0], "UC":[0,0,0,0,0,0],
         "GA":[0,0,0,0,0,0], "GU": [0,0,0,0,0,0], "GG":[0,0,0,0,0,0], "GC":[0,0,0,0,0,0],
         "CA":[0,0,0,0,0,0], "CU": [0,0,0,0,0,0], "CG":[0,0,0,0,0,0], "CC":[0,0,0,0,0,0]}

bases_per_transcript = {"AA":[[],[],[],[],[],[]], "AU": [[],[],[],[],[],[]], "AG":[[],[],[],[],[],[]], "AC":[[],[],[],[],[],[]],
         "UA":[[],[],[],[],[],[]], "UU": [[],[],[],[],[],[]], "UG":[[],[],[],[],[],[]], "UC":[[],[],[],[],[],[]],
         "GA":[[],[],[],[],[],[]], "GU": [[],[],[],[],[],[]], "GG":[[],[],[],[],[],[]], "GC":[[],[],[],[],[],[]],
         "CA":[[],[],[],[],[],[]], "CU": [[],[],[],[],[],[]], "CG":[[],[],[],[],[],[]], "CC":[[],[],[],[],[],[]]}

#simply count the occurrence of dinucleotides in transcriptome
bases_transcriptome = {"AA":0, "AU": 0, "AG":0, "AC":0,
         "UA":0, "UU": 0, "UG":0, "UC":0,
         "GA":0, "GU": 0, "GG":0, "GC":0,
         "CA":0, "CU": 0, "CG":0, "CC":0}
distributions={"AA":[], "AU": [], "AG":[], "AC":[],
         "UA":[], "UU": [], "UG":[], "UC":[],
         "GA":[], "GU": [], "GG":[], "GC":[],
         "CA":[], "CU": [], "CG":[], "CC":[]}

for i, file in enumerate([library_file_1cP,
             library_file_2cP,
             library_file_3cP,
             library_file_1OH,
             library_file_2OH,
             library_file_3OH]):
    with open(file, 'r') as read_file:
        current_transcript_id = ""
        bases_temp = {"AA":0, "AU": 0, "AG": 0, "AC": 0,
                 "UA": 0, "UU": 0, "UG": 0, "UC": 0,
                 "GA": 0, "GU": 0, "GG": 0, "GC": 0,
                 "CA": 0, "CU": 0, "CG": 0, "CC": 0}
        total_count=0
        for line in read_file:
            columns = line.rstrip().split()
            base, neighboring_base, transcript_id, raw_count, normalized_count = columns[3], columns[4], columns[5], int(columns[6]), columns[7]
            if base != "N" and neighboring_base != "N":
                bases[base+neighboring_base][i]+=raw_count
                if normalized_count!="n.p.":
                    distributions[base + neighboring_base].append(float(normalized_count))

                if i == 0:
                    bases_transcriptome[base + neighboring_base] += 1
                if transcript_id == current_transcript_id:
                    bases_temp[base + neighboring_base] += raw_count
                    total_count += raw_count
                else:
                    if normalized_count=="n.p.":
                        continue
                    if current_transcript_id!="":
                        for dinucleotide in bases_temp:
                            if total_count==0:
                                print(base, neighboring_base, transcript_id, raw_count)
                            val=bases_temp[dinucleotide]/total_count*100
                            bases_per_transcript[dinucleotide][i].append(val)
                            if val > 80.0:
                                print(i, current_transcript_id, total_count)
                        bases_temp = {"AA": 0, "AU": 0, "AG": 0, "AC": 0,
                                      "UA": 0, "UU": 0, "UG": 0, "UC": 0,
                                      "GA": 0, "GU": 0, "GG": 0, "GC": 0,
                                      "CA": 0, "CU": 0, "CG": 0, "CC": 0}
                        total_count = 0
                    current_transcript_id=transcript_id
                    bases_temp[base+neighboring_base]+=raw_count
                    total_count+=raw_count

        for dinucleotide in bases_temp:
            val = bases_temp[dinucleotide] / total_count * 100
            bases_per_transcript[dinucleotide][i].append(val)

# for key in distributions:
#     plt.hist(distributions[key], histtype="step", label=key)
# plt.yscale("log")
# plt.legend(ncol=3,  prop={'size': 6}, loc='upper center', bbox_to_anchor=(0.5, -0.1))
# plt.tight_layout()
# plt.savefig("./Dinucleotide_distributions.pdf")

total_count=bases_transcriptome["AA"]+bases_transcriptome["AU"]+bases_transcriptome["AG"]+bases_transcriptome["AC"]+\
            bases_transcriptome["UA"]+bases_transcriptome["UU"]+bases_transcriptome["UG"]+bases_transcriptome["UC"]+\
            bases_transcriptome["GA"]+bases_transcriptome["GU"]+bases_transcriptome["GG"]+bases_transcriptome["GC"]+\
            bases_transcriptome["CA"]+bases_transcriptome["CU"]+bases_transcriptome["CG"]+bases_transcriptome["CC"]
# print(total_count, bases_transcriptome)

dinucleotides=[bases_transcriptome["AA"]/total_count*100, bases_transcriptome["AU"]/total_count*100,bases_transcriptome["AG"]/total_count*100,bases_transcriptome["AC"]/total_count*100,\
            bases_transcriptome["UA"]/total_count*100,bases_transcriptome["UU"]/total_count*100,bases_transcriptome["UG"]/total_count*100,bases_transcriptome["UC"]/total_count*100,\
            bases_transcriptome["GA"]/total_count*100,bases_transcriptome["GU"]/total_count*100,bases_transcriptome["GG"]/total_count*100,bases_transcriptome["GC"]/total_count*100,\
            bases_transcriptome["CA"]/total_count*100,bases_transcriptome["CU"]/total_count*100,bases_transcriptome["CG"]/total_count*100,bases_transcriptome["CC"]/total_count*100]
nucleotides=[(bases_transcriptome["AA"]/total_count*100 + bases_transcriptome["AU"]/total_count*100 + bases_transcriptome["AG"]/total_count*100 + bases_transcriptome["AC"]/total_count*100),\
             (bases_transcriptome["UA"]/total_count*100+bases_transcriptome["UU"]/total_count*100+bases_transcriptome["UG"]/total_count*100+bases_transcriptome["UC"]/total_count*100),\
             (bases_transcriptome["GA"]/total_count*100+bases_transcriptome["GU"]/total_count*100+bases_transcriptome["GG"]/total_count*100+bases_transcriptome["GC"]/total_count*100),\
             (bases_transcriptome["CA"]/total_count*100+bases_transcriptome["CU"]/total_count*100+bases_transcriptome["CG"]/total_count*100+bases_transcriptome["CC"]/total_count*100)]
# print(dinucleotides)
"""cP libraries Pb-treated"""

total_count_1cP=bases["AA"][0]+bases["AU"][0]+bases["AG"][0]+bases["AC"][0]+\
            bases["UA"][0]+bases["UU"][0]+bases["UG"][0]+bases["UC"][0]+\
            bases["GA"][0]+bases["GU"][0]+bases["GG"][0]+bases["GC"][0]+\
            bases["CA"][0]+bases["CU"][0]+bases["CG"][0]+bases["CC"][0]
total_count_2cP=bases["AA"][1]+bases["AU"][1]+bases["AG"][1]+bases["AC"][1]+\
            bases["UA"][1]+bases["UU"][1]+bases["UG"][1]+bases["UC"][1]+\
            bases["GA"][1]+bases["GU"][1]+bases["GG"][1]+bases["GC"][1]+\
            bases["CA"][1]+bases["CU"][1]+bases["CG"][1]+bases["CC"][1]
base_A=np.mean([((bases["AA"][0]+bases["AU"][0]+bases["AG"][0]+bases["AC"][0])/total_count_1cP*100),((bases["AA"][1]+bases["AU"][1]+bases["AG"][1]+bases["AC"][1])/total_count_3cP*100)])
neighboring_base_A=np.mean([((bases["AA"][0]+bases["UA"][0]+bases["GA"][0]+bases["CA"][0])/total_count_1cP*100),((bases["AA"][1]+bases["UA"][1]+bases["GA"][1]+bases["CA"][1])/total_count_3cP*100)])

base_U=np.mean([((bases["UA"][0]+bases["UU"][0]+bases["UG"][0]+bases["UC"][0])/total_count_1cP*100),((bases["UA"][1]+bases["UU"][1]+bases["UG"][1]+bases["UC"][1])/total_count_3cP*100)])
neighboring_base_U=np.mean([((bases["AU"][0]+bases["UU"][0]+bases["GU"][0]+bases["CU"][0])/total_count_1cP*100),((bases["AU"][1]+bases["UU"][1]+bases["GU"][1]+bases["CU"][1])/total_count_3cP*100)])

base_G=np.mean([((bases["GA"][0]+bases["GU"][0]+bases["GG"][0]+bases["GC"][0])/total_count_1cP*100),((bases["GA"][1]+bases["GU"][1]+bases["GG"][1]+bases["GC"][1])/total_count_3cP*100)])
neighboring_base_G=np.mean([((bases["AG"][0]+bases["UG"][0]+bases["GG"][0]+bases["CG"][0])/total_count_1cP*100),((bases["AG"][1]+bases["UG"][1]+bases["GG"][1]+bases["CG"][1])/total_count_3cP*100)])

base_C=np.mean([((bases["CA"][0]+bases["CU"][0]+bases["CG"][0]+bases["CC"][0])/total_count_1cP*100),((bases["CA"][1]+bases["CU"][1]+bases["CG"][1]+bases["CC"][1])/total_count_3cP*100)])
neighboring_base_C=np.mean([((bases["AC"][0]+bases["UC"][0]+bases["GC"][0]+bases["CC"][0])/total_count_1cP*100),((bases["AC"][1]+bases["UC"][1]+bases["GC"][1]+bases["CC"][1])/total_count_3cP*100)])

nucleotides_cP_pb=[base_A, base_U, base_G, base_C]

neighboring_bases_1cP_pb=[bases["AA"][0]/total_count_1cP*100,
                         bases["AU"][0]/total_count_1cP*100,
                         bases["AG"][0]/total_count_1cP*100,
                         bases["AC"][0]/total_count_1cP*100,
                         bases["UA"][0]/total_count_1cP*100,
                         bases["UU"][0]/total_count_1cP*100,
                         bases["UG"][0]/total_count_1cP*100,
                         bases["UC"][0]/total_count_1cP*100,
                         bases["GA"][0]/total_count_1cP*100,
                         bases["GU"][0]/total_count_1cP*100,
                         bases["GG"][0]/total_count_1cP*100,
                         bases["GC"][0]/total_count_1cP*100,
                         bases["CA"][0]/total_count_1cP*100,
                         bases["CU"][0]/total_count_1cP*100,
                         bases["CG"][0]/total_count_1cP*100,
                         bases["CC"][0]/total_count_1cP*100]
neighboring_bases_2cP_pb=[bases["AA"][1]/total_count_2cP*100,
                         bases["AU"][1]/total_count_2cP*100,
                         bases["AG"][1]/total_count_2cP*100,
                         bases["AC"][1]/total_count_2cP*100,
                         bases["UA"][1]/total_count_2cP*100,
                         bases["UU"][1]/total_count_2cP*100,
                         bases["UG"][1]/total_count_2cP*100,
                         bases["UC"][1]/total_count_2cP*100,
                         bases["GA"][1]/total_count_2cP*100,
                         bases["GU"][1]/total_count_2cP*100,
                         bases["GG"][1]/total_count_2cP*100,
                         bases["GC"][1]/total_count_2cP*100,
                         bases["CA"][1]/total_count_2cP*100,
                         bases["CU"][1]/total_count_2cP*100,
                         bases["CG"][1]/total_count_2cP*100,
                         bases["CC"][1]/total_count_2cP*100]

dinucleotides_cP_pb=[np.mean([neighboring_bases_1cP_pb[0],neighboring_bases_2cP_pb[0]]),
np.mean([neighboring_bases_1cP_pb[1],neighboring_bases_2cP_pb[1]]),
np.mean([neighboring_bases_1cP_pb[2],neighboring_bases_2cP_pb[2]]),
np.mean([neighboring_bases_1cP_pb[3],neighboring_bases_2cP_pb[3]]),
np.mean([neighboring_bases_1cP_pb[4],neighboring_bases_2cP_pb[4]]),
np.mean([neighboring_bases_1cP_pb[5],neighboring_bases_2cP_pb[5]]),
np.mean([neighboring_bases_1cP_pb[6],neighboring_bases_2cP_pb[6]]),
np.mean([neighboring_bases_1cP_pb[7],neighboring_bases_2cP_pb[7]]),
np.mean([neighboring_bases_1cP_pb[8],neighboring_bases_2cP_pb[8]]),
np.mean([neighboring_bases_1cP_pb[9],neighboring_bases_2cP_pb[9]]),
np.mean([neighboring_bases_1cP_pb[10],neighboring_bases_2cP_pb[10]]),
np.mean([neighboring_bases_1cP_pb[11],neighboring_bases_2cP_pb[11]]),
np.mean([neighboring_bases_1cP_pb[12],neighboring_bases_2cP_pb[12]]),
np.mean([neighboring_bases_1cP_pb[13],neighboring_bases_2cP_pb[13]]),
np.mean([neighboring_bases_1cP_pb[14],neighboring_bases_2cP_pb[14]]),
np.mean([neighboring_bases_1cP_pb[15],neighboring_bases_2cP_pb[15]])]

dinucleotide_AA_cP_pb=bases_per_transcript["AA"][0]+bases_per_transcript["AA"][1]
dinucleotide_AU_cP_pb=bases_per_transcript["AU"][0]+bases_per_transcript["AU"][1]
dinucleotide_AG_cP_pb=bases_per_transcript["AG"][0]+bases_per_transcript["AG"][1]
dinucleotide_AC_cP_pb=bases_per_transcript["AC"][0]+bases_per_transcript["AC"][1]
dinucleotide_UA_cP_pb=bases_per_transcript["UA"][0]+bases_per_transcript["UA"][1]
dinucleotide_UU_cP_pb=bases_per_transcript["UU"][0]+bases_per_transcript["UU"][1]
dinucleotide_UG_cP_pb=bases_per_transcript["UG"][0]+bases_per_transcript["UG"][1]
dinucleotide_UC_cP_pb=bases_per_transcript["UC"][0]+bases_per_transcript["UC"][1]
dinucleotide_GA_cP_pb=bases_per_transcript["GA"][0]+bases_per_transcript["GA"][1]
dinucleotide_GU_cP_pb=bases_per_transcript["GU"][0]+bases_per_transcript["GU"][1]
dinucleotide_GG_cP_pb=bases_per_transcript["GG"][0]+bases_per_transcript["GG"][1]
dinucleotide_GC_cP_pb=bases_per_transcript["GC"][0]+bases_per_transcript["GC"][1]
dinucleotide_CA_cP_pb=bases_per_transcript["CA"][0]+bases_per_transcript["CA"][1]
dinucleotide_CU_cP_pb=bases_per_transcript["CU"][0]+bases_per_transcript["CU"][1]
dinucleotide_CG_cP_pb=bases_per_transcript["CG"][0]+bases_per_transcript["CG"][1]
dinucleotide_CC_cP_pb=bases_per_transcript["CC"][0]+bases_per_transcript["CC"][1]

# xcP1 = ["2',3'-cP", "next to"]
# xcP2 = ["end", "2',3'-cP end"]
# xcP = [f"{x1}\n{x2}" for x1, x2, in zip(xcP1,xcP2)]
#
# y_A = np.array([base_A, neighboring_base_A])
# y_U = np.array([base_U, neighboring_base_U])
# y_G = np.array([base_G, neighboring_base_G])
# y_C = np.array([base_C, neighboring_base_C])
# barWidth = 0.7
#
# fig, (ax1, ax2) = plt.subplots(1, 2,  sharey=True, figsize=(5.7,5))#sharex=True,
# # fig = plt.figure(figsize=(4,5))
# # ax = fig.gca()
# ax1.bar(xcP, y_C, width=barWidth, bottom=y_G + y_A + y_U, color="#67a9cf", label = "C")
# ax1.bar(xcP, y_U, width=barWidth, bottom=y_G + y_A, color="#2166ac", label = "U")
# ax1.bar(xcP, y_A, width=barWidth, bottom=y_G, color="#f4a582", label = "A")
# ax1.bar(xcP, y_G, width=barWidth, color="#ca0020", label = "G")
# ax1.tick_params(axis='x', which='major', labelsize=14)
# ax1.tick_params(axis='x', which='minor', labelsize=14)
# ax1.set_ylabel("nucleotides at cleavage sites [%]")
# ax1.set_title("2'3'-cP")
# # plt.ylim(top=102)
# # plt.xticks(rotation = 45)
# # plt.legend()
# # plt.title("2'3'-cP")
# # plt.tight_layout()
# # plt.savefig("./Nucleotide_compositon_cP_pb_barplot.pdf")
# # plt.close()

"""OH libraries Pb-treated"""
total_count_1OH=bases["AA"][3]+bases["AU"][3]+bases["AG"][3]+bases["AC"][3]+\
            bases["UA"][3]+bases["UU"][3]+bases["UG"][3]+bases["UC"][3]+\
            bases["GA"][3]+bases["GU"][3]+bases["GG"][3]+bases["GC"][3]+\
            bases["CA"][3]+bases["CU"][3]+bases["CG"][3]+bases["CC"][3]
total_count_2OH=bases["AA"][4]+bases["AU"][4]+bases["AG"][4]+bases["AC"][4]+\
            bases["UA"][4]+bases["UU"][4]+bases["UG"][4]+bases["UC"][4]+\
            bases["GA"][4]+bases["GU"][4]+bases["GG"][4]+bases["GC"][4]+\
            bases["CA"][4]+bases["CU"][4]+bases["CG"][4]+bases["CC"][4]
base_A=np.mean([((bases["AA"][3]+bases["AU"][3]+bases["AG"][3]+bases["AC"][3])/total_count_1OH*100),((bases["AA"][4]+bases["AU"][4]+bases["AG"][4]+bases["AC"][4])/total_count_2OH*100)])
neighboring_base_A=np.mean([((bases["AA"][3]+bases["UA"][3]+bases["GA"][3]+bases["CA"][3])/total_count_1OH*100),((bases["AA"][4]+bases["UA"][4]+bases["GA"][4]+bases["CA"][4])/total_count_2OH*100)])

base_U=np.mean([((bases["UA"][3]+bases["UU"][3]+bases["UG"][3]+bases["UC"][3])/total_count_1OH*100),((bases["UA"][4]+bases["UU"][4]+bases["UG"][4]+bases["UC"][4])/total_count_2OH*100)])
neighboring_base_U=np.mean([((bases["AU"][3]+bases["UU"][3]+bases["GU"][3]+bases["CU"][3])/total_count_1OH*100),((bases["AU"][4]+bases["UU"][4]+bases["GU"][4]+bases["CU"][4])/total_count_2OH*100)])

base_G=np.mean([((bases["GA"][3]+bases["GU"][3]+bases["GG"][3]+bases["GC"][3])/total_count_1OH*100),((bases["GA"][4]+bases["GU"][4]+bases["GG"][4]+bases["GC"][4])/total_count_2OH*100)])
neighboring_base_G=np.mean([((bases["AG"][3]+bases["UG"][3]+bases["GG"][3]+bases["CG"][3])/total_count_1OH*100),((bases["AG"][4]+bases["UG"][4]+bases["GG"][4]+bases["CG"][4])/total_count_2OH*100)])

base_C=np.mean([((bases["CA"][3]+bases["CU"][3]+bases["CG"][3]+bases["CC"][3])/total_count_1OH*100),((bases["CA"][4]+bases["CU"][4]+bases["CG"][4]+bases["CC"][4])/total_count_2OH*100)])
neighboring_base_C=np.mean([((bases["AC"][3]+bases["UC"][3]+bases["GC"][3]+bases["CC"][3])/total_count_1OH*100),((bases["AC"][4]+bases["UC"][4]+bases["GC"][4]+bases["CC"][4])/total_count_2OH*100)])

nucleotides_OH_pb=[base_A, base_U, base_G, base_C]

neighboring_bases_1OH_pb=[bases["AA"][3]/total_count_1OH*100,
                         bases["AU"][3]/total_count_1OH*100,
                         bases["AG"][3]/total_count_1OH*100,
                         bases["AC"][3]/total_count_1OH*100,
                         bases["UA"][3]/total_count_1OH*100,
                         bases["UU"][3]/total_count_1OH*100,
                         bases["UG"][3]/total_count_1OH*100,
                         bases["UC"][3]/total_count_1OH*100,
                         bases["GA"][3]/total_count_1OH*100,
                         bases["GU"][3]/total_count_1OH*100,
                         bases["GG"][3]/total_count_1OH*100,
                         bases["GC"][3]/total_count_1OH*100,
                         bases["CA"][3]/total_count_1OH*100,
                         bases["CU"][3]/total_count_1OH*100,
                         bases["CG"][3]/total_count_1OH*100,
                         bases["CC"][3]/total_count_1OH*100]
neighboring_bases_2OH_pb=[bases["AA"][4]/total_count_2OH*100,
                         bases["AU"][4]/total_count_2OH*100,
                         bases["AG"][4]/total_count_2OH*100,
                         bases["AC"][4]/total_count_2OH*100,
                         bases["UA"][4]/total_count_2OH*100,
                         bases["UU"][4]/total_count_2OH*100,
                         bases["UG"][4]/total_count_2OH*100,
                         bases["UC"][4]/total_count_2OH*100,
                         bases["GA"][4]/total_count_2OH*100,
                         bases["GU"][4]/total_count_2OH*100,
                         bases["GG"][4]/total_count_2OH*100,
                         bases["GC"][4]/total_count_2OH*100,
                         bases["CA"][4]/total_count_2OH*100,
                         bases["CU"][4]/total_count_2OH*100,
                         bases["CG"][4]/total_count_2OH*100,
                         bases["CC"][4]/total_count_2OH*100]

dinucleotides_OH_pb=[np.mean([neighboring_bases_1OH_pb[0],neighboring_bases_2OH_pb[0]]),
np.mean([neighboring_bases_1OH_pb[1],neighboring_bases_2OH_pb[1]]),
np.mean([neighboring_bases_1OH_pb[2],neighboring_bases_2OH_pb[2]]),
np.mean([neighboring_bases_1OH_pb[3],neighboring_bases_2OH_pb[3]]),
np.mean([neighboring_bases_1OH_pb[4],neighboring_bases_2OH_pb[4]]),
np.mean([neighboring_bases_1OH_pb[5],neighboring_bases_2OH_pb[5]]),
np.mean([neighboring_bases_1OH_pb[6],neighboring_bases_2OH_pb[6]]),
np.mean([neighboring_bases_1OH_pb[7],neighboring_bases_2OH_pb[7]]),
np.mean([neighboring_bases_1OH_pb[8],neighboring_bases_2OH_pb[8]]),
np.mean([neighboring_bases_1OH_pb[9],neighboring_bases_2OH_pb[9]]),
np.mean([neighboring_bases_1OH_pb[10],neighboring_bases_2OH_pb[10]]),
np.mean([neighboring_bases_1OH_pb[11],neighboring_bases_2OH_pb[11]]),
np.mean([neighboring_bases_1OH_pb[12],neighboring_bases_2OH_pb[12]]),
np.mean([neighboring_bases_1OH_pb[13],neighboring_bases_2OH_pb[13]]),
np.mean([neighboring_bases_1OH_pb[14],neighboring_bases_2OH_pb[14]]),
np.mean([neighboring_bases_1OH_pb[15],neighboring_bases_2OH_pb[15]])]

dinucleotide_AA_OH_pb=bases_per_transcript["AA"][3]+bases_per_transcript["AA"][4]
dinucleotide_AU_OH_pb=bases_per_transcript["AU"][3]+bases_per_transcript["AU"][4]
dinucleotide_AG_OH_pb=bases_per_transcript["AG"][3]+bases_per_transcript["AG"][4]
dinucleotide_AC_OH_pb=bases_per_transcript["AC"][3]+bases_per_transcript["AC"][4]
dinucleotide_UA_OH_pb=bases_per_transcript["UA"][3]+bases_per_transcript["UA"][4]
dinucleotide_UU_OH_pb=bases_per_transcript["UU"][3]+bases_per_transcript["UU"][4]
dinucleotide_UG_OH_pb=bases_per_transcript["UG"][3]+bases_per_transcript["UG"][4]
dinucleotide_UC_OH_pb=bases_per_transcript["UC"][3]+bases_per_transcript["UC"][4]
dinucleotide_GA_OH_pb=bases_per_transcript["GA"][3]+bases_per_transcript["GA"][4]
dinucleotide_GU_OH_pb=bases_per_transcript["GU"][3]+bases_per_transcript["GU"][4]
dinucleotide_GG_OH_pb=bases_per_transcript["GG"][3]+bases_per_transcript["GG"][4]
dinucleotide_GC_OH_pb=bases_per_transcript["GC"][3]+bases_per_transcript["GC"][4]
dinucleotide_CA_OH_pb=bases_per_transcript["CA"][3]+bases_per_transcript["CA"][4]
dinucleotide_CU_OH_pb=bases_per_transcript["CU"][3]+bases_per_transcript["CU"][4]
dinucleotide_CG_OH_pb=bases_per_transcript["CG"][3]+bases_per_transcript["CG"][4]
dinucleotide_CC_OH_pb=bases_per_transcript["CC"][3]+bases_per_transcript["CC"][4]

# xOH1 = ["next to", "5'-OH"]
# xOH2 = ["5'-OH end", "end"]
# xOH = [f"{x1}\n{x2}" for x1, x2, in zip(xOH1,xOH2)]
# y_A = np.array([base_A, neighboring_base_A])
# y_U = np.array([base_U, neighboring_base_U])
# y_G = np.array([base_G, neighboring_base_G])
# y_C = np.array([base_C, neighboring_base_C])
# barWidth = 0.7
#
# # fig = plt.figure(figsize=(4,5))
# # ax = fig.gca()
# ax2.bar(xOH, y_C, width=barWidth, bottom=y_G + y_A + y_U, color="#67a9cf", label = "C")
# ax2.bar(xOH, y_U, width=barWidth, bottom=y_G + y_A, color="#2166ac", label = "U")
# ax2.bar(xOH, y_A, width=barWidth, bottom=y_G, color="#f4a582", label = "A")
# ax2.bar(xOH, y_G, width=barWidth, color="#ca0020", label = "G")
# # plt.ylabel("nucleotides at cleavage sites [%]")
# # plt.xticks(rotation = 45)
# ax2.tick_params(axis='x', which='major', labelsize=14)
# ax2.tick_params(axis='x', which='minor', labelsize=14)
# plt.ylim(top=102)
# plt.legend()
# plt.title("5'-OH")
# plt.tight_layout()
# plt.savefig("./Nucleotide_compositon_both_pb_barplot.pdf")
# plt.close()


"""cP libraries H2O"""
total_count_3cP=bases["AA"][2]+bases["AU"][2]+bases["AG"][2]+bases["AC"][2]+\
            bases["UA"][2]+bases["UU"][2]+bases["UG"][2]+bases["UC"][2]+\
            bases["GA"][2]+bases["GU"][2]+bases["GG"][2]+bases["GC"][2]+\
            bases["CA"][2]+bases["CU"][2]+bases["CG"][2]+bases["CC"][2]

base_A=(bases["AA"][2]+bases["AU"][2]+bases["AG"][2]+bases["AC"][2])/total_count_3cP*100
neighboring_base_A=(bases["AA"][2]+bases["UA"][2]+bases["GA"][2]+bases["CA"][2])/total_count_3cP*100

base_U=(bases["UA"][2]+bases["UU"][2]+bases["UG"][2]+bases["UC"][2])/total_count_3cP*100
neighboring_base_U=(bases["AU"][2]+bases["UU"][2]+bases["GU"][2]+bases["CU"][2])/total_count_3cP*100

base_G=(bases["GA"][2]+bases["GU"][2]+bases["GG"][2]+bases["GC"][2])/total_count_3cP*100
neighboring_base_G=(bases["AG"][2]+bases["UG"][2]+bases["GG"][2]+bases["CG"][2])/total_count_3cP*100

base_C=(bases["CA"][2]+bases["CU"][2]+bases["CG"][2]+bases["CC"][2])/total_count_3cP*100
neighboring_base_C=(bases["AC"][2]+bases["UC"][2]+bases["GC"][2]+bases["CC"][2])/total_count_3cP*100

nucleotides_cP_h2o=[base_A, base_U, base_G, base_C]
dinucleotides_cP_h2o=[bases["AA"][2]/total_count_3cP*100,
                         bases["AU"][2]/total_count_3cP*100,
                         bases["AG"][2]/total_count_3cP*100,
                         bases["AC"][2]/total_count_3cP*100,
                         bases["UA"][2]/total_count_3cP*100,
                         bases["UU"][2]/total_count_3cP*100,
                         bases["UG"][2]/total_count_3cP*100,
                         bases["UC"][2]/total_count_3cP*100,
                         bases["GA"][2]/total_count_3cP*100,
                         bases["GU"][2]/total_count_3cP*100,
                         bases["GG"][2]/total_count_3cP*100,
                         bases["GC"][2]/total_count_3cP*100,
                         bases["CA"][2]/total_count_3cP*100,
                         bases["CU"][2]/total_count_3cP*100,
                         bases["CG"][2]/total_count_3cP*100,
                         bases["CC"][2]/total_count_3cP*100]

dinucleotide_AA_cP_h2o=bases_per_transcript["AA"][2]
dinucleotide_AU_cP_h2o=bases_per_transcript["AU"][2]
dinucleotide_AG_cP_h2o=bases_per_transcript["AG"][2]
dinucleotide_AC_cP_h2o=bases_per_transcript["AC"][2]
dinucleotide_UA_cP_h2o=bases_per_transcript["UA"][2]
dinucleotide_UU_cP_h2o=bases_per_transcript["UU"][2]
dinucleotide_UG_cP_h2o=bases_per_transcript["UG"][2]
dinucleotide_UC_cP_h2o=bases_per_transcript["UC"][2]
dinucleotide_GA_cP_h2o=bases_per_transcript["GA"][2]
dinucleotide_GU_cP_h2o=bases_per_transcript["GU"][2]
dinucleotide_GG_cP_h2o=bases_per_transcript["GG"][2]
dinucleotide_GC_cP_h2o=bases_per_transcript["GC"][2]
dinucleotide_CA_cP_h2o=bases_per_transcript["CA"][2]
dinucleotide_CU_cP_h2o=bases_per_transcript["CU"][2]
dinucleotide_CG_cP_h2o=bases_per_transcript["CG"][2]
dinucleotide_CC_cP_h2o=bases_per_transcript["CC"][2]

# y_A = np.array([base_A, neighboring_base_A])
# y_U = np.array([base_U, neighboring_base_U])
# y_G = np.array([base_G, neighboring_base_G])
# y_C = np.array([base_C, neighboring_base_C])
#
# fig = plt.figure(figsize=(3,6))
# ax = fig.gca()
# ax.bar(xcP, y_C, width=barWidth, bottom=y_G + y_A + y_U, color='y', label = "C")
# ax.bar(xcP, y_U, width=barWidth, bottom=y_G + y_A, color='lightblue', label = "U")
# ax.bar(xcP, y_A, width=barWidth, bottom=y_G, color='g', label = "A")
# ax.bar(xcP, y_G, width=barWidth, color='lightgrey', label = "G")
# plt.ylabel("Nucleotide composition (%)")
# # plt.xticks(rotation = 45)
# plt.legend()
# plt.title("cP libraries H2O")
# plt.tight_layout()
# plt.savefig("./Nucleotide_compositon_cP_h2o_barplot.pdf")
# plt.close()

"""OH libraries H2O"""
total_count_3OH=bases["AA"][5]+bases["AU"][5]+bases["AG"][5]+bases["AC"][5]+\
            bases["UA"][5]+bases["UU"][5]+bases["UG"][5]+bases["UC"][5]+\
            bases["GA"][5]+bases["GU"][5]+bases["GG"][5]+bases["GC"][5]+\
            bases["CA"][5]+bases["CU"][5]+bases["CG"][5]+bases["CC"][5]

base_A=(bases["AA"][5]+bases["AU"][5]+bases["AG"][5]+bases["AC"][5])/total_count_3OH*100
neighboring_base_A=(bases["AA"][5]+bases["UA"][5]+bases["GA"][5]+bases["CA"][5])/total_count_3OH*100

base_U=(bases["UA"][5]+bases["UU"][5]+bases["UG"][5]+bases["UC"][5])/total_count_3OH*100
neighboring_base_U=(bases["AU"][5]+bases["UU"][5]+bases["GU"][5]+bases["CU"][5])/total_count_3OH*100

base_G=(bases["GA"][5]+bases["GU"][5]+bases["GG"][5]+bases["GC"][5])/total_count_3OH*100
neighboring_base_G=(bases["AG"][5]+bases["UG"][5]+bases["GG"][5]+bases["CG"][5])/total_count_3OH*100

base_C=(bases["CA"][5]+bases["CU"][5]+bases["CG"][5]+bases["CC"][5])/total_count_3OH*100
neighboring_base_C=(bases["AC"][5]+bases["UC"][5]+bases["GC"][5]+bases["CC"][5])/total_count_3OH*100

nucleotides_OH_h2o=[base_A, base_U, base_G, base_C]

dinucleotides_OH_h2o=[bases["AA"][5]/total_count_3OH*100,
                         bases["AU"][5]/total_count_3OH*100,
                         bases["AG"][5]/total_count_3OH*100,
                         bases["AC"][5]/total_count_3OH*100,
                         bases["UA"][5]/total_count_3OH*100,
                         bases["UU"][5]/total_count_3OH*100,
                         bases["UG"][5]/total_count_3OH*100,
                         bases["UC"][5]/total_count_3OH*100,
                         bases["GA"][5]/total_count_3OH*100,
                         bases["GU"][5]/total_count_3OH*100,
                         bases["GG"][5]/total_count_3OH*100,
                         bases["GC"][5]/total_count_3OH*100,
                         bases["CA"][5]/total_count_3OH*100,
                         bases["CU"][5]/total_count_3OH*100,
                         bases["CG"][5]/total_count_3OH*100,
                         bases["CC"][5]/total_count_3OH*100]

dinucleotide_AA_OH_h2o=bases_per_transcript["AA"][5]
dinucleotide_AU_OH_h2o=bases_per_transcript["AU"][5]
dinucleotide_AG_OH_h2o=bases_per_transcript["AG"][5]
dinucleotide_AC_OH_h2o=bases_per_transcript["AC"][5]
dinucleotide_UA_OH_h2o=bases_per_transcript["UA"][5]
dinucleotide_UU_OH_h2o=bases_per_transcript["UU"][5]
dinucleotide_UG_OH_h2o=bases_per_transcript["UG"][5]
dinucleotide_UC_OH_h2o=bases_per_transcript["UC"][5]
dinucleotide_GA_OH_h2o=bases_per_transcript["GA"][5]
dinucleotide_GU_OH_h2o=bases_per_transcript["GU"][5]
dinucleotide_GG_OH_h2o=bases_per_transcript["GG"][5]
dinucleotide_GC_OH_h2o=bases_per_transcript["GC"][5]
dinucleotide_CA_OH_h2o=bases_per_transcript["CA"][5]
dinucleotide_CU_OH_h2o=bases_per_transcript["CU"][5]
dinucleotide_CG_OH_h2o=bases_per_transcript["CG"][5]
dinucleotide_CC_OH_h2o=bases_per_transcript["CC"][5]

# y_A = np.array([base_A, neighboring_base_A])
# y_U = np.array([base_U, neighboring_base_U])
# y_G = np.array([base_G, neighboring_base_G])
# y_C = np.array([base_C, neighboring_base_C])
#
# fig = plt.figure(figsize=(3,6))
# ax = fig.gca()
# ax.bar(xOH, y_C, width=barWidth, bottom=y_G + y_A + y_U, color='y', label = "C")
# ax.bar(xOH, y_U, width=barWidth, bottom=y_G + y_A, color='lightblue', label = "U")
# ax.bar(xOH, y_A, width=barWidth, bottom=y_G, color='g', label = "A")
# ax.bar(xOH, y_G, width=barWidth, color='lightgrey', label = "G")
# plt.ylabel("Nucleotide composition (%)")
# # plt.xticks(rotation = 45)
# plt.legend()
# plt.title("OH libraries H2O")
# plt.tight_layout()
# plt.savefig("./Nucleotide_compositon_OH_h2o_barplot.pdf")
# plt.close()



"""plot occurrence in transcriptome"""
fig = plt.figure(figsize=(10,4))
ax = fig.gca()
ax.bar(["AA","AU","AG","AC","UA","UU","UG","UC","GA","GU","GG","GC","CA","CU","CG","CC"], dinucleotides, color="#0571b0")
#ax.axhline(6.25, linestyle='--', color="green")
plt.ylabel("occurrence \n in transcriptome [%]")
plt.xlim(left=-0.5, right=15.5)
plt.savefig("./Dinucleotides_transcriptome.pdf")
plt.close()


# fig = plt.figure(figsize=(5,4))
# ax = fig.gca()
# ax.bar(["A","U","G","C"], nucleotides, color="#0571b0")
# # ax.axhline(25, linestyle='--', color="green")
# plt.ylabel("occurrence \n in transcriptome [%]")
# plt.ylim(bottom=0)
# plt.tight_layout()
# plt.savefig("./Nucleotides_transcriptome.pdf")
# plt.close()


# barWidth = 0.5
# fig = plt.figure(figsize=(3,5))
# ax = fig.gca()
# ax.bar(1, nucleotides[3], width=barWidth, bottom=nucleotides[2] + nucleotides[0] + nucleotides[1], color="#67a9cf", label = "C")
# ax.bar(1, nucleotides[1], width=barWidth, bottom=nucleotides[2] + nucleotides[0], color="#2166ac", label = "U")
# ax.bar(1, nucleotides[0], width=barWidth, bottom=nucleotides[2], color="#f4a582", label = "A")
# ax.bar(1, nucleotides[2], width=barWidth, color="#ca0020", label = "G")
# ax.tick_params(axis='x', which='both', bottom=False, top=False,labelbottom=False)
# plt.ylabel("occurrence in transcriptome [%]")
# plt.ylim(top=102)
# plt.legend()
# plt.tight_layout()
# plt.savefig("./Nucleotides_transcriptome_stacked_barplot.pdf")
# plt.close()

"""plot occurrence at cleavage site with boxplots"""
width=0.4
labels=["AA","AU","AG","AC","UA","UU","UG","UC","GA","GU","GG","GC","CA","CU","CG","CC"]
x=np.arange(len(labels))

fig = plt.figure(figsize=(10,4))
ax = fig.gca()
ax.bar(x-width/2, dinucleotides_cP_pb, width,  label="2'3'-cP", color="#f4a582")
ax.bar(x+width/2, dinucleotides_OH_pb, width,  label="5'-OH", color="#0571b0")
medianprops = dict(color="#ca0020")
ax.boxplot([dinucleotide_AA_cP_pb,
            dinucleotide_AU_cP_pb,
            dinucleotide_AG_cP_pb,
            dinucleotide_AC_cP_pb,
            dinucleotide_UA_cP_pb,
            dinucleotide_UU_cP_pb,
            dinucleotide_UG_cP_pb,
            dinucleotide_UC_cP_pb,
            dinucleotide_GA_cP_pb,
            dinucleotide_GU_cP_pb,
            dinucleotide_GG_cP_pb,
            dinucleotide_GC_cP_pb,
            dinucleotide_CA_cP_pb,
            dinucleotide_CU_cP_pb,
            dinucleotide_CG_cP_pb,
            dinucleotide_CC_cP_pb], positions=x-width/2, widths=0.3, showfliers=False, medianprops=medianprops,whis=1.5)
ax.boxplot([dinucleotide_AA_OH_pb,
            dinucleotide_AU_OH_pb,
            dinucleotide_AG_OH_pb,
            dinucleotide_AC_OH_pb,
            dinucleotide_UA_OH_pb,
            dinucleotide_UU_OH_pb,
            dinucleotide_UG_OH_pb,
            dinucleotide_UC_OH_pb,
            dinucleotide_GA_OH_pb,
            dinucleotide_GU_OH_pb,
            dinucleotide_GG_OH_pb,
            dinucleotide_GC_OH_pb,
            dinucleotide_CA_OH_pb,
            dinucleotide_CU_OH_pb,
            dinucleotide_CG_OH_pb,
            dinucleotide_CC_OH_pb], positions=x+width/2, widths=0.3, showfliers=False, medianprops=medianprops,whis=1.5)
plt.xticks(np.arange(0,16), labels)
plt.xlim(left=-0.5, right=15.5)
plt.ylabel("cleaved dinucleotide [%]")
plt.legend()
plt.savefig("./Dinucleotides_cleavage_sites_Pb_boxplots.pdf", bbox_inches='tight')
plt.close()


"""plot occurrence at cleavage site"""
width=0.4
labels=["AA","AU","AG","AC","UA","UU","UG","UC","GA","GU","GG","GC","CA","CU","CG","CC"]
x=np.arange(len(labels))

fig = plt.figure(figsize=(10,4))
ax = fig.gca()
ax.bar(x-width/2, dinucleotides_cP_pb, width,  label="2'3'-cP", color="#f4a582")
ax.bar(x+width/2, dinucleotides_OH_pb, width,  label="5'-OH", color="#0571b0")
plt.xticks(np.arange(0,16), labels)
plt.xlim(left=-0.5, right=15.5)
plt.ylabel("cleaved dinucleotide [%]")
plt.legend()
plt.savefig("./Dinucleotides_cleavage_sites_Pb.pdf", bbox_inches='tight')
plt.close()



# fig = plt.figure(figsize=(10,4))
# ax = fig.gca()
# ax.bar(x-width/2, dinucleotides_cP_h2o, width,  label="2'3'-cP $\mathregular{Pb^{2\!+}}\!(-)$", color="#f4a582")
# ax.bar(x+width/2, dinucleotides_OH_h2o, width,  label="5'-OH $\mathregular{Pb^{2\!+}}\!(-)$", color="#0571b0")
# plt.xticks(np.arange(0,16), labels)
# plt.xlim(left=-0.5, right=15.5)
# plt.ylabel("cleaved dinucleotide [%]")
# plt.legend()
# plt.savefig("./Dinucleotides_cleavage_sites_H2O.pdf", bbox_inches='tight')
# plt.close()


"""plot occurrence at cleavage site with boxplots H2O"""

fig = plt.figure(figsize=(10,4))
ax = fig.gca()
ax.bar(x-width/2, dinucleotides_cP_h2o, width,  label="2'3'-cP $\mathregular{Pb^{2\!+}}\!(-)$", color="#f4a582")
ax.bar(x+width/2, dinucleotides_OH_h2o, width,  label="5'-OH $\mathregular{Pb^{2\!+}}\!(-)$", color="#0571b0")
medianprops = dict(color="#ca0020")
ax.boxplot([dinucleotide_AA_cP_h2o,
            dinucleotide_AU_cP_h2o,
            dinucleotide_AG_cP_h2o,
            dinucleotide_AC_cP_h2o,
            dinucleotide_UA_cP_h2o,
            dinucleotide_UU_cP_h2o,
            dinucleotide_UG_cP_h2o,
            dinucleotide_UC_cP_h2o,
            dinucleotide_GA_cP_h2o,
            dinucleotide_GU_cP_h2o,
            dinucleotide_GG_cP_h2o,
            dinucleotide_GC_cP_h2o,
            dinucleotide_CA_cP_h2o,
            dinucleotide_CU_cP_h2o,
            dinucleotide_CG_cP_h2o,
            dinucleotide_CC_cP_h2o], positions=x-width/2, widths=0.3, showfliers=False, medianprops=medianprops,whis=1.5)#, whis=[10, 90])
ax.boxplot([dinucleotide_AA_OH_h2o,
            dinucleotide_AU_OH_h2o,
            dinucleotide_AG_OH_h2o,
            dinucleotide_AC_OH_h2o,
            dinucleotide_UA_OH_h2o,
            dinucleotide_UU_OH_h2o,
            dinucleotide_UG_OH_h2o,
            dinucleotide_UC_OH_h2o,
            dinucleotide_GA_OH_h2o,
            dinucleotide_GU_OH_h2o,
            dinucleotide_GG_OH_h2o,
            dinucleotide_GC_OH_h2o,
            dinucleotide_CA_OH_h2o,
            dinucleotide_CU_OH_h2o,
            dinucleotide_CG_OH_h2o,
            dinucleotide_CC_OH_h2o], positions=x+width/2, widths=0.3, showfliers=False, medianprops=medianprops,whis=1.5)#,whis=[10, 90])
plt.xticks(np.arange(0,16), labels)
plt.xlim(left=-0.5, right=15.5)
plt.ylabel("cleaved dinucleotide [%]")
plt.legend()
plt.savefig("./Dinucleotides_cleavage_sites_H2O_with_boxplots.pdf", bbox_inches='tight')
plt.close()
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import sys
import statistics
import random
plt.rcParams.update({'font.size':19})
plt.rcParams.update({'axes.titlesize':19})
plt.rc('legend', fontsize=18)
plt.rc('xtick', labelsize=19)
plt.rc('ytick', labelsize=19)

signal_file = sys.argv[1]
benchmark_gff_file = sys.argv[2]
output_file_name = sys.argv[3]
output_folder = sys.argv[4]


"""read in signal_file"""
signal = {}
with open(signal_file) as file:
    for line in file:
        line=line.rstrip()
        columns = line.split("\t")
        base, neighboring_base, transcript_id, nc, position, structure, library, library_type = columns[3], columns[4], columns[5], columns[7], int(columns[8]), columns[9], columns[10], columns[11]
        if transcript_id not in signal:
            signal[transcript_id] = {}
        if position not in signal[transcript_id]:
            signal[transcript_id][position]={}
        if "structure" not in signal[transcript_id][position]:
            signal[transcript_id][position]["structure"] = structure
        if "CA" not in signal[transcript_id][position]:
            if base == "C" and neighboring_base == "A":
                signal[transcript_id][position]["CA"] = True
            else:
                signal[transcript_id][position]["CA"] = False
        signal[transcript_id][position][library]=nc

"""read in benchmark set"""
benchmark_transcripts = []

with open(benchmark_gff_file, 'r') as file:
    for line in file:
        line = line.rstrip()
        columns = line.split("\t")
        transcript_id = columns[8]
        benchmark_transcripts.append(transcript_id)


"""data for the reference set"""
ref_zero_zero_all = 0
ref_zero_zero_unpaired = 0
ref_cP_zero_all = 0
ref_cP_zero_unpaired = 0
ref_OH_zero_all = 0
ref_OH_zero_unpaired = 0
joint_occ_reference_set = []
joint_occ_reference_set_CA = []
OH_joint_occ_reference_set = []
OH_joint_occ_reference_set_CA = []
cP_joint_occ_reference_set = []
cP_joint_occ_reference_set_CA = []

h2o_ref_zero_zero_all = 0
h2o_ref_zero_zero_unpaired = 0
h2o_ref_cP_zero_all = 0
h2o_ref_cP_zero_unpaired = 0
h2o_ref_OH_zero_all = 0
h2o_ref_OH_zero_unpaired = 0
h2o_joint_occ_reference_set = []
h2o_joint_occ_reference_set_CA = []
h2o_OH_joint_occ_reference_set = []
h2o_OH_joint_occ_reference_set_CA = []
h2o_cP_joint_occ_reference_set = []
h2o_cP_joint_occ_reference_set_CA = []


pos_unpaired_counter = 0
pos_counter = 0
CA_counter = 0
for transcript_id in signal:
    if transcript_id in benchmark_transcripts:   #save the data for the reference set
        for pos in signal[transcript_id]:
            pos_counter += 1
            CA = False
            if signal[transcript_id][pos]["CA"]:
                CA = True
                CA_counter += 1

            signal_1_OH=signal[transcript_id][pos]["1-OH"]
            signal_2_OH=signal[transcript_id][pos]["2-OH"]
            signal_3_OH = signal[transcript_id][pos]["3-OH"]        # h2o sample
            signal_1_cP=signal[transcript_id][pos]["1-cP"]
            signal_2_cP = signal[transcript_id][pos]["2-cP"]        
            signal_3_cP=signal[transcript_id][pos]["3-cP"]          # h2o sample
            structure=signal[transcript_id][pos]["structure"]
            p_unpaired=None
            if structure == "0":
                p_unpaired=1.0
                pos_unpaired_counter += 1
            elif structure == "1":
                p_unpaired=0.0

            """Pb treated libraries"""
            if all(elem != "NA" for elem in [signal_1_OH,signal_2_OH,signal_1_cP,signal_2_cP]):         # valid signal in all libraries
                if not CA:
                    joint_occ_reference_set.append([statistics.mean([float(signal_1_cP), float(signal_2_cP)]), statistics.mean([float(signal_1_OH), float(signal_2_OH)]), p_unpaired])
                else:
                    joint_occ_reference_set_CA.append([statistics.mean([float(signal_1_cP), float(signal_2_cP)]),
                                                    statistics.mean([float(signal_1_OH), float(signal_2_OH)]),
                                                    p_unpaired])

                if all(float(elem) == 0.0 for elem in [signal_1_OH, signal_2_OH, signal_1_cP, signal_2_cP]):
                    ref_zero_zero_all += 1
                    if structure == "0":
                        ref_zero_zero_unpaired += 1

            if all(elem != "NA" for elem in [signal_1_cP,signal_2_cP]):                                 # valid signal in all cP libraries
                if not CA:
                    cP_joint_occ_reference_set.append([statistics.mean([float(signal_1_cP), float(signal_2_cP)]), p_unpaired])

                else:
                    cP_joint_occ_reference_set_CA.append(
                        [statistics.mean([float(signal_1_cP), float(signal_2_cP)]), p_unpaired])

                if all(float(elem) == 0.0 for elem in [signal_1_cP, signal_2_cP]):
                    ref_cP_zero_all += 1
                    if structure == "0":
                        ref_cP_zero_unpaired += 1


            if all(elem != "NA" for elem in [signal_1_OH,signal_2_OH]):                                 # valid signal in all OH libraries
                if not CA:
                    OH_joint_occ_reference_set.append([statistics.mean([float(signal_1_OH), float(signal_2_OH)]), p_unpaired])
                else:
                    OH_joint_occ_reference_set_CA.append(
                        [statistics.mean([float(signal_1_OH), float(signal_2_OH)]), p_unpaired])

                if all(float(elem) == 0.0 for elem in [signal_1_OH, signal_2_OH]):
                    ref_OH_zero_all += 1
                    if structure == "0":
                        ref_OH_zero_unpaired += 1

            """H2O libraries"""
            if all(elem != "NA" for elem in [signal_3_OH, signal_3_cP]):         # valid signal in both libraries
                if not CA:
                    h2o_joint_occ_reference_set.append([float(signal_3_cP), float(signal_3_OH), p_unpaired])
                else:
                    h2o_joint_occ_reference_set_CA.append([float(signal_3_cP), float(signal_3_OH), p_unpaired])

                if all(float(elem) == 0.0 for elem in [signal_3_OH, signal_3_cP]):
                    h2o_ref_zero_zero_all += 1
                    if structure == "0":
                        h2o_ref_zero_zero_unpaired += 1

            if signal_3_cP != "NA":                                 # valid signal in cP library
                if not CA:
                    h2o_cP_joint_occ_reference_set.append([float(signal_3_cP), p_unpaired])
                else:
                    h2o_cP_joint_occ_reference_set_CA.append([float(signal_3_cP), p_unpaired])

                if float(signal_3_cP) == 0.0:
                    h2o_ref_cP_zero_all += 1
                    if structure == "0":
                        h2o_ref_cP_zero_unpaired += 1

            if signal_3_OH != "NA":                                 # valid signal in OH library
                if not CA:
                    h2o_OH_joint_occ_reference_set.append([float(signal_3_OH), p_unpaired])
                else:
                    h2o_OH_joint_occ_reference_set_CA.append([float(signal_3_OH), p_unpaired])

                if float(signal_3_OH) == 0.0:
                    h2o_ref_OH_zero_all += 1
                    if structure == "0":
                        h2o_ref_OH_zero_unpaired += 1


s=[]
p=[]
cP_joint_occ_reference_set_shuffeled=[]
for elem in cP_joint_occ_reference_set:
    s.append(elem[0])
    p.append(elem[1])
random.shuffle(s)
for i, elem in enumerate(s):
    cP_joint_occ_reference_set_shuffeled.append([elem, p[i]])

"""for the reference set Pb  - NOT CA"""
np_joint_occ_reference_set = np.array(joint_occ_reference_set)  # convert to np.array

"""define the bins manually, index[0]=Signal_cP, [1]=Signal_OH, [2]=unpairedness"""

bins_reference_set = np.array([[0.0, 0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 7.0], [0.0, 0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 7.0], [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]], dtype=object)
hist_reference_set, binedges_reference_set = np.histogramdd(np_joint_occ_reference_set, bins=bins_reference_set, normed=False)
np.set_printoptions(suppress=True)

"""for the reference set Pb  - CA"""
np_joint_occ_reference_set_CA = np.array(joint_occ_reference_set_CA)
bins_reference_set_CA = bins_reference_set
hist_reference_set_CA, binedges_reference_set_CA = np.histogramdd(np_joint_occ_reference_set_CA, bins=bins_reference_set_CA, normed=False)
np.set_printoptions(suppress=True)

"""do this for signal_cP only - NOT CA"""
nP_cP_joint_occ_reference_set = np.array(cP_joint_occ_reference_set)
cP_bins_reference_set = np.array([[0.0, 0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 7.0], [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]], dtype=object)
cP_hist_reference_set, cP_binedges_reference_set = np.histogramdd(nP_cP_joint_occ_reference_set, bins=cP_bins_reference_set, normed=False)

"""do this for signal_cP only - CA"""
nP_cP_joint_occ_reference_set_CA = np.array(cP_joint_occ_reference_set_CA)
cP_bins_reference_set_CA = cP_bins_reference_set
cP_hist_reference_set_CA, cP_binedges_reference_set_CA = np.histogramdd(nP_cP_joint_occ_reference_set_CA, bins=cP_bins_reference_set_CA, normed=False)

"""do this for shuffeled signal_cP only - NOT CA"""
nP_cP_joint_occ_reference_set_shuffeled = np.array(cP_joint_occ_reference_set_shuffeled)  # convert to np.
cP_hist_reference_set_shuffeled, cP_binedges_reference_set_shuffeled = np.histogramdd(nP_cP_joint_occ_reference_set_shuffeled, bins=cP_bins_reference_set, normed=False)

"""do this for signal_OH only - NOT CA"""
nP_OH_joint_occ_reference_set = np.array(OH_joint_occ_reference_set)
OH_bins_reference_set = cP_bins_reference_set
OH_hist_reference_set, OH_binedges_reference_set = np.histogramdd(nP_OH_joint_occ_reference_set, bins=OH_bins_reference_set, normed=False)

"""do this for signal_OH only - CA"""
nP_OH_joint_occ_reference_set_CA = np.array(OH_joint_occ_reference_set_CA)
OH_bins_reference_set_CA = cP_bins_reference_set
OH_hist_reference_set_CA, OH_binedges_reference_set_CA = np.histogramdd(nP_OH_joint_occ_reference_set_CA, bins=OH_bins_reference_set_CA, normed=False)


"""for the reference set H2O - NOT CA"""
h2o_np_joint_occ_reference_set = np.array(h2o_joint_occ_reference_set)
h2o_bins_reference_set = np.array([[0.0, 0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 7.0], [0.0, 0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 7.0], [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]], dtype=object)
h2o_hist_reference_set, h2o_binedges_reference_set = np.histogramdd(h2o_np_joint_occ_reference_set, bins=h2o_bins_reference_set, normed=False)
np.set_printoptions(suppress=True)

"""for the reference set H2O - CA"""
h2o_np_joint_occ_reference_set_CA = np.array(h2o_joint_occ_reference_set_CA)
h2o_bins_reference_set_CA = h2o_bins_reference_set
h2o_hist_reference_set_CA, h2o_binedges_reference_set_CA = np.histogramdd(h2o_np_joint_occ_reference_set_CA, bins=h2o_bins_reference_set_CA, normed=False)
np.set_printoptions(suppress=True)

"""do this for signal_cP only - NOT CA"""
h2o_nP_cP_joint_occ_reference_set = np.array(h2o_cP_joint_occ_reference_set)
h2o_cP_bins_reference_set = np.array([[0.0, 0.25, 0.5, 0.75, 1.0, 2.0, 3.0, 7.0], [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]], dtype=object)
h2o_cP_hist_reference_set, h2o_cP_binedges_reference_set = np.histogramdd(h2o_nP_cP_joint_occ_reference_set, bins=h2o_cP_bins_reference_set, normed=False)

"""do this for signal_cP only - CA"""
h2o_nP_cP_joint_occ_reference_set_CA = np.array(h2o_cP_joint_occ_reference_set_CA)
h2o_cP_bins_reference_set_CA = h2o_cP_bins_reference_set
h2o_cP_hist_reference_set_CA, h2o_cP_binedges_reference_set_CA = np.histogramdd(h2o_nP_cP_joint_occ_reference_set_CA, bins=h2o_cP_bins_reference_set_CA, normed=False)

"""do this for signal_OH only - NOT CA"""
h2o_nP_OH_joint_occ_reference_set = np.array(h2o_OH_joint_occ_reference_set)
h2o_OH_bins_reference_set = h2o_cP_bins_reference_set
h2o_OH_hist_reference_set, h2o_OH_binedges_reference_set = np.histogramdd(h2o_nP_OH_joint_occ_reference_set, bins=h2o_OH_bins_reference_set, normed=False)

"""do this for signal_OH only - CA"""
h2o_nP_OH_joint_occ_reference_set_CA = np.array(h2o_OH_joint_occ_reference_set_CA)
h2o_OH_bins_reference_set_CA = h2o_cP_bins_reference_set
h2o_OH_hist_reference_set_CA, h2o_OH_binedges_reference_set_CA = np.histogramdd(h2o_nP_OH_joint_occ_reference_set_CA, bins=h2o_OH_bins_reference_set_CA, normed=False)


def plot_histos_with_scatter(joint_occ, sample_name):
    x = []
    y = []
    for entry in joint_occ:
        x.append(entry[0])
        y.append(entry[1])

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]

    # start with a rectangular Figure
    fig = plt.figure(figsize=(8, 8))

    ax_scatter = plt.axes(rect_scatter)
    ax_scatter.tick_params(direction='in', top=True, right=True)
    ax_histx = plt.axes(rect_histx)
    ax_histx.tick_params(direction='in', labelbottom=False)
    ax_histy = plt.axes(rect_histy)
    ax_histy.tick_params(direction='in', labelleft=False)

    # the scatter plot:
    ax_scatter.scatter(x, y, marker='.', color='#280592')

    # now determine nice limits by hand:
    binwidth_x = 0.5
    binwidth_y = 0.1
    lim_x = np.ceil(np.abs(x).max() / binwidth_x) * binwidth_x
    lim_y = np.ceil(np.abs(y).max() / binwidth_y) * binwidth_y
    ax_scatter.set_xlim((0, lim_x))
    ax_scatter.set_ylim((0, lim_y))
    ax_scatter.set_xlabel('Probing Signal')
    ax_scatter.set_ylabel('P(unpaired)')
    ax_scatter.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], minor=False)
    ax_scatter.set_axisbelow(True)
    ax_scatter.grid(True)

    # the histograms
    bins_x = np.arange(0, lim_x + binwidth_x, binwidth_x)
    bins_y = np.arange(0, lim_y + binwidth_y, binwidth_y)
    ax_histx.hist(x, bins=bins_x, color='#A333AE')
    ax_histy.hist(y, bins=bins_y, orientation='horizontal', color='#F5954F')
    ax_histx.set_ylabel('count')
    ax_histy.set_xlabel('count')

    ax_histx.set_xlim(ax_scatter.get_xlim())
    ax_histy.set_ylim(ax_scatter.get_ylim())

    #fig.suptitle(output)
    plt.savefig(output_folder+"/"+sample_name + "_scatter_with_histos.pdf")
    plt.close()

def plot3DScatter(data, sample_name):
    """ Setup a 3D figure and plot points """
    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    ax1.plot(data[:, 0], data[:, 1], data[:, 2], 'k.', alpha=0.3)
    ax1.set_xlabel('S(cP)')
    ax1.set_ylabel('S(OH)')
    ax1.set_zlabel('P(unpaired)')
    # plt.savefig(output_folder+"/"+sample_name+"3dscatter.pdf")
    plt.show()
    plt.close()

def calcCondProb_histo(s_cP, s_OH, histo, bins):
    num_datapoints = np.sum(histo)
    cP_bin_length = len(bins[0])
    OH_bin_length = len(bins[1])

    """find corresponding bins in the histogram"""
    cP_bin = np.digitize([s_cP], bins[0][0:cP_bin_length-1]) - 1      #make it 0-based because we use it as list index
    OH_bin = np.digitize([s_OH], bins[1][0:OH_bin_length-1]) - 1     #make it 0-based because we use it as list index

    p_joint_up_cP_OH = float(np.sum(histo[cP_bin[0]][OH_bin[0]][5:10]) / num_datapoints)  # P(unpaired, S(cP), S(OH))
    p_joint_cP_OH = float(np.sum(histo[cP_bin[0]][OH_bin[0]]) / num_datapoints)  # P(S(cP),S(OH))
    if p_joint_cP_OH==0.00:
        unpaired_prob = 0.00                                        # avoid division by zero if bin is empty

    else:
        unpaired_prob = float(p_joint_up_cP_OH / p_joint_cP_OH)

    return unpaired_prob


def calcCondProb_single_sample_histo(signal, histo, bins):
    num_datapoints = np.sum(histo)
    bin_count=len(bins[0])

    """find corresponding bins in the histogram"""

    signal_bin = np.digitize([signal], bins[0][0:bin_count-1]) - 1                  #we need t exclude the upper bound from bins here so that max count do not form seperate bin that does not exist in histo
    p_joint = float(np.sum(histo[signal_bin[0]][5:10]) / num_datapoints)  # P(unpaired, S)
    p_S = float(np.sum(histo[signal_bin[0]]) / num_datapoints)            # P(S)
    unpaired_prob = float(p_joint / p_S)

    return unpaired_prob

def plotHist(histo, bins, sample_name):
    x3 = []  # cP-Signal on x-axis
    y3 = []  # OH-Signal on y-axis
    z3 = np.zeros((len(bins[0]) - 1) * (len(bins[1]) - 1))
    dx = []  # the bin width
    dy = []  # the bin width
    dz = []  # the cond. prob to be unpaired given both signals

    """calculate: P(unpaired | S(cP), S(OH)) = P(unpaired, S(cP), S(OH))/P(S(cP),S(OH))
    loop over the bin combinations of S(cP) S(OH)"""
    for signal_cP_bin in range(0, len(bins[0]) - 1):
        for signal_OH_bin in range(0, len(bins[1]) - 1):
            x3.append(bins[0][signal_cP_bin])
            y3.append(bins[1][signal_OH_bin])
            dx.append(bins[0][signal_cP_bin + 1] - bins[0][signal_cP_bin])
            dy.append(bins[1][signal_OH_bin + 1] - bins[1][signal_OH_bin])
            """take left bin edge as exemplary value for the bin to calc cond. prob"""
            s_cP = bins[0][signal_cP_bin]
            s_OH = bins[1][signal_OH_bin]
            qi=calcCondProb_histo(s_cP, s_OH, histo, bins)
            dz.append(qi)

    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    ax1.bar3d(x3, y3, z3, dx, dy, dz, color='orange', edgecolor='black', alpha=0.8)
    ax1.set_xlabel('S(cP)')
    ax1.set_ylabel('S(OH)')
    ax1.set_zlabel('P(unpaired)')
    plt.savefig(output_folder+"/"+sample_name+"_qi_estimation_hist.pdf")
    # plt.show()
    plt.close()


def fsigmoid(x, k, a, b):
    return 1.0 / (1.0 + np.exp(-k*(x-a)))+b

def fsigmoid2D(data_tuple, k, m, n, a, b):
    (x, y) = data_tuple
    return 1.0 / (1.0 + np.exp(-k*(m*x+n*y-a)))+b

def fitEstimationCurve_single_sample(histo, bins, zero, zero_up):
    x_vals = []
    y_vals = []

    """use bincenters as values for the fit"""
    bin_left_edge = bins[0]
    bincenters = np.array([0.5 * (bin_left_edge[i] + bin_left_edge[i+1]) for i in range(len(bin_left_edge)-1)])
    """loop over the bins"""
    for val in bincenters:
        qi=calcCondProb_single_sample_histo(val, histo, bins)
        x_vals.append(val)
        y_vals.append(qi)

    xdata = np.array(x_vals)
    ydata = np.array(y_vals)

    sigma = np.ones(len(xdata))

    popt, pcov = curve_fit(fsigmoid, xdata, ydata, sigma=sigma)
    return popt

ref_cP_fit_params = fitEstimationCurve_single_sample(cP_hist_reference_set, cP_bins_reference_set, ref_cP_zero_all, ref_cP_zero_unpaired)
ref_cP_fit_params_CA = fitEstimationCurve_single_sample(cP_hist_reference_set_CA, cP_bins_reference_set_CA, ref_cP_zero_all, ref_cP_zero_unpaired)
ref_OH_fit_params = fitEstimationCurve_single_sample(OH_hist_reference_set, OH_bins_reference_set, ref_OH_zero_all, ref_OH_zero_unpaired)
ref_OH_fit_params_CA = fitEstimationCurve_single_sample(OH_hist_reference_set_CA, OH_bins_reference_set_CA, ref_OH_zero_all, ref_OH_zero_unpaired)
h2o_ref_cP_fit_params = fitEstimationCurve_single_sample(h2o_cP_hist_reference_set, h2o_cP_bins_reference_set, h2o_ref_cP_zero_all, h2o_ref_cP_zero_unpaired)
h2o_ref_cP_fit_params_CA = fitEstimationCurve_single_sample(h2o_cP_hist_reference_set_CA, h2o_cP_bins_reference_set_CA, h2o_ref_cP_zero_all, h2o_ref_cP_zero_unpaired)
h2o_ref_OH_fit_params = fitEstimationCurve_single_sample(h2o_OH_hist_reference_set, h2o_OH_bins_reference_set, h2o_ref_OH_zero_all, h2o_ref_OH_zero_unpaired)
h2o_ref_OH_fit_params_CA = fitEstimationCurve_single_sample(h2o_OH_hist_reference_set_CA, h2o_OH_bins_reference_set_CA, h2o_ref_OH_zero_all, h2o_ref_OH_zero_unpaired)
ref_cP_fit_params_shuffeled = fitEstimationCurve_single_sample(cP_hist_reference_set_shuffeled, cP_bins_reference_set, ref_cP_zero_all, ref_cP_zero_unpaired)


def fitEstimationCurve(histo, bins, zero_zero, zero_zero_up):
    x_vals = []
    y_vals = []
    z_vals = []

    """use bincenters as values for the fit"""
    cP_bin_left_edge = bins[0]
    OH_bin_left_edge = bins[1]
    cP_bincenters = np.array([0.5 * (cP_bin_left_edge[i] + cP_bin_left_edge[i+1]) for i in range(len(cP_bin_left_edge)-1)])
    OH_bincenters = np.array([0.5 * (OH_bin_left_edge[i] + OH_bin_left_edge[i+1]) for i in range(len(OH_bin_left_edge)-1)])

    """loop over the bin combinations of S(cP) S(OH)"""
    for s_cP in cP_bincenters:
        for s_OH in OH_bincenters:
            qi=calcCondProb_histo(s_cP, s_OH, histo, bins)
            if qi != 0.00:                                              # exclude bins without values
                x_vals.append(s_cP)
                y_vals.append(s_OH)
                z_vals.append(qi)

    xdata = np.array(x_vals)
    ydata = np.array(y_vals)
    zdata = np.array(z_vals)

    sigma = np.ones(len(xdata))

    popt, pcov = curve_fit(fsigmoid2D, (xdata, ydata), zdata, bounds=([-np.inf, 0.5, 0.5, -np.inf, -1.0 ], [+np.inf, 1.5, 1.5, +np.inf, 1.0]), sigma=sigma)
    return popt

fit_params_2D_reference_set = fitEstimationCurve(hist_reference_set, bins_reference_set, ref_zero_zero_all, ref_zero_zero_unpaired)
fit_params_2D_reference_set_CA = fitEstimationCurve(hist_reference_set_CA, bins_reference_set_CA, ref_zero_zero_all, ref_zero_zero_unpaired)
h2o_fit_params_2D_reference_set = fitEstimationCurve(h2o_hist_reference_set, h2o_bins_reference_set, h2o_ref_zero_zero_all, h2o_ref_zero_zero_unpaired)
h2o_fit_params_2D_reference_set_CA = fitEstimationCurve(h2o_hist_reference_set_CA, h2o_bins_reference_set_CA, h2o_ref_zero_zero_all, h2o_ref_zero_zero_unpaired)


def plotEstimationCurve_single_sample(histo, bins, sample_name, fit_params, zero, zero_up, color, colorfit, label):
    x_vals = []
    y_vals = []


    bin_left_edge = bins[0]
    bincenters = np.array([0.5 * (bin_left_edge[i] + bin_left_edge[i+1]) for i in range(len(bin_left_edge)-1)])


    """loop over the bins"""
    for val in bincenters:
        qi=calcCondProb_single_sample_histo(val, histo, bins)
        x_vals.append(val)
        y_vals.append(qi)
    binwidth=np.array([(bin_left_edge[i+1] - bin_left_edge[i]) for i in range(len(bin_left_edge)-1)])
    plt.bar(bincenters, fsigmoid(bincenters, *fit_params), width=binwidth, edgecolor='lightgray', color='None')

    xdata = np.array(x_vals)
    ydata = np.array(y_vals)

    """plot fitting curve"""
    x = np.linspace(0, 7, 100)
    y = fsigmoid(x, *fit_params)

    eq=r'$\frac{1}{1+ e^{-k\times(x-a)}}+b$'
    plt.plot(x, y,color=colorfit)

    plt.plot(xdata,ydata,'o', color=color, label=label)
    plt.ylim(bottom=0, top=1.05)
    plt.xlabel(r'$S$')
    plt.ylabel(r'$p(S)$')
    plt.title(sample_name)


def plot2DEstimationCurve(histo, bins, sample_name, fit_params, zero_zero, zero_zero_up):
    x_vals = []
    y_vals = []
    z_vals = []

    """use bincenters as values for the fit"""
    cP_bin_left_edge = bins[0]
    OH_bin_left_edge = bins[1]
    cP_bincenters = np.array([0.5 * (cP_bin_left_edge[i] + cP_bin_left_edge[i+1]) for i in range(len(cP_bin_left_edge)-1)])
    OH_bincenters = np.array([0.5 * (OH_bin_left_edge[i] + OH_bin_left_edge[i+1]) for i in range(len(OH_bin_left_edge)-1)])

    """loop over the bin combinations of S(cP) S(OH)"""
    for s_cP in cP_bincenters:
        for s_OH in OH_bincenters:
            qi=calcCondProb_histo(s_cP, s_OH, histo, bins)
            if qi != 0.00:
                x_vals.append(s_cP)
                y_vals.append(s_OH)
                z_vals.append(qi)

    xdata = np.array(x_vals)
    ydata = np.array(y_vals)
    zdata = np.array(z_vals)

    """for the fit"""
    x = np.linspace(0, 7, 25)
    y = np.linspace(0, 7, 25)
    X, Y = np.meshgrid(x, y)
    Z = fsigmoid2D((X, Y), *fit_params)

    fig = plt.figure(figsize = (8, 6))
    ax = fig.add_subplot(projection='3d')
    ax.scatter(xdata, ydata, zdata,marker='o', color="#ca0020")
    ax.plot_wireframe(X, Y, Z, color="#0571b0")
    ax.set_xlabel(r'$S^{cP}$', labelpad=12)
    ax.set_ylabel(r'$S^{OH}$',labelpad=14)
    ax.set_title(sample_name)
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel(r'$p(S^{cP},S^{OH})$', rotation=90, labelpad=15)
    ax.set_zlim(bottom=0.0, top=1.0)
    ax.view_init(19, -146)
    plt.savefig(output_folder+"/"+sample_name+"_qi_estimation_curve2D.svg")
    plt.close()


#####################################################################################################
"""do stuff"""

plotHist(hist_reference_set, bins_reference_set, "both_reference_NOT_CA")
plotHist(hist_reference_set_CA, bins_reference_set_CA, "both_reference_CA")

fig = plt.figure(figsize=(5,5))
plotEstimationCurve_single_sample(cP_hist_reference_set, cP_bins_reference_set, "", ref_cP_fit_params, ref_cP_zero_all, ref_cP_zero_unpaired, "#0571b0", "#0571b0", "not CA")
plotEstimationCurve_single_sample(cP_hist_reference_set_CA, cP_bins_reference_set_CA, "2'3'-cP", ref_cP_fit_params_CA, ref_cP_zero_all, ref_cP_zero_unpaired, "#f4a582","#f08556", "CA")
plt.legend(loc="upper left", bbox_to_anchor=(0.417,0.3))
plt.tight_layout()
plt.savefig(output_folder + "/" + "cP_reference" + "_estimation_curve.pdf")
plt.close()

fig = plt.figure(figsize=(5,5))
plotEstimationCurve_single_sample(OH_hist_reference_set, OH_bins_reference_set, "", ref_OH_fit_params, ref_OH_zero_all, ref_OH_zero_unpaired, "#0571b0", "#0571b0", "not CA")
plotEstimationCurve_single_sample(OH_hist_reference_set_CA, OH_bins_reference_set_CA, "5'-OH", ref_OH_fit_params_CA, ref_OH_zero_all, ref_OH_zero_unpaired, "#f4a582","#f08556", "CA")
plt.legend(loc="upper left", bbox_to_anchor=(0.417,0.3))
plt.tight_layout()
plt.savefig(output_folder + "/" + "OH_reference" + "_estimation_curve.pdf")
plt.close()

fig = plt.figure(figsize=(5,5))
plotEstimationCurve_single_sample(cP_hist_reference_set_shuffeled, cP_bins_reference_set, "", ref_cP_fit_params_shuffeled, ref_cP_zero_all, ref_cP_zero_unpaired,"#0571b0", "#0571b0", "shuffeled data")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig(output_folder + "/" + "shuffeled" + "_estimation_curve.pdf")
plt.close()

plotHist(h2o_hist_reference_set, h2o_bins_reference_set, "h2o_both_reference_NOT_CA")
plotHist(h2o_hist_reference_set_CA, h2o_bins_reference_set_CA, "h2o_both_reference_CA")
fig = plt.figure(figsize = (5, 5))
plotEstimationCurve_single_sample(h2o_cP_hist_reference_set, h2o_cP_bins_reference_set, "", h2o_ref_cP_fit_params, h2o_ref_cP_zero_all, h2o_ref_cP_zero_unpaired, "#0571b0", "#0571b0", "2'3'-cP $\mathregular{Pb^{2\!+}}\!(-)$")
plotEstimationCurve_single_sample(h2o_OH_hist_reference_set, h2o_OH_bins_reference_set, "not CA", h2o_ref_OH_fit_params, h2o_ref_OH_zero_all, h2o_ref_OH_zero_unpaired, "#f4a582","#f08556", "5'-OH $\mathregular{Pb^{2\!+}}\!(-)$")
plt.legend(loc="upper left", bbox_to_anchor=(0.15,0.4))
plt.tight_layout()
plt.savefig(output_folder + "/" + "h2o_reference" + "_NOT_CA_estimation_curve.pdf")
plt.close()

fig = plt.figure(figsize = (5, 5))
plotEstimationCurve_single_sample(h2o_cP_hist_reference_set_CA, h2o_cP_bins_reference_set_CA, "", h2o_ref_cP_fit_params_CA, h2o_ref_cP_zero_all, h2o_ref_cP_zero_unpaired, "#0571b0", "#0571b0", "2'3'-cP $\mathregular{Pb^{2\!+}}\!(-)$")
plotEstimationCurve_single_sample(h2o_OH_hist_reference_set_CA, h2o_OH_bins_reference_set_CA, "CA", h2o_ref_OH_fit_params_CA, h2o_ref_OH_zero_all, h2o_ref_OH_zero_unpaired, "#f4a582","#f08556", "5'-OH $\mathregular{Pb^{2\!+}}\!(-)$")
plt.legend(loc="upper left", bbox_to_anchor=(0.15,0.4))
plt.tight_layout()
plt.savefig(output_folder + "/" + "h2o_reference" + "_CA_estimation_curve.pdf")
plt.close()

plt.rcParams.update({'font.size':21})
plt.rcParams.update({'axes.titlesize':21})
plt.rcParams['xtick.major.pad']= 1
plt.rc('legend', fontsize=21)
plt.rc('xtick', labelsize=21)
plt.rc('ytick', labelsize=21)

plot2DEstimationCurve(hist_reference_set, bins_reference_set, "not CA", fit_params_2D_reference_set, ref_zero_zero_all, ref_zero_zero_unpaired)
plot2DEstimationCurve(hist_reference_set_CA, bins_reference_set_CA, "CA", fit_params_2D_reference_set_CA, ref_zero_zero_all, ref_zero_zero_unpaired)

plot2DEstimationCurve(h2o_hist_reference_set, h2o_bins_reference_set, "h2o_both_reference_NOT_CA", h2o_fit_params_2D_reference_set, h2o_ref_zero_zero_all, h2o_ref_zero_zero_unpaired)
plot2DEstimationCurve(h2o_hist_reference_set_CA, h2o_bins_reference_set_CA, "h2o_both_reference_CA", h2o_fit_params_2D_reference_set_CA, h2o_ref_zero_zero_all, h2o_ref_zero_zero_unpaired)



with open(output_folder+"/"+output_file_name, 'w') as output_file:
    for transcript_id in benchmark_transcripts:
        for pos in signal[transcript_id]:
            CA = False
            if signal[transcript_id][pos]["CA"]:
                CA = True

            signal_1_OH = signal[transcript_id][pos]["1-OH"]
            signal_2_OH = signal[transcript_id][pos]["2-OH"]
            signal_3_OH = signal[transcript_id][pos]["3-OH"]
            signal_1_cP = signal[transcript_id][pos]["1-cP"]
            signal_2_cP = signal[transcript_id][pos]["2-cP"]
            signal_3_cP = signal[transcript_id][pos]["3-cP"]
            qi_both = "-"
            qi_cP = "-"
            qi_OH = "-"
            h2o_qi_both = "-"
            h2o_qi_cP = "-"
            h2o_qi_OH = "-"

            if all(elem != "NA" for elem in [signal_1_OH,signal_2_OH,signal_1_cP,signal_2_cP]):         # valid signal in both libraries
                if not CA:
                    qi_both = fsigmoid2D((statistics.mean([float(signal_1_cP), float(signal_2_cP)]), statistics.mean([float(signal_1_OH), float(signal_2_OH)])), *fit_params_2D_reference_set)
                else:
                    qi_both = fsigmoid2D((statistics.mean([float(signal_1_cP), float(signal_2_cP)]),
                                          statistics.mean([float(signal_1_OH), float(signal_2_OH)])),
                                         *fit_params_2D_reference_set_CA)

            if all(elem != "NA" for elem in [signal_1_cP,signal_2_cP]):                                 # valid signal in cP libraries
                if not CA:
                    qi_cP = fsigmoid(statistics.mean([float(signal_1_cP), float(signal_2_cP)]), *ref_cP_fit_params)
                else:
                    qi_cP = fsigmoid(statistics.mean([float(signal_1_cP), float(signal_2_cP)]), *ref_cP_fit_params_CA)

            if all(elem != "NA" for elem in [signal_1_OH,signal_2_OH]):                                 # valid signal in OH libraries
                if not CA:
                    qi_OH=fsigmoid(statistics.mean([float(signal_1_OH), float(signal_2_OH)]), *ref_OH_fit_params)
                else:
                    qi_OH = fsigmoid(statistics.mean([float(signal_1_OH), float(signal_2_OH)]), *ref_OH_fit_params_CA)

            if all(elem != "NA" for elem in [signal_3_OH, signal_3_cP]):                                # valid signal in both H2O libraries
                if not CA:
                    h2o_qi_both=fsigmoid2D((float(signal_3_cP), float(signal_3_OH)), *h2o_fit_params_2D_reference_set)
                else:
                    h2o_qi_both = fsigmoid2D((float(signal_3_cP), float(signal_3_OH)), *h2o_fit_params_2D_reference_set_CA)

            if signal_3_cP != "NA":                                                                     # valid signal in cP H2O library
                if not CA:
                    h2o_qi_cP=fsigmoid(float(signal_3_cP), *h2o_ref_cP_fit_params)
                else:
                    h2o_qi_cP = fsigmoid(float(signal_3_cP), *h2o_ref_cP_fit_params_CA)

            if signal_3_OH != "NA":                                                                     # valid signal in OH H2O library
                if not CA:
                    h2o_qi_OH=fsigmoid(float(signal_3_OH), *h2o_ref_OH_fit_params)
                else:
                    h2o_qi_OH = fsigmoid(float(signal_3_OH), *h2o_ref_OH_fit_params_CA)


            output_file.write(transcript_id + "\t" + str(pos) + "\t" + signal_1_cP + "\t" + signal_2_cP + "\t" + signal_3_cP + "\t" + signal_1_OH + "\t" + signal_2_OH + "\t" + signal_3_OH + "\t" + str(qi_cP) + "\t" + str(qi_OH) + "\t" + str(qi_both) + "\t" + str(h2o_qi_cP) + "\t" + str(h2o_qi_OH) + "\t" + str(h2o_qi_both) + "\n")

    """add rRNA domains as special case"""
    for domain in [("16S_ribosomal_RNA_domain_I", 1, 559), ("16S_ribosomal_RNA_domain_II", 560, 913),
                   ("16S_ribosomal_RNA_domain_III", 914, 1397), ("16S_ribosomal_RNA_domain_IV", 1398, 1542),
                   ("23S_ribosomal_RNA_domain_I", 22, 534), ("23S_ribosomal_RNA_domain_II", 587, 1268),
                   ("23S_ribosomal_RNA_domain_III", 1283, 1654), ("23S_ribosomal_RNA_domain_IV", 1655, 2016),
                   ("23S_ribosomal_RNA_domain_V", 2030, 2634), ("23S_ribosomal_RNA_domain_VI", 2637, 2898)]:
        transcript_id_full, start, end = domain[0], domain[1], domain[2]
        transcript_id = transcript_id_full.split("_")[0] + "_" + transcript_id_full.split("_")[1] + "_" + \
                        transcript_id_full.split("_")[2]
        i = 0
        for pos in range(start, end + 1):
            i += 1
            CA = False
            if signal[transcript_id][pos]["CA"]:
                CA = True
            signal_1_OH = signal[transcript_id][pos]["1-OH"]
            signal_2_OH = signal[transcript_id][pos]["2-OH"]
            signal_3_OH = signal[transcript_id][pos]["3-OH"]
            signal_1_cP = signal[transcript_id][pos]["1-cP"]
            signal_2_cP = signal[transcript_id][pos]["2-cP"]
            signal_3_cP = signal[transcript_id][pos]["3-cP"]
            qi_both = "-"
            qi_cP = "-"
            qi_OH = "-"
            h2o_qi_both = "-"
            h2o_qi_cP = "-"
            h2o_qi_OH = "-"

            if all(elem != "NA" for elem in
                   [signal_1_OH, signal_2_OH, signal_1_cP, signal_2_cP]):  # valid signal in both libraries
                if not CA:
                    qi_both = fsigmoid2D((statistics.mean([float(signal_1_cP), float(signal_2_cP)]),
                                          statistics.mean([float(signal_1_OH), float(signal_2_OH)])),
                                         *fit_params_2D_reference_set)
                else:
                    qi_both = fsigmoid2D((statistics.mean([float(signal_1_cP), float(signal_2_cP)]),
                                          statistics.mean([float(signal_1_OH), float(signal_2_OH)])),
                                         *fit_params_2D_reference_set_CA)

            if all(elem != "NA" for elem in [signal_1_cP, signal_2_cP]):  # valid signal in cP libraries
                if not CA:
                    qi_cP = fsigmoid(statistics.mean([float(signal_1_cP), float(signal_2_cP)]), *ref_cP_fit_params)
                else:
                    qi_cP = fsigmoid(statistics.mean([float(signal_1_cP), float(signal_2_cP)]), *ref_cP_fit_params_CA)

            if all(elem != "NA" for elem in [signal_1_OH, signal_2_OH]):  # valid signal in OH libraries
                if not CA:
                    qi_OH = fsigmoid(statistics.mean([float(signal_1_OH), float(signal_2_OH)]), *ref_OH_fit_params)
                else:
                    qi_OH = fsigmoid(statistics.mean([float(signal_1_OH), float(signal_2_OH)]), *ref_OH_fit_params_CA)

            if all(elem != "NA" for elem in [signal_3_OH, signal_3_cP]):  # valid signal in both H2O libraries
                if not CA:
                    h2o_qi_both = fsigmoid2D((float(signal_3_cP), float(signal_3_OH)), *h2o_fit_params_2D_reference_set)
                else:
                    h2o_qi_both = fsigmoid2D((float(signal_3_cP), float(signal_3_OH)),
                                             *h2o_fit_params_2D_reference_set_CA)

            if signal_3_cP != "NA":  # valid signal in cP H2O library
                if not CA:
                    h2o_qi_cP = fsigmoid(float(signal_3_cP), *h2o_ref_cP_fit_params)
                else:
                    h2o_qi_cP = fsigmoid(float(signal_3_cP), *h2o_ref_cP_fit_params_CA)

            if signal_3_OH != "NA":  # valid signal in OH H2O library
                if not CA:
                    h2o_qi_OH = fsigmoid(float(signal_3_OH), *h2o_ref_OH_fit_params)
                else:
                    h2o_qi_OH = fsigmoid(float(signal_3_OH), *h2o_ref_OH_fit_params_CA)

            output_file.write(transcript_id_full + "\t" + str(
                i) + "\t" + signal_1_cP + "\t" + signal_2_cP + "\t" + signal_3_cP + "\t" + signal_1_OH + "\t" + signal_2_OH + "\t" + signal_3_OH + "\t" + str(qi_cP) + "\t" + str(qi_OH) + "\t" + str(qi_both) + "\t" + str(h2o_qi_cP) + "\t" + str(h2o_qi_OH) + "\t" + str(h2o_qi_both)+ "\n")
import sys
import matplotlib.pyplot as plt
import numpy as np
import math
plt.rcParams.update({'font.size': 19})
plt.rcParams.update({'axes.titlesize': 19})
plt.rc('xtick', labelsize=19)
plt.rc('ytick', labelsize=19)

signal_file = sys.argv[1]
output_folder = sys.argv[2]

signal = {}
benchmark_transcripts = set()

counter= 0

with open(signal_file) as file:
    for line in file:
        line=line.rstrip()
        columns = line.split("\t")
        transcript_id, nc, position, structure, sequence, library, library_type = columns[5],columns[7],int(columns[8]),columns[9],columns[3],columns[10], columns[11]
        if library not in signal:
            signal[library] = {}
        if "type" not in signal[library]:
            signal[library]["type"] = library_type
        if transcript_id not in signal[library]:
            signal[library][transcript_id] = [[],[],[],[]]      #0: position; 1: normalized counts; 2: sequence; 3: structure
        if structure != "-":                                    #reference structure is known for this transcript
            benchmark_transcripts.add(transcript_id)
            signal[library][transcript_id][0].append(position)
            signal[library][transcript_id][1].append(nc)
            signal[library][transcript_id][2].append(sequence)
            signal[library][transcript_id][3].append(int(structure))
            counter+=1

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order + 1)
    half_window = (window_size - 1) // 2
    # precompute coefficients
    b = np.mat([[k ** i for i in order_range] for k in range(-half_window, half_window + 1)])
    m = np.linalg.pinv(b).A[deriv] * rate ** deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs(y[1:half_window + 1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window - 1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode='valid')



# counter = 0
# for transcript_id in benchmark_transcripts:
#     transcript_length = len(signal["1-cP"][transcript_id][0])
#     segment_counts=math.ceil(transcript_length/100)
#     segment_ends=[]
#     if transcript_length < 100:
#         segment_ends.append(transcript_length)
#     else:
#         for i in range(1,segment_counts+1):
#             segment_ends.append(i*100)
#     for segment in segment_ends:
#         fig = plt.figure(figsize=(10,4))
#         ax = fig.gca()
#         fig2 = plt.figure(figsize=(10,4))
#         ax2 = fig2.gca()
#         cmap = plt.get_cmap('turbo')
#         counter2 = 0
#         for library in signal:
#             positions_up = []
#             vals_up = []
#             positions_p = []
#             vals_p = []
#             vals=[]
#             positions = []
#             for i, nc in enumerate(signal[library][transcript_id][1]):
#                 if signal[library][transcript_id][0][i] >= segment - 100 + 1 and signal[library][transcript_id][0][
#                     i] < segment + 1:
#                     structure = signal[library][transcript_id][3][i]
#                     if structure == 0:
#                         positions_up.append(signal[library][transcript_id][0][i])
#                         vals_up.append(-0.5)
#                     elif structure == 1:
#                         positions_p.append(signal[library][transcript_id][0][i])
#                         vals_p.append(-0.5)
#
#                     if nc != "NA":
#                         nc = float(nc)
#                         positions.append(signal[library][transcript_id][0][i])
#                         vals.append(nc)
#
#             color = ""
#
#             if library.split("-")[0]=="1":
#                 color = 'blue'
#             elif library.split("-")[0]=="2":
#                 color = 'green'
#             elif library.split("-")[0]=="3":
#                 color = 'red'
#
#             if signal[library]["type"]=="cP":
#                 ax.plot(positions, vals, '-', marker = ".", color=color, alpha=0.5, label=library)
#             else:
#                 ax2.plot(positions, vals, '-', marker = ".", color=color, alpha=0.5, label=library)
#
#         for plot in [ax, ax2]:
#             plot.plot(positions_up, vals_up, 'o', color="lightpink", label="unpaired")
#             plot.plot(positions_p, vals_p, '.', color="lightpink", label="paired")
#             if transcript_length < 100:
#                 plot.set_xticks(np.arange(1, transcript_length + 1, 5))
#                 plot.set_xticklabels(np.arange(1, transcript_length + 1, 5), rotation='vertical')
#                 plot.set_xlim(left=0, right=transcript_length + 1)
#                 plot.set_ylim(top=7.2, bottom=-0.8)
#             else:
#                 plot.set_xticks(np.arange(segment-100+1, segment+2, 5))
#                 plot.set_xticklabels(np.arange(segment-100+1, segment+2, 5), rotation='vertical')
#                 plot.set_xlim(left=segment-100, right=segment+1)
#                 plot.set_ylim(top=7.2, bottom=-0.8)
#
#             plot.set_xlabel("Nucleotide Position")
#             plot.set_ylabel("Normalized Probing Signal")
#             plot.legend()
#             plot.set_title(transcript_id)
#             plot.spines["top"].set_visible(False)
#             plot.spines["right"].set_visible(False)
#
#         plt.tight_layout()
#         fig.savefig(output_folder + "/"+ transcript_id + "_" + str(segment) + "_cP.pdf")
#         fig2.savefig(output_folder + "/"+ transcript_id + "_" + str(segment) + "_OH.pdf")
#         plt.close('all')

fig3, (ax3, ax4) = plt.subplots(2, 1, sharex=True,  sharey=True, figsize=(10,6))

vals_1_cP=[]
vals_3_cP=[]
positions_cP=[]

vals_1_OH=[]
vals_2_OH=[]
positions_OH=[]

for library in ["1-cP", "2-cP"]:
    transcript_id="tRNA_Leu_CAG_1"
    positions_up = []
    vals_up = []
    positions_p = []
    vals_p = []

    for i, nc in enumerate(signal[library][transcript_id][1]):
        if nc != "NA":
            structure = signal[library][transcript_id][3][i]
            if structure == 0:
                positions_up.append(signal[library][transcript_id][0][i])
                vals_up.append(float(nc))
            elif structure == 1:
                positions_p.append(signal[library][transcript_id][0][i])
                vals_p.append(float(nc))

            if library == "1-cP":
                vals_1_cP.append(float(nc))
                positions_cP.append(signal[library][transcript_id][0][i])
            elif library == "2-cP":
                vals_3_cP.append(float(nc))

    marker = "o"
    color_p="#0571b0"
    color_up="#f4a582"

    ax3.plot(positions_up, vals_up, marker=marker, linestyle="none", mfc=color_up, mec=color_up, alpha=0.6, label="unpaired")
    ax3.plot(positions_p, vals_p, marker=marker, linestyle="none", mfc=color_p, mec=color_p, alpha=0.5, label="paired")

means_cP=np.mean([vals_1_cP, vals_3_cP], axis=0)
ycP= savitzky_golay(means_cP, 5, 3)
ax3.plot(positions_cP,ycP, color='black', linewidth=0.75)

for library in ["1-OH", "2-OH"]:
    transcript_id="tRNA_Leu_CAG_1"
    positions_up = []
    vals_up = []
    positions_p = []
    vals_p = []

    for i, nc in enumerate(signal[library][transcript_id][1]):
        if nc != "NA":
            structure = signal[library][transcript_id][3][i]
            if structure == 0:
                positions_up.append(signal[library][transcript_id][0][i])
                vals_up.append(float(nc))
            elif structure == 1:
                positions_p.append(signal[library][transcript_id][0][i])
                vals_p.append(float(nc))

            if library == "1-OH":
                vals_1_OH.append(float(nc))
                positions_OH.append(signal[library][transcript_id][0][i])
            elif library == "2-OH":
                vals_2_OH.append(float(nc))

    marker = "o"
    color_p = "#0571b0"
    color_up = "#f4a582"

    ax4.plot(positions_up, vals_up, marker=marker, linestyle="none", mfc=color_up, mec=color_up, alpha=0.6,
             label="unpaired")
    ax4.plot(positions_p, vals_p, marker=marker, linestyle="none", mfc=color_p, mec=color_p, alpha=0.5, label="paired")
    if library == "1-OH":
        ax4.legend(loc="upper right")

means_OH=np.mean([vals_1_OH, vals_2_OH], axis=0)
yOH = savitzky_golay(means_OH, 5, 3)
ax4.plot(positions_OH,yOH, color='black', linewidth=0.75)

ax3.text(0, 3.0, "2'3'-cP", fontsize=19)
ax4.text(0, 3.0, "5'-OH", fontsize=19)
ax4.set_xlabel("nucleotide position")
ax3.set_ylabel(r'$S$')
ax4.set_ylabel(r'$S$')
xticks=[1,10,20,30,40,50,60,70,80]
ax3.set_xticks(xticks)
ax4.set_xticks(xticks)

fig3.tight_layout()
fig3.savefig(output_folder + "/" + "tRNA_Leu_CAG_1_profile_Pb.pdf")
plt.close()


"""H2O samples"""


fig4, (ax5, ax6) = plt.subplots(2, 1, sharex=True,  sharey=True, figsize=(10,6))
for library in ["2-cP", "3-OH"]:
    transcript_id="tRNA_Leu_CAG_1"
    positions= []
    vals= []

    for i, nc in enumerate(signal[library][transcript_id][1]):
        if nc != "NA":
            positions.append(signal[library][transcript_id][0][i])
            vals.append(float(nc))

    marker = "o"
    color_p = "#0571b0"
    color_up = "#ca0020"

    if library == "2-cP":
        ax5.plot(positions_cP, means_cP, linestyle="-", linewidth=2, color=color_up, alpha=1,
                 label="$\mathregular{Pb^{2\!+}}\!(+)$")
        ax5.plot(positions, vals,  linestyle="-", linewidth=2, color=color_p, alpha=0.75,
                     label="$\mathregular{Pb^{2\!+}}\!(-)$")

    elif library == "3-OH":
        ax6.plot(positions_OH, means_OH, linestyle="-", linewidth=2,  color=color_up, alpha=1,
                 label="$\mathregular{Pb^{2\!+}}\!(+)$")
        ax6.plot(positions, vals,  linestyle="-", linewidth=2, color=color_p, alpha=0.75,
                 label="$\mathregular{Pb^{2\!+}}\!(-)$")
        ax6.legend(loc="upper right")

ax5.text(0, 4.9, "2'3'-cP", fontsize=19)
ax6.text(0, 4.9, "5'-OH", fontsize=19)
ax6.set_xlabel("nucleotide position")
ax5.set_ylabel(r'$S$')
ax6.set_ylabel(r'$S$')
xticks=[1,10,20,30,40,50,60,70,80]
ax5.set_xticks(xticks)
ax6.set_xticks(xticks)
fig4.tight_layout()
fig4.savefig(output_folder + "/" + "tRNA_Leu_CAG_1_profile_H2O.pdf")
plt.close()



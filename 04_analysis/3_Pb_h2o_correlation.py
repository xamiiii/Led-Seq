import numpy as np
import matplotlib.pyplot as plt
import sys
plt.rcParams.update({'font.size': 19})
plt.rcParams.update({'axes.titlesize': 19})
plt.rc('xtick', labelsize=19)
plt.rc('ytick', labelsize=19)

signal_file = sys.argv[1]


x_cP=[]
y_cP=[]
x_OH=[]
y_OH=[]

with open(signal_file) as file:
    for line in file:
        line=line.rstrip()
        columns = line.split(",")
        transcript_id, position, h2o, pb, library = columns[0],columns[1],columns[3],columns[4],columns[5]
        if library=="P":
            if h2o != "NA" and pb != "NA":
                x_cP.append(float(h2o))
                y_cP.append(float(pb))
        elif library=="OH":
            if h2o != "NA" and pb != "NA":
                x_OH.append(float(h2o))
                y_OH.append(float(pb))

x_cP_log=[]
y_cP_log=[]
x_OH_log=[]
y_OH_log=[]
"log transformation"
for x in x_cP:
    x_cP_log.append(np.log(x+1))
for y in y_cP:
    y_cP_log.append(np.log(y+1))
for x in x_OH:
    x_OH_log.append(np.log(x+1))
for y in y_OH:
    y_OH_log.append(np.log(y+1))



m, b = np.polyfit(x_cP_log, y_cP_log, 1)
n, k = np.polyfit(x_OH_log, y_OH_log, 1)

y_regression_cP=[]
y_regression_OH=[]
for x in x_cP_log:
    y_regression_cP.append((m*x)+b)
for x in x_OH_log:
    y_regression_OH.append((n*x)+k)

fig = plt.figure(figsize=(5,4.7))
ax = fig.gca()
ax.plot(x_cP_log, y_cP_log, color="#0571b0", linestyle="none", markersize=2, marker=".")
ax.plot(x_cP_log, y_regression_cP, linestyle="-", color="#ca0020")
ax.set_ylabel(r'$S_{\mathregular{Pb^{2\!+}}\!(+)}$ (log transformed)')
ax.set_xlabel(r'$S_{\mathregular{Pb^{2\!+}}\!(-)}$ (log transformed)')
plt.tight_layout()
plt.savefig("./cP_pb_h2o_scatter.pdf")
plt.close()


fig = plt.figure(figsize=(5,4.7))
ax = fig.gca()
ax.plot(x_OH_log, y_OH_log, color="#0571b0",linestyle="none", markersize=2, marker=".")
ax.plot(x_OH_log, y_regression_OH, linestyle="-",color="#ca0020")
ax.set_ylabel(r'$S_{\mathregular{Pb^{2\!+}}\!(+)}$ (log transformed)')
ax.set_xlabel(r'$S_{\mathregular{Pb^{2\!+}}\!(-)}$ (log transformed)')
plt.tight_layout()
plt.savefig("./OH_pb_h2o_scatter.pdf")
plt.close()
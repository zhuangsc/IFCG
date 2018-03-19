#MIT License
#
#Copyright (c) 2018 Sicong Zhuang
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.


#!/usr/bin/env python3

import numpy as np
import numpy.linalg as nlin
import os, sys, math, re, statistics
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
import statistics as st


MAT = (
"af_shell8",
"cfd2",
"ecology2",
"consph",
"G2_circuit",
"G3_circuit",
"nd24k",
"thermal2",
" ",
"Mean",
#" ",
#" ",
#" ",
)

VER = ("cg_alg1", "cg_alg4", "cg_alg7", "cg_alg3", "cg_alg4_ifcg", "cg_alg4_ifcg_v2")#, "cg_alg4_v3")
VER = ("cg_alg1", "cg_alg4_ifcg", "cg_alg4_ifcg_centinel", "cg_alg4_ifcg_v2", "cg_alg4_ifcg_v2_centinel")
#COLOR = ('r', 'g', 'b', 'k', 'y', 'm', 'c')
COLOR = ('0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1') #Grayscale
TASKS = 32
FUSE = 20
#PROCESS = [1, 3, 6, 12]
PROCESS = [1, 2, 4, 8, 16]
#PROCESS = [1, 8, 16, 24]
#PROCESS = [ p+1 for p in range(16) ]

VER_NAME = ("PCG", "Pipelined", "Gropp", "Chronopoulos", "IFCG", "IFCG2")
VER_NAME = ("PCG", "ifcg_waiton", "ifcg_sentinel", "ifcg_v2_waiton", "ifcg_v2_sentinel")
PREC0 = "1E-8"
PREC = "1E-7"
CONV_CRIT = float(PREC) * 2


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: {} [dir]".format(sys.argv[0]))
        sys.exit(1)
    ROOT = os.path.abspath(sys.argv[1])

    VELAPSE = list()

    for mat in MAT:
        ELAPSE = list()
        belp = 0
        for ver in VER:
            PELAPSE = list()
            for process in PROCESS:
                if mat == " " or mat == "Mean":
                    PELAPSE.append(np.inf)
                    continue
                DIR = "{}_{}_{!r}_{}_{}".format(ver, TASKS, process, PREC0, FUSE)
                conv = 0
                mat_elp = list()
                logf = "{}_{}_1.log".format(ver, mat)
                logp = os.path.join(ROOT, DIR, logf)
                with open(logp, "r") as log:
                    for line in log:
                        conv = line.split()[1]
                        mat_elp.append(float(line.split()[2]))
                        convnc = float(PREC0)
                        if float(conv) <= convnc:
                            break
                elp = 0
                if process == 1:
                    if ver == "cg_alg1":
                        belp = sum(mat_elp)
                elp = belp/sum(mat_elp)
                #if mat == "G2_circuit" and process == 16:
                #    print(ver, sum(mat_elp))

                if float(conv) <= CONV_CRIT:
                    PELAPSE.append(elp)
                else:
                    PELAPSE.append(np.inf)
            ELAPSE.append(PELAPSE)
        VELAPSE.append(ELAPSE)

    for p, process in enumerate(PROCESS):
        for v, ver in enumerate(VER):
            tmp = 0
            for m, mat in enumerate(MAT[:-2]):
                tmp += VELAPSE[m][v][p]
            tmp /= len(MAT[:-2])
            VELAPSE[-1][v][p] = tmp

    for imat, mat in enumerate(MAT):
        if mat == " ":
            continue
        for iver, ver in enumerate(VER):
            print(mat, ver)
            print(VELAPSE[imat][iver])
#    sys.exit(1)

    print("speedups")
    label_size = 12
    mpl.rcParams['xtick.labelsize'] = label_size 
    mpl.rcParams['ytick.labelsize'] = label_size 

#    fig = plt.figure(figsize=(4, 4))
#    gs1 = gridspec.GridSpec(4, 4)
#    gs1.update(wspace=0.325, hspace=0.35)
#    #fig,axes = plt.subplots(4, 3, sharex=False, sharey=False)
#    bar_width = 0.15
#    index = np.arange(len(PROCESS))
#    pname = "{}/cg_speedup_bar.pdf".format(ROOT)
#    with PdfPages(pname) as pdf:
#        for imat, mat in enumerate(MAT):
#            x = int(imat/3)
#            y = imat%3
#            ax1 = plt.subplot(4, 3, imat+1)#gs1[imat])
#            for iver, ver in enumerate(VER):
#                ax1.bar(index+bar_width*iver, VELAPSE[imat][iver], bar_width, alpha=1, label="{}".format(VER_NAME[iver]), color=COLOR[iver], align='edge')
#            ax1.set_aspect(0.3, adjustable='box')
#            ax1.set_title("{}".format(mat), fontsize=8)
#            ax1.set_xticks(index+bar_width)
#            ax1.set_xticklabels(PROCESS)
#            ax1.set_ylim([0,14])
#        #axes[-1][-1].axis('off')
#        #axes[-1][-2].axis('off')
#        #axes[-1][-3].axis('off')
#        fig.text(0.5, 0.11, 'Processes', ha='center', fontsize=10)
#        fig.text(0.005, 0.6, 'Speed-up', va='center', rotation='vertical', fontsize=10)
#        plt.subplots_adjust(wspace=0.11,hspace=0.01)
#        plt.legend(bbox_to_anchor=(0.7, 1.0), loc=1, ncol=6, borderaxespad=0., fontsize=9)
#        #plt.tight_layout()
#        pdf.savefig()
#        plt.close()


#    fig,axes = plt.subplots(4, 3, sharex=False, sharey=False)
#    bar_width = 0.13
#    index = np.arange(len(PROCESS))
#    pname = "{}/cg_speedup_bar.pdf".format(ROOT)
#    with PdfPages(pname) as pdf:
#        for imat, mat in enumerate(MAT):
#            x = int(imat/3)
#            y = imat%3
#            for iver, ver in enumerate(VER):
#                axes[x][y].bar(index+bar_width*iver, VELAPSE[imat][iver], bar_width, alpha=1, label="{}".format(VER_NAME[iver]), color=COLOR[iver], align='edge')
#            axes[x][y].autoscale(enable=False, axis='both', tight=None)
#            axes[x][y].set_aspect(0.3, adjustable='box')
#            #axes[x][y].set_aspect(1./axes[x][y].get_data_ratio(), adjustable='box')
#            axes[x][y].set_title("{}".format(mat), fontsize=8)
#            axes[x][y].set_xticks(index+bar_width)
#            axes[x][y].set_xticklabels(PROCESS)
#            axes[x][y].set_ylim([0,14])
#        axes[-1][-1].axis('off')
#        axes[-1][-2].axis('off')
#        axes[-1][-3].axis('off')
#        fig.text(0.5, 0.11, 'Processes', ha='center', fontsize=10)
#        fig.text(0.155, 0.6, 'Speed-up', va='center', rotation='vertical', fontsize=10)
#        plt.legend(bbox_to_anchor=(0.5, 0.5), loc=4, ncol=3, borderaxespad=0., fontsize=8)
#        plt.subplots_adjust(wspace=-0.45, hspace=0.40)
#        pdf.savefig(bbox_inches='tight', pad_inches=0.5)
#        plt.close()

    fig,axes = plt.subplots(2,5)
    fig.set_size_inches(18.5, 8.5)
    bar_width = 0.14
    index = np.arange(len(PROCESS))
    pname = "{}/cg_speedup_bar.pdf".format(ROOT)
    with PdfPages(pname) as pdf:
        for imat, mat in enumerate(MAT):
            x = int(imat/5)
            y = imat%5
            for iver, ver in enumerate(VER):
                axes[x][y].bar(index+bar_width*iver, VELAPSE[imat][iver], bar_width, alpha=1, label="{}".format(VER_NAME[iver]), color=COLOR[iver])
            axes[x][y].set_title("{}".format(mat), fontsize=12)
            axes[x][y].set_xticks(index+bar_width)
            axes[x][y].set_xticklabels(PROCESS)
            axes[x][y].set_ylim(ymin=0)
            axes[x][y].set_ylim([0,12])
        axes[-1][-2].axis('off')
        #fig.text(0.5, 0.04, 'Processes', ha='center', fontsize=10)
        #fig.text(0.04, 0.5, 'Speed-up', va='center', rotation='vertical', fontsize=10)
        #plt.legend(bbox_to_anchor=(0, 0), loc='upper left', ncol=1, fontsize=8)
        plt.legend(bbox_to_anchor=(-0.25, 1.0), loc=1, borderaxespad=0., fontsize=12)
        pdf.savefig(dpi=100)
        plt.close()

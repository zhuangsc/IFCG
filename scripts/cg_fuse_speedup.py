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
from matplotlib.backends.backend_pdf import PdfPages


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

#VER = ("cg_alg1", "cg_alg3", "cg_alg4", "cg_alg7", "cg_alg4_v2")
VER = ("cg_alg4_v2")#("cg_alg1", "cg_alg3", "cg_alg4", "cg_alg7", "cg_alg4_v2")
VER = ("cg_alg4_ifcg_v2")
#VER_NAME = ("PCG", "Chronopoulos", "Pipelined", "Gropp", "Pipelined V2")
#COLOR = ('r', 'g', 'b', 'k', 'y', 'm', 'c')
COLOR = ('0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1') #Grayscale
TASKS = 32
FUSE = 100
PROCESS = [1, 2, 4, 8, 16]
#PROCESS = [ p+1 for p in range(16) ]
FUSE = (1, 5, 20, 50, 80, 100, 200)

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
        for fuse in FUSE:
            PELAPSE = list()
            for process in PROCESS:
                if mat == " " or mat == "Mean":
                    PELAPSE.append(np.inf)
                    continue
                DIR = "{}_{}_{!r}_{}_{}".format(VER, TASKS, process, PREC0, fuse)
                conv = 0
                mat_elp = list()
                logf = "{}_{}.log".format(VER, mat)
                logp = os.path.join(ROOT, DIR, logf)
                with open(logp, "r") as log:
                    for line in log:
                        conv = line.split()[1]
                        mat_elp.append(float(line.split()[2]))
                        if float(conv) <= float(PREC):
                            break
                elp = 0
                if process == 1:
                    if fuse == 1:
                        belp = sum(mat_elp)
                elp = belp/sum(mat_elp)
                if process == 16:
                    print(mat, fuse, elp, belp, sum(mat_elp))

                if float(conv) <= CONV_CRIT:
                    PELAPSE.append(elp)
                else:
                    PELAPSE.append(np.inf)
            ELAPSE.append(PELAPSE)
        VELAPSE.append(ELAPSE)

    for p, process in enumerate(PROCESS):
        for f, fuse in enumerate(FUSE):
            tmp = 0
            for m, mat in enumerate(MAT[:-3]):
                tmp += VELAPSE[m][f][p]
            tmp /= len(MAT[:-3])
            VELAPSE[-1][f][p] = tmp

#    for imat, mat in enumerate(MAT):
#        if mat == "Mean":
#            for ifuse, fuse in enumerate(FUSE):
#                print(mat, fuse)
#                print(VELAPSE[imat][ifuse])
#    sys.exit(1)

    print("speedups")
    label_size = 9
    mpl.rcParams['xtick.labelsize'] = label_size 
    mpl.rcParams['ytick.labelsize'] = label_size 

#    fig,axes = plt.subplots(4, 3, sharex=False, sharey=False)
#    plt.tight_layout()
#    bar_width = 0.1
#    index = np.arange(len(PROCESS))
#    pname = "{}/cg_speedup_bar.pdf".format(ROOT)
#    with PdfPages(pname) as pdf:
#        for imat, mat in enumerate(MAT):
#            x = int(imat/3)
#            y = imat%3
#            for ifuse, fuse in enumerate(FUSE):
#                axes[x][y].bar(index+bar_width*ifuse, VELAPSE[imat][ifuse], bar_width, alpha=1, label="fuse {}".format(fuse), color=COLOR[ifuse])
#            axes[x][y].autoscale(enable=False, axis='both', tight=None)
#            axes[x][y].set_aspect(0.3, adjustable='box')
#            axes[x][y].set_title("{}".format(mat), fontsize=8)
#            axes[x][y].set_ylim([0,14])
#            axes[x][y].set_xticklabels(PROCESS)
#        axes[-1][-1].axis('off')
#        axes[-1][-2].axis('off')
#        axes[-1][-3].axis('off')
#        fig.text(0.5, 0.11, 'Processes', ha='center', fontsize=10)
#        fig.text(0.12, 0.6, 'Speed-up', va='center', rotation='vertical', fontsize=10)
#        plt.legend(bbox_to_anchor=(1.5, 1.0), loc=0, ncol=4, borderaxespad=0., fontsize=8)
#        plt.subplots_adjust(wspace=-0.45, hspace=0.40)
#        pdf.savefig(bbox_inches='tight', pad_inches=0.5)
#        plt.close()


    print("speedups")
    label_size = 12
    mpl.rcParams['xtick.labelsize'] = label_size 
    mpl.rcParams['ytick.labelsize'] = label_size 

    fig,axes = plt.subplots(2, 5, sharex=False, sharey=False)
    fig.set_size_inches(18.5, 8.5)
    bar_width = 0.12
    index = np.arange(len(PROCESS))
    pname = "{}/cg_speedup_bar.pdf".format(ROOT)
    with PdfPages(pname) as pdf:
        for imat, mat in enumerate(MAT):
            x = int(imat/5)
            y = imat%5
            for ifuse, fuse in enumerate(FUSE):
                axes[x][y].bar(index+bar_width*ifuse, VELAPSE[imat][ifuse], bar_width, alpha=1, label="fuse {}".format(fuse), color=COLOR[ifuse])
            axes[x][y].set_title("{}".format(mat), fontsize=12)
            axes[x][y].set_ylim([0,14])
            axes[x][y].set_xticklabels(PROCESS)
        axes[-1][-2].axis('off')
        #fig.text(0.5, 0.04, 'Processes', ha='center', fontsize=10)
        #fig.text(0.04, 0.5, 'Speed-up', va='center', rotation='vertical', fontsize=10)
        plt.legend(bbox_to_anchor=(-0.29, 1.0), loc=1, borderaxespad=0., fontsize=14)
        pdf.savefig()
        plt.close()

#    for imat, mat in enumerate(MAT):
#        pname = "{}/cg_speedup_{}.pdf".format(ROOT, mat)
#        with PdfPages(pname) as pdf:
#            for ifuse, fuse in enumerate(FUSE):
#                plt.plot(range(len(PROCESS)), VELAPSE[imat][ifuse], label="fuse {}".format(fuse), marker='o')
#            plt.title("{} TASKS {}".format(mat, TASKS))
#            plt.legend(loc=0)
#            plt.xticks(range(len(PROCESS)), PROCESS)
#            plt.ylim([1,PROCESS[-1]])
#            plt.xlabel("Processes")
#            plt.ylabel("Speed-up")
#            pdf.savefig()
#            plt.close()

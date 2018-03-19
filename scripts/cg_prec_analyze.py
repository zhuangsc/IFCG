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
import matplotlib.pyplot as plt
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
"Mean",
)

VER = ("cg_alg1", "cg_alg3", "cg_alg4", "cg_alg7", "cg_alg4_ifcg", "cg_alg4_ifcg_v2")
VER = ("cg_alg4_ifcg", "cg_alg4_ifcg_centinel", "cg_alg4_ifcg_v2", "cg_alg4_ifcg_v2_centinel")
COLOR = ('0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1') #Grayscale
#COLOR = ('r', 'g', 'b', 'k', 'y', 'm', 'c')
TASKS = 32
PROCESS = 16
FUSE = 20
VER_NAME = ("PCG", "Chronopoulos", "Pipelined", "Gropp", "IFCG", "IFCG2")
VER_NAME = ("ifcg_waiton", "ifcg_sentinel", "ifcg_v2_waiton", "ifcg_v2_sentinel")

PREC0 = "1E-8"
PREC = "1E-7"
CONV_CRIT = float(PREC) * 2

def ver_extract(ver, VELAPSE, VITER):
    ELAPSE = list()
    CONV = list()
    ITER = list()
    for mat in MAT:
        if mat == "Mean":
            ELAPSE.append(st.mean(ELAPSE))
            ITER.append(st.mean(ITER))
            continue

        DIR = "{}_{}_{}_{}_{}".format(ver, TASKS, PROCESS, PREC0, FUSE)
        conv = 0
        iters = 0
        mat_elp = list()
        fuse_iter = 0
        logf = "{}_{}_1.log".format(ver, mat)
        logp = os.path.join(ROOT, DIR, logf)
        with open(logp, "r") as log:
            for line in log:
                conv = float(line.split()[1])
                iters = int(line.split()[0])
                mat_elp.append(float(line.split()[2]))
#                if ver == "cg_alg4_at":
#                    fuse_iter += int(line.split()[3])
                convnc = float(PREC0)
                if mat == "nd24k":
                    convnc = float(PREC)
                if conv <= convnc:
                    break
        if conv <= CONV_CRIT:
            ELAPSE.append(sum(mat_elp))
            if ver == "cg_alg4_ifcg" or ver == "cg_alg4_ifcg_v2" or ver == "cg_alg4_ifcg_centinel" or ver == "cg_alg4_ifcg_v2_centinel":
                ITER.append(iters*20+20)
            elif ver == "cg_alg4_at":
                ITER.append(fuse_iter)
            else:
                ITER.append(iters+1)
        else:
            print("INF: {} {} {}".format(mat, PREC, ver))
            ELAPSE.append(np.inf)
            ITER.append(np.inf)
    VELAPSE.append(ELAPSE)
    VITER.append(ITER)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: {} [dir]".format(sys.argv[0]))
        sys.exit(1)
    ROOT = os.path.abspath(sys.argv[1])

    VELAPSE = list()
    VCONV = list()
    VITER = list()

    for ver in VER:
        ver_extract(ver, VELAPSE, VITER)

    MVITER = np.array(VITER)
    #print(MVITER.transpose())
    #print(MAT)
    #sys.exit(0)
    print("plotting")
    width = 0.1
    index = np.arange(len(MAT))
    pname = "{}/cg_elp_{}_{}_{}.pdf".format(ROOT, TASKS, PROCESS, PREC)
    with PdfPages(pname) as pdf:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for i, ver in enumerate(VER):
            ax.bar(index+width*i, VELAPSE[i], width=0.1, label=VER_NAME[i], color=COLOR[i])
        ax.set_ylabel(r'Elapse time ($\mu$s)')
        plt.title("TASKS {} PROCESS {} PREC {}".format(TASKS, PROCESS, PREC))
        ax.set_xticks(index+width)
        ax.set_xticklabels(MAT, rotation=45, fontsize=8)
        ax.legend(loc=0, fontsize=10)
        pdf.savefig()
        plt.close()

    pname = "{}/cg_iter_{}_{}_{}.pdf".format(ROOT, TASKS, PROCESS, PREC)
    with PdfPages(pname) as pdf:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for i, ver in enumerate(VER):
            ax.bar(index+width*i, VITER[i], width=0.1, label=VER_NAME[i], color=COLOR[i])
        ax.set_ylabel('Iterations')
        plt.title("TASKS {} PROCESS {} PREC {}".format(TASKS, PROCESS, PREC))
        ax.set_xticks(index+width)
        ax.set_xticklabels(MAT, rotation=45, fontsize=8)
        ax.legend(loc=0, fontsize=10)
        pdf.savefig()
        plt.close()

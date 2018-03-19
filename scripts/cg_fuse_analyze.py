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
"Mean"
)

VER = ("cg_alg4_v2")#("cg_alg1", "cg_alg3", "cg_alg4", "cg_alg7", "cg_alg4_v2")
#VER_NAME = ("", "Chronopoulos", "Pipelined", "Gropp", "V2")
#COLOR = ('r', 'g', 'b', 'k', 'y', 'm', 'c')
COLOR = ('0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1') #Grayscale
TASKS = 32
PROCESS = 16
FUSE = (1, 2, 5, 20, 50, 80, 100, 200)

PREC = "1E-8"
CONV_CRIT = float(PREC) * 2

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: {} [dir]".format(sys.argv[0]))
        sys.exit(1)
    ROOT = os.path.abspath(sys.argv[1])

    VELAPSE = list()
    VCONV = list()
    VITER = list()

    for fuse in FUSE:
        ELAPSE = list()
        CONV = list()
        ITER = list()
        for mat in MAT:
            if mat == "Mean":
                ELAPSE.append(st.mean(ELAPSE))
                ITER.append(st.mean(ITER))
                continue

            DIR = "{}_{}_{}_{}_{}".format(VER, TASKS, PROCESS, PREC, fuse)
            conv = 0
            iters = 0
            mat_elp = list()
            logf = "{}_{}.log".format(VER, mat)
            logp = os.path.join(ROOT, DIR, logf)
            with open(logp, "r") as log:
                for line in log:
                    conv = float(line.split()[1])
                    iters = int(line.split()[0])
                    mat_elp.append(float(line.split()[2]))
            if conv <= CONV_CRIT:
                ELAPSE.append(sum(mat_elp))
                ITER.append(iters*fuse)
            else:
                ELAPSE.append(np.inf)
                ITER.append(np.inf)
        VELAPSE.append(ELAPSE)
        VITER.append(ITER)


    print(MAT)
    print(VELAPSE)
    print("elp")
    width = 0.1
    index = np.arange(len(MAT))
    pname = "{}/cg_elp_{}_{}.pdf".format(ROOT, TASKS, PROCESS)
    with PdfPages(pname) as pdf:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for i, fuse in enumerate(FUSE):
            ax.bar(index+width*i, VELAPSE[i], width=0.1, label="fuse {}".format(fuse), color=COLOR[i])
        ax.set_ylabel(r'Elapse time ($\mu$s)')
        plt.title("TASKS {} PROCESS {}".format(TASKS, PROCESS))
        ax.set_xticks(index+width)
        ax.set_xticklabels(MAT, rotation=45, fontsize=8)
        ax.legend(loc=0)
        pdf.savefig()
        plt.close()

    print("iter")
    pname = "{}/cg_iter_{}_{}.pdf".format(ROOT, TASKS, PROCESS)
    with PdfPages(pname) as pdf:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for i, fuse in enumerate(FUSE):
            ax.bar(index+width*i, VITER[i], width=0.1, label="fuse {}".format(fuse), color=COLOR[i])
        ax.set_ylabel('Iterations')
        plt.title("TASKS {} PROCESS {}".format(TASKS, PROCESS))
        ax.set_xticks(index+width)
        ax.set_xticklabels(MAT, rotation=45, fontsize=8)
        ax.legend(loc=0)
        pdf.savefig()
        plt.close()

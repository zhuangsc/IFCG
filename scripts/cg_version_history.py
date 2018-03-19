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
)

#VER = ("cg_alg1", "cg_alg3", "cg_alg4", "cg_alg7", "cg_alg4_ifcg", "cg_alg4_ifcg_v2")
VER = ("cg_alg1", "cg_alg4", "cg_alg4_ifcg", "cg_alg4_ifcg_v2")
VER = ("cg_alg4_ifcg", "cg_alg4_ifcg_centinel", "cg_alg4_ifcg_v2", "cg_alg4_ifcg_v2_centinel")
COLOR = ('b', 'g', 'r', 'k', 'y', 'm', 'c')
#COLOR = ('0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1') #Grayscale
MARKER = ("o", ".", ",", "v", "^", "<", ">", "1", "2", "3", "4", "s")
TASKS = 32
PROCESS = 16
FUSE = 20
VER_NAME = ("PCG", "Pipelined", "IFCG", "IFCG2")
VER_NAME = ("ifcg_waiton", "ifcg_alltask", "ifcg_v2_waiton", "ifcg_v2_alltask")
#VER_NAME = ("PCG", "Chronopoulos", "Pipelined", "Gropp", "IFCG", "IFCG2")
TOT_ITER = 3500

PREC = "1E-16"
CONV_CRIT = float(PREC) * 2

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: {} [dir]".format(sys.argv[0]))
        sys.exit(1)
    ROOT = os.path.abspath(sys.argv[1])


    VELAPSE = list()
    VCONV = list()
    VITER = list()
    VHIST = list()

    for ver in VER:
        ELAPSE = list()
        CONV = list()
        ITER = list()
        HIST = list()
        for mat in MAT:
            DIR = "{}_{}_{}_{}_{}".format(ver, TASKS, PROCESS, PREC, FUSE)
            conv = 0
            iters = 0
            mat_elp = list()
            mat_res = list()
            logf = "{}_{}.log".format(ver, mat)
            logp = os.path.join(ROOT, DIR, logf)
            with open(logp, "r") as log:
                for line in log:
                    conv = float(line.split()[1])
                    iters = int(line.split()[0])
                    mat_elp.append(float(line.split()[2]))
                    if ver == "cg_alg4_ifcg" or ver == "cg_alg4_ifcg_v2" or ver == "cg_alg4_ifcg_centinel" or ver == "cg_alg4_ifcg_v2_centinel":
                        if iters*FUSE >= TOT_ITER:
                            break
                        for z in range(FUSE):
                            mat_res.append(math.log10(conv))
                    else:
                        if iters >= TOT_ITER:
                            break
                        mat_res.append(math.log10(conv))
            if ver == "cg_alg4_ifcg" or ver == "cg_alg4_ifcg_v2" or ver == "cg_alg4_ifcg_centinel" or ver == "cg_alg4_ifcg_v2_centinel":
                rep = iters*FUSE+FUSE
            else:
                rep = 1
            HIST.append(mat_res)
            ELAPSE.append(sum(mat_elp))
            ITER.append(len(mat_res))
        VELAPSE.append(ELAPSE)
        VITER.append(ITER)
        VHIST.append(HIST)

    for i, mat in enumerate(MAT):
        print("{}".format(mat))
        pname = "{}/cg_history_{}.pdf".format(ROOT, mat)
        with PdfPages(pname) as pdf:
            for j, ver in enumerate(VER_NAME):
                plt.plot(range(VITER[j][i]), VHIST[j][i], label=ver, color=COLOR[j])#, linestyle='--', )
            #plt.title('CG hist {} TASKS {} PROCESS {}'.format(mat, TASKS, PROCESS))
            plt.title(mat)
            plt.ylabel(r"$log_{10} ||b-Ax||_2/||b||_2$")
            plt.xlabel("Iterations")
            plt.legend(loc=0, fontsize=9)
            pdf.savefig()
            plt.close()

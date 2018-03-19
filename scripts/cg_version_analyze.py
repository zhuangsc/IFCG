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

VER = ("cg_alg4_bp", "cg_alg4_ifcg_bp", "cg_alg4_ifcg_v2_bp", "cg_alg4", "cg_alg4_ifcg", "cg_alg4_ifcg_v2")
#VER = ("cg_alg4_ifcg_bp", "cg_alg4_ifcg_v2_bp", "cg_alg4_bp")
COLOR = ('0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1') #Grayscale
#COLOR = ('r', 'g', 'b', 'k', 'y', 'm', 'c')
TASKS = 32
PROCESS = 16
FUSE = 100
VER_NAME = ("Pipelined CG BP", "IFCG BP", "IFCG2 BP", "Pipelined CG", "IFCG", "IFCG2")

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
        logf = "{}_{}.log".format(ver, mat)
        logp = os.path.join(ROOT, DIR, logf)
        with open(logp, "r") as log:
            for line in log:
                conv = float(line.split()[1])
                iters = int(line.split()[0])
                mat_elp.append(float(line.split()[2]))
                if ver == "cg_alg4_at":
                    fuse_iter += int(line.split()[3])
                if conv <= float(PREC):
                    break

        if conv <= CONV_CRIT:
            if ver == "cg_alg4_bp":
                belp = sum(mat_elp)
            elp = sum(mat_elp)/belp
            ELAPSE.append(elp)
            if ver == "cg_alg4_ifcg" or ver == "cg_alg4_ifcg_v2" or ver == "cg_alg4_ifcg_bp" or ver == "cg_alg4_ifcg_v2_bp":
                ITER.append(iters*FUSE+FUSE)
            elif ver == "cg_alg4_at":
                ITER.append(fuse_iter)
            else:
                ITER.append(iters+1)
        else:
            print("INF: {} {}".format(mat, PREC))
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
            logf = "{}_{}.log".format(ver, mat)
            logp = os.path.join(ROOT, DIR, logf)
            with open(logp, "r") as log:
                for line in log:
                    conv = float(line.split()[1])
                    iters = int(line.split()[0])
                    mat_elp.append(float(line.split()[2]))
                    if ver == "cg_alg4_at":
                        fuse_iter += int(line.split()[3])
                    if conv <= float(PREC):
                        break

            if conv <= CONV_CRIT:
                if ver == "cg_alg4_bp":
                    belp = sum(mat_elp)
                elp = sum(mat_elp)/belp
                #ELAPSE.append(elp)
                ELAPSE.append(sum(mat_elp))
                if ver == "cg_alg4_ifcg" or ver == "cg_alg4_ifcg_v2" or ver == "cg_alg4_ifcg_bp" or ver == "cg_alg4_ifcg_v2_bp":
                    ITER.append(iters*FUSE+FUSE)
                elif ver == "cg_alg4_at":
                    ITER.append(fuse_iter)
                else:
                    ITER.append(iters+1)
            else:
                print("INF: {} {}".format(mat, PREC))
                ELAPSE.append(np.inf)
                ITER.append(np.inf)
        VELAPSE.append(ELAPSE)
        VITER.append(ITER)

        #ver_extract(ver, VELAPSE, VITER)

    baseline = VELAPSE[0][:]
    for iver, ver in enumerate(VER):
        for i in range(len(VELAPSE[iver])):
            VELAPSE[iver][i] /= baseline[i]
        print(ver)
        print(VELAPSE[iver])

    #sys.exit(1)
    print("plotting")
    width = 0.1
    index = np.arange(len(MAT))
    pname = "{}/cg_elp_{}_{}_{}.pdf".format(ROOT, TASKS, PROCESS, PREC)
    with PdfPages(pname) as pdf:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for i, ver in enumerate(VER):
            ax.bar(index+width*i, VELAPSE[i], width=0.1, label=VER_NAME[i], color=COLOR[i])
        #ax.set_ylabel(r'Elapse time ($\mu$s)')
        #plt.title("TASKS {} PROCESS {} PREC {}".format(TASKS, PROCESS, PREC))
        ax.set_xticks(index+width)
        ax.set_xticklabels(MAT, rotation=45, fontsize=8)
        ax.set_ylim([0,1.5])
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
        #plt.title("TASKS {} PROCESS {} PREC {}".format(TASKS, PROCESS, PREC))
        ax.set_xticks(index+width)
        ax.set_xticklabels(MAT, rotation=45, fontsize=8)
        ax.legend(loc=0, fontsize=10)
        pdf.savefig()
        plt.close()

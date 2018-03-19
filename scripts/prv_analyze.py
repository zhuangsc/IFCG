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


#! /usr/bin/env python3

import sys, os, re, shutil

WD = "9200002" ### Work Descriptor ID event
ITER = 20
RTT = 17

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("USAGE: {} [prv]".format(sys.argv[0]))
        sys.exit(1)
    PRV = os.path.abspath(sys.argv[1])
    BASE_DIR = os.path.dirname(PRV)
    PRV_BASE = os.path.basename(PRV)

    BPRV = "fuse_" + PRV_BASE
    BPRV = os.path.join(BASE_DIR, BPRV)

    COM_ROOT = os.path.splitext(PRV_BASE)[0]
    ROW = COM_ROOT + ".row"
    BROW = "fuse_" + ROW
    BROW = os.path.join(BASE_DIR, BROW)
    PCF = COM_ROOT + ".pcf"
    BPCF = "fuse_" + PCF
    BPCF = os.path.join(BASE_DIR, BPCF)

    max_id = 0
    with open(PRV, "r") as prv:
        for line in prv:
            match = re.search("{}:[0-9]+".format(WD), line)
            if match:
                #print(match.group(0))
                wd_id = int(match.group(0).split(":")[1]) ### use group(0) since same state would not appear twice in a line
                if max_id < wd_id:
                    max_id = wd_id
    TITER = ( max_id - RTT ) / ITER
    print("max_id {} RTT {}".format(max_id, RTT))
    print("tasks per iter: {}".format(TITER))

    with open(PRV, "r") as prv:
        with open(BPRV, "w") as bprv:
            for line in prv:
                match = re.search("{}:[0-9]+".format(WD), line)
                if match:
                    idx0 = int(match.group(0).split(":")[1])
                    idx1 = (idx0-RTT)/TITER
                    idx1 = int(idx1)
#                    print("idx0: {} idx1: {}".format(idx0, idx1))
                    line = re.sub("{}:[0-9]+".format(WD), "{}:{}".format(WD, idx1), line)
#                    print(line)
                bprv.write(line)

    shutil.copy2(ROW, BROW)
    shutil.copy2(PCF, BPCF)

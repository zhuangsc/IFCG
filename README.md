Source code of the Iteration-Fusing Conjugate Gradient (IFCG)

How to run:

   ./dcg [bm] [it] [precision] [correction] [iter_fuse] [rep] [orth_fac] [HB_MAT] [FULL?] [CG_VER] [cglog] [B]

    bm: block size

    it: number of iterations

    precision: target convergence (e.g. 1E-8)

    correction: frequency of residual correction (e.g. every 50 iterations)

    iter_fuse: FUSE parameter

    rep: execution repetition

    orth_fac: 0, no effect

    HB_PATH: path to the matrix file (.rb)

    FULL: 1

    CG_VER: the version of CG (specified in cg_main.c)

    cglog (optional): 0/1 toggle for logfile output

    B (optional): path to a right hand side file. If not specified, randomly generates a B vector and writes to a file RHS.dat

Reference:

Sicong Zhuang and Marc Casas. 2017. Iteration-fusing conjugate gradient. In Proceedings of the International Conference on Supercomputing (ICS '17). 
ACM, New York, NY, USA, Article 21, 10 pages. DOI: https://doi.org/10.1145/3079079.3079091

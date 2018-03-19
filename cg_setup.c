//MIT License
//
//Copyright (c) 2018 Sicong Zhuang
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.


#include "cg_setup.h"

#include <stdlib.h>
#include <stdio.h>

extern void *Ahbh;
extern void *Ahb;
extern void *Acsr;
extern void *preconditioner;
extern css **S;
extern csn **N;
extern cs *Acs;
extern char *rhsfname;
extern char *aname;

extern css *SA;
extern csn *NA;
extern cs *AA;

void gen_cg_rhs(int m, int n, double *b)
{
	unsigned int seed = 1591613054;
	srand(seed);
	double range = (double) 1;

	int i;
	for ( i = 0; i < m*n; i++ ) {
		b[i] = ((double)rand() / (double)RAND_MAX ) * range - range/2;
	}
}

int cg_setup(int n, int bm, double **x, double **rhs) 
{
	AA = malloc(sizeof(cs));
	SA = malloc(sizeof(css));
	NA = malloc(sizeof(csn));
	hbmat_t *MA = malloc(sizeof(hbmat_t));
	hb_sym_diag(Acsr, bm, MA);
	AA->nz = -1;
	AA->nzmax = MA->elemc;
	AA->n = AA->m = MA->m;
	AA->p = MA->vptr;
	AA->i = MA->vpos;
	AA->x = MA->vval;
	SA = cs_schol(AA, 0);
	NA = cs_chol(AA, SA);

	Ahbh = hb2hbh(Acsr, bm, 1);
	hbmat_t *lAhb = Acsr;
	int elemc = lAhb->elemc;
	int bs = (n+bm-1)/bm;
	preconditioner = malloc(bs * sizeof(hbmat_t));
	hb_sym_diag_block(Acsr, bm, preconditioner);
	hbmat_t *precond = preconditioner;
	Acs = malloc(bs*sizeof(cs));
	S = malloc(bs*sizeof(css*));
	N = malloc(bs*sizeof(csn*));
	for (int i = 0; i < bs; i++ ) {
		hbmat_t *Apre = &precond[i];
		cs *Acspre = &Acs[i];
		Acspre->nz = -1;
		Acspre->nzmax = Apre->elemc;
		Acspre->n = Acspre->m = Apre->m;
		Acspre->p = Apre->vptr;
		Acspre->i = Apre->vpos;
		Acspre->x = Apre->vval;
		S[i] = cs_schol(Acspre, 0);
//		N[i] = cs_chol_ilu(Acspre, S[i]);
		N[i] = cs_chol(Acspre, S[i]);
	}
	for (int i = 0; i < bs; i++) {
		free(precond[i].vptr);
		free(precond[i].vpos);
		free(precond[i].vval);
	}
	free(precond);

	*x = malloc(n * sizeof(double));
//	double *lx = malloc(2 * n * sizeof(double));
//	if ( lx == NULL ) {
//		fprintf(stderr, "err: could not allocate x\n");
//		return 1;
//	}
//	x[0] = &lx[0];
//	x[1] = &lx[n];

	/* Construct or read the rhs (B) */
	
	double *lrhs;
	if ( rhsfname != NULL ) {
		printf("read rhs\n");
		lrhs = *rhs = (double*) calloc(n, sizeof(double));
		if ( lrhs == NULL ) {
			fprintf(stderr, "err: could not allocate rhs\n");
			return 3;
		}

		FILE *rhsf;
		if ( rhsf = fopen(rhsfname, "r") ) {
			read_mm2dense(rhsf, n, 1, lrhs);
			fclose(rhsf);
		} else {
			fprintf(stderr, "setup: could not open %s\n", rhsfname);
			return 6;
		}
	} else {
		printf("gen rhs\n");

		lrhs = *rhs = malloc(n*sizeof(double));
		gen_cg_rhs(n, 1, lrhs);
		fprint_dense2mm("RHS.dat", "RHS", n, 1, lrhs, n);
	}

	return 0;

}

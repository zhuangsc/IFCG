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


#ifndef __CG_AUX_H__
#define __CG_AUX_H__

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "vector.h"
#include "csparse.h"
#include "hb_io.h"

#define FP_SQRT 	sqrt
#define FP_RAND 	drand48
#define FP_SEED 	srand48
#define FP_ABS 		fabs	
#define FP_EXP 		frexp
#define FP_LOG10 	log10
#define FP_POW 		pow

#define FP_SCANSPEC	scan_dconspec

#ifdef INTEL_MKL

#include "mkl.h"
#define BLAS_cp(n, dx, incx, dy, incy) 						cblas_dcopy(n, dx, incx, dy, incy)
#define BLAS_dot(n, dx, incx, dy, incy) 					cblas_ddot(n, dx, incx, dy, incy) 
#define BLAS_axpy(n, da, dx, incx, dy, incy) 				cblas_daxpy(n, da, dx, incx, dy, incy)
#define SBLAS_csrmv(trans, m, n, alpha, matdescra, avval, avpos, avptr, avptr1, Bptr, beta, Cptr) \
	mkl_dcsrmv(trans, &m, &n, &alpha, matdescra, avval, avpos, avptr, avptr1, Bptr, &beta, Cptr)

#elif defined LAPACK

#include "cblas.h"
#define BLAS_cp(n, dx, incx, dy, incy) 						cblas_dcopy(n, dx, incx, dy, incy)
#define BLAS_dot(n, dx, incx, dy, incy) 					cblas_ddot(n, dx, incx, dy, incy) 
#define BLAS_axpy(n, da, dx, incx, dy, incy) 				cblas_daxpy(n, da, dx, incx, dy, incy)
#define SBLAS_csrmv(trans, m, n, alpha, matdescra, avval, avpos, avptr, avptr1, Bptr, beta, Cptr) \
	manual_csrmv(trans, m, n, alpha, avval, avpos, avptr, Bptr, beta, Cptr)

#endif


typedef struct strhbmat {
	int m, n;
	int elemc;
	int *vptr;
	int *vpos;
	void *vval;
	int *vdiag;
	int *udiagc;
	int b;
	int type;
	struct strhbmat *orig;
	struct strhbmat *trans;
	struct strhbmat *hyper;
	int orig_row;
	int orig_col;
	int *e_tree;
	int FACT;

	/*
	 * The following for hyper-matrix only
	 */
	int *vptr_pool;
	int *vpos_pool;
	void *vval_pool;
	int vptr_unit;
	int vpos_unit;
	int vval_unit;
	int vptr_pp;
	int vpos_pp;
	int vval_pp;
//	pthread_mutex_t* mtx;

} hbmat_t;

extern const char *scan_dconspec;
extern const char *scan_sconspec;


void hb_read_double(char *input_file, int *m, int *n, int *elemc, int **vptr, int **vpos, double **vval);
void hb_reset(hbmat_t *A);
void one2zero(hbmat_t* in_matrix);

void hb_sym_expand(hbmat_t *A, hbmat_t *B);

void hb_init_basic(hbmat_t *A, hbmat_t *B);

void hb_free(hbmat_t *A);

void* __hb2hbh_block(int I, int J, hbmat_t *A, int b, hbmat_t *Bp) ;

hbmat_t* hb2hbh(hbmat_t *A, int b, int is_csr);

void hb_sym_diag_block(hbmat_t *src_mat, int bsze, hbmat_t *diagb);

int read_mm2dense(FILE *f, int m, int n, double *A);

void print_dense2mm(FILE *f, const char *name, int m, int n, const double *A, int lda);

void fprint_dense2mm(const char *fname, const char *name, int m, int n, const double *A, int lda);


static inline void __attribute__((always_inline)) bblas_dcopy(int p, int bm, int bn, int m, int n, double *X, double *Y) 
{
	int i;
	for ( i=0; i<m; i+=bm ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			__t_copy(p, c, d, m, n, X, Y, j*m+i, j*m+i);
//			__t_copy(p, c, d, m, n, &X[j*m+i], &Y[j*m+i]);
		}
	}
}

static inline void __attribute__((always_inline)) hbsblas_dcsrmv(int p, int b, double alpha, hbmat_t *Ahbh, double *B, double beta, double *C) 
{
	int M = Ahbh->m;
	int N = Ahbh->n;
	int *vptr = Ahbh->vptr;
	int *vpos = Ahbh->vpos;
	hbmat_t **vval = Ahbh->vval;
	int offs = vptr[0] == 0 ? 0 : 1; //Detect zero/one based

	int cmaj = 1;
	char *trans = "N";
	char *matdescra = "GLNC";
	double fp_one = 1.0;

	int I;
	for ( I = 0; I < M; ++I ) {
		double *Cptr = &C[I*b];
		int first = 1;
		int J;
		for ( J = vptr[I]; J < vptr[I+1]; J++ ) {
			hbmat_t *A = vval[J];
			int icol = vpos[J];
			double *Bptr = &B[icol*b];
			double *avval = A->vval;
			int *avpos = A->vpos;
			int *avptr = A->vptr;
			int m = A->m;
			int n = A->n;
			if ( first ) {
				#pragma omp task in(B[icol*b:icol*b+n-1]) out(C[I*b:I*b+m-1]) no_copy_deps label(csrmv_hbh) priority(p)
				SBLAS_csrmv(trans, m, n, alpha, matdescra, avval, avpos, avptr, avptr+1, Bptr, beta, Cptr);
//				mkl_dcsrmv(trans, &m, &n, &alpha, matdescra, avval, avpos, avptr, avptr+1, Bptr, &beta, Cptr);
				first = 0;
			} else {
				#pragma omp task in(B[icol*b:icol*b+n-1]) out(C[I*b:I*b+m-1]) no_copy_deps label(csrmv_hbh) priority(p)
//				#pragma omp task in(B[icol*b;n]) out(C[I*b;m]) no_copy_deps label(csrmv_hbh) priority(p)
				SBLAS_csrmv(trans, m, n, alpha, matdescra, avval, avpos, avptr, avptr+1, Bptr, fp_one, Cptr);
//				mkl_dcsrmv(trans, &m, &n, &alpha, matdescra, avval, avpos, avptr, avptr+1, Bptr, &fp_one, Cptr);
			}
		}
	}
}

static inline void __attribute__((always_inline)) bsblas_dcholsolv2(int p, int b, int m, css **S, csn **N, double *B, double *x)
{
	int idx;
	int i;
	for ( i = 0, idx = 0; i < m; i+=b, idx++) {
		int bs = b < m-i ? b : m-i;
		css *sptr = S[idx];
		csn *nptr = N[idx];
		double *bptr = &B[i];
		double *xptr = &x[i];
		#pragma omp task in(B[i:i+bs-1]) out(x[i:i+bs-1]) label(dcholsolv2)
		cs_cholsol2(bs, sptr, nptr, bptr, xptr);
	}
}

static inline void __attribute__((always_inline)) bsblas_dcholsolv2_seq(int p, int b, int m, css **S, csn **N, double *B, double *x)
{
	int idx;
	int i;
	for ( i = 0, idx = 0; i < m; i+=b, idx++) {
		int bs = b < m-i ? b : m-i;
		css *sptr = S[idx];
		csn *nptr = N[idx];
		double *bptr = &B[i];
		double *xptr = &x[i];
		cs_cholsol2(bs, sptr, nptr, bptr, xptr);
	}
}

static inline void __attribute__((always_inline)) dcholsolv2_blk(int p, int m, css *S, csn *N, double *B, double *x)
{
	#pragma omp task in(B[0:m-1]) out(x[0:m-1]) label(dcholsolv2_blk)
	cs_cholsol2(m, S, N, B, x);
}

static inline void __attribute__((always_inline)) dcholsolv2_nested(int p, int b, int m, css **S, csn **N, double *B, double *x)
{
	#pragma omp task in(B[0:m-1]) out(x[0:m-1]) label(dcholsolv2_nested)
	{
		int idx;
		int i;
		for ( i = 0, idx = 0; i < m; i+=b, idx++) {
			int bs = b < m-i ? b : m-i;
			css *sptr = S[idx];
			csn *nptr = N[idx];
			double *bptr = &B[i];
			double *xptr = &x[i];
			#pragma omp task label(dcholsolv2_in) //in([bs]bptr) out([bs]xptr) label(dcholsolv2_in)
			cs_cholsol2(bs, sptr, nptr, bptr, xptr);
		}
		#pragma omp taskwait
	}
}

static inline void __attribute__((always_inline)) bblas_ddot(int p, int bm, int bn, int m, int n, double *X, double *Y, double *result) 
{
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int idx;
		int i;
		for ( i=0, idx=0; i<m; i+=bm, ++idx ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;

			__t_dot(p, c, d, m, n, X, Y, j*m+i, j*m+i, result);
//			__t_dot(p, c, d, m, n, &X[j*m+i], &Y[j*m+i], result);
		}
		result += bn;
	}
}

static inline void __attribute__((always_inline)) bblas_dcpaxpy_comb(int bm, int bn, int m, int n, double alpha, double *Anum, double *Aden, double *X1, double *X2, double *Y1, double *Y2, double *Z1, double *Z2)
{
	int i;
	for ( i=0; i<m; i+=bm ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			__t_cpaxpy_comb(c, d, m, n, alpha, &Anum[j], &Aden[j], &X1[j*m+i], &X2[j*m+i], &Y1[j*m+i], &Y2[j*m+i], &Z1[j*m+i], &Z2[j*m+i]);
		}
	}
}

static inline void __attribute__((always_inline)) bblas_extm_daxpy(int p, int bm, int bn, int m, int n, double *SAnum, double *SAden, double *X, double *Y, double *Z) 
{
	int i;
	for ( i=0; i<m; i+=bm ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			__t_extm_axpy(c, d, m, n, &SAnum[j], &SAden[j], &X[j*m+i], &Y[j*m+i], &Z[j*m+i], p);
		}
	}
}

static inline __attribute__((always_inline)) void cg_ddot2(int p, int bm, int bn, int m, int n, double *X, double *Y, double *result, double *A, double *B, double *result2) 
{
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int idx;
		int i;
		for ( i=0, idx=0; i<m; i+=bm, ++idx ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;
			_cg_dot2(p, c, d, m, n, X, Y, j*m+i, j*m+i, result, A, B, j*m+i, j*m+i, result2);
		}
		result += bn;
		result2 += bn;
	}
}

static inline void __attribute__((always_inline)) bblas_ddot_array(int p, int bm, int bn, int m, int n, double *X, double *Y, double *result) 
{
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int idx;
		int i;
		for ( i=0, idx=0; i<m; i+=bm, ++idx ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;
			__t_dot_array(p, c, d, m, n, X, Y, j*m+i, j*m+i, result, idx);
		}
		result += bn;
	}
}

static inline __attribute__((always_inline)) void cg_ddot2_array(int p, int bm, int bn, int m, int n, double *X, double *Y, double *result, double *A, double *B, double *result2) 
{
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int idx;
		int i;
		for (i=0, idx=0; i<m; i+=bm, idx++) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;
			_cg_dot2_array(p, c, d, m, n, X, Y, j*m+i, j*m+i, result, idx, A, B, j*m+i, j*m+i, result2, idx);
		}
		result += bn;
		result2 += bn;
	}
}


#pragma omp task in(X[initx:initx+bm-1]) out(Y[inity:inity+bm-1]) priority(p) label(dcopy) no_copy_deps
void __t_copy(int p, int bm, int bn, int m, int n, double *X, double *Y, int initx, int inity);

#pragma omp task in(X[initx:initx+bm-1], Y[inity:inity+bm-1]) concurrent(result[0:bn-1]) no_copy_deps priority(p) label(ddot)
void __t_dot(int p, int bm, int bn, int m, int n, double *X, double *Y, int initx, int inity, double *result);

#pragma omp task in(X[initx:initx+bm-1], Y[inity:inity+bm-1]) concurrent(result[0:(m+bm-1)/bm-1]) no_copy_deps priority(p) label(ddot_array)
void __t_dot_array(int p, int bm, int bn, int m, int n, double *X, double *Y, int initx, int inity, double *result, int initr);

#pragma omp task in(X1[0:bm-1], X2[0:bm-1], Anum[0:bn-1], Aden[0:bn-1], Y1[0:bm-1], Y2[0:bm-1]) out(Z1[0:bm-1], Z2[0:bm-1]) no_copy_deps priority(1) label(dcpaxpy_comb)
void __t_cpaxpy_comb(int bm, int bn, int m, int n, double alpha, double *Anum, double *Aden, double *X1, double *X2, double *Y1, double *Y2, double *Z1, double *Z2);

#pragma omp task in(X[0:bm-1], Y[0:bm-1], SAnum[0:bn-1], SAden[0:bn-1]) out(Z[0:bm-1]) no_copy_deps priority(p) label(extm_axpy)
void __t_extm_axpy(int bm, int bn, int m, int n, double *SAnum, double *SAden, double *X, double *Y, double *Z, int p);

#pragma omp task in(X[initx:initx+bm-1], Y[inity:inity+bm-1], A[inita:inita+bm-1], B[initb:initb+bm-1]) concurrent([bn]result, [bn]result2) no_copy_deps priority(p) label(cg_dot2)
void _cg_dot2(int p, int bm, int bn, int m, int n, double *X, double *Y, int initx, int inity, double *result, double *A, double *B, int inita, int initb, double *result2);

#pragma omp task in(X[initx:initx+bm-1], Y[inity:inity+bm-1], A[inita:inita+bm-1], B[initb:initb+bm-1]) concurrent(result[0:(m+bm-1)/bm-1], result2[0:(m+bm-1)/bm-1]) no_copy_deps priority(p) label(ddot2_array)
void _cg_dot2_array(int p, int bm, int bn, int m, int n, double *X, double *Y, int initx, int inity, double *result, int initr, double *A, double *B, int inita, int initb, double *result2, int initr2);


void manual_csrmv(char *trans, int m, int n, double alpha, double *avval, int *avpos, int *avptr, double *Bptr, double beta, double *Cptr);

#endif //__CG_AUX_H__

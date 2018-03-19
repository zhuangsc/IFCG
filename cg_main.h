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


#ifndef __CG_MAIN_H__
#define __CG_MAIN_H__

#include <float.h>
#include "cg_config.h"
#include "cg_setup.h"
#include "cg_aux.h"

struct timeval start, stop;

static inline __attribute__((always_inline)) void dump_info(char *name, int k, double *residuals, unsigned int *elapse)
{
	FILE *log = fopen(name, "w");
	for ( int i = 0; i <= k; i++ ) {
		fprintf(log, "%d %E %u\n", i, residuals[i], elapse[i]);
	}
	fclose(log);
}

static inline __attribute__((always_inline)) void dump_info3(char *name, int k, double *residuals, unsigned int *elapse, int *fuses)
{
	FILE *log = fopen(name, "w");
	for ( int i = 0; i <= k; i++ ) {
		fprintf(log, "%d %E %u %d\n", i, residuals[i], elapse[i], fuses[i]);
	}
	fclose(log);
}

static inline __attribute__((always_inline)) void dump_info4(char *name, int k, double *residuals, unsigned int *elapse, int *fuses, int *iters)
{
	FILE *log = fopen(name, "w");
	for ( int i = 0; i <= k; i++ ) {
		fprintf(log, "%d %E %u %d %d\n", i, residuals[i], elapse[i], fuses[i], iters[i]);
	}
	fclose(log);
}


static inline __attribute__((always_inline)) void start_timer()
{
	gettimeofday(&start, NULL);
}

static inline __attribute__((always_inline)) void stop_timer(unsigned int *elp)
{
	gettimeofday(&stop, NULL);
	*elp = (stop.tv_sec - start.tv_sec) * 1e6 + stop.tv_usec - start.tv_usec;
}



/* Commutative implementation */
#pragma omp task in([bm]X, [bm]Y, [bm]A, [bm]B) commutative([bn]result, [bn]result2) no_copy_deps priority(p) label(cg_dot2)
void _cg_dot2_commutative(int p, int bm, int bn, int m, int n, double *X, double *Y, double *result, double *A, double *B, double *result2) 
{
	double fp_one = 1.0;
	int i_one = 1;
	double local_result[bn];
	for ( int j=0; j<bn; ++j ) {
		local_result[j] = BLAS_dot(bm, X, i_one, Y, i_one);
		X += m;
		Y += m;
	}

	double local_result2[bn];
	int j;
	for ( int j=0; j<bn; ++j ) {
		local_result2[j] = BLAS_dot(bm, A, i_one, B, i_one);
		A += m;
		B += m;
	}

	BLAS_axpy(bn, fp_one, local_result, i_one, result, i_one);
	BLAS_axpy(bn, fp_one, local_result2, i_one, result2, i_one);
}


static inline __attribute__((always_inline)) void cg_ddot2_commutative(int p, int bm, int bn, int m, int n, double *X, double *Y, double *result, double *A, double *B, double *result2) 
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

			_cg_dot2_commutative(p, c, d, m, n, &X[j*m+i], &Y[j*m+i], result, &A[j*m+i], &B[j*m+i], result2);
		}
		result += bn;
		result2 += bn;
	}
}


#endif //__CG_MAIN_H__

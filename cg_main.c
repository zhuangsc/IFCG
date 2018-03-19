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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

#include "cg_main.h"
#include "csparse.h"

void *Ahbh;
void *Ahb;
void *Acsr;
int n;
int bm;
int cgit;
double prec;
int correction;
int iter_fuse;
//int is_precond;
int rep;
char *aname;
char *rhsfname;
double orth_fac;
int cglog_level;
int cg_ver;

double *rhs;
double *x;

void *preconditioner;
css **S;
csn **N;
cs *Acs;

css *SA;
csn *NA;
cs *AA;

double fp_one = 1.0;
double fp_mone = -1.0;
double fp_nought = 0.0;

int main(int argc, char *argv[])
{
	if ( cg_config(argc, argv) ) {
		return 1;
	}

	if ( cg_setup(n, bm, &x, &rhs) ) {
		return 2;
	}

	hbmat_t *Atmp = Acsr;
	int dim = Atmp->m;
	double *pool = malloc(20*dim*sizeof(double));
	#pragma omp register ([20*dim]pool)

	for ( int i = 0; i < rep; i++ ) {
		memset(pool, 0, 20*dim*sizeof(double));
		switch ( cg_ver ) {
			case 0: 
				printf("ALG1 PCG\n");
				CG_ALG1(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, orth_fac, cglog_level);  //Algorithm 1 original PCG
				break;
			case 1:
				printf("ALG3 Chronopoulos\n");
				CG_ALG3(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, orth_fac, cglog_level);  //Algorithm 3 Chronopoulos
				break;
			case 2:
				printf("ALG4 Pipelined\n");
				CG_ALG4(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, orth_fac, cglog_level);  //Algorithm 4 Pipelined
				break;
			case 3:
				printf("ALG7 Gropp\n");
				CG_ALG7(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, orth_fac, cglog_level);  //Algorithm 7 Gropp
				break;
			case 4:
				printf("ALG4 IFCG\n");
				CG_IFCG(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, iter_fuse, orth_fac, cglog_level);  //Algorithm 4 IFCG
				break;
			case 5:
				printf("ALG4 IFCG V2\n");
				CG_IFCG_V2(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, iter_fuse, orth_fac, cglog_level);  //Algorithm 4 IFCG V2
				break;
			case 6:
				printf("ALG4 IFCG CENTINEL\n");
				CG_IFCG_Centinel(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, iter_fuse, orth_fac, cglog_level);  //Algorithm 4 IFCG Centinel
				break;
			case 7:
				printf("ALG4 IFCG V2 CENTINEL\n");
				CG_IFCG_V2_Centinel(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, iter_fuse, orth_fac, cglog_level);  //Algorithm 4 IFCG V2 Centinel
				break;
			//TODO BPs
			case 10:
				printf("ALG4 IFCG BP\n");
				CG_IFCG_BP(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, iter_fuse, orth_fac, cglog_level, pool);  //Algorithm 4 IFCG BP
				break;
			case 11:
				printf("ALG4 IFCG V2 BP\n");
				CG_IFCG_V2_BP(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, iter_fuse, orth_fac, cglog_level, pool);  //Algorithm 4 IFCG V2 BP
				break;
			case 12: 
				printf("ALG1 PCG BP\n");
				CG_ALG1_BP(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, orth_fac, cglog_level, pool);  //Algorithm 1 original PCG BP
				break;
			case 13:
				printf("ALG3 Chronopoulos BP\n");
				CG_ALG3_BP(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, orth_fac, cglog_level, pool);  //Algorithm 3 Chronopoulos BP
				break;
			case 14:
				printf("ALG4 Pipelined BP\n");
				CG_ALG4_BP(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, orth_fac, cglog_level, pool);  //Algorithm 4 Pipelined BP
				break;
			case 15:
				printf("ALG7 Gropp BP\n");
				CG_ALG7_BP(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, orth_fac, cglog_level, pool);  //Algorithm 7 Gropp BP
				break;
			//TODO ILUs
			case 20:
				printf("ALG4 IFCG ILU BP\n");
				CG_IFCG_ILU_BP(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, iter_fuse, orth_fac, cglog_level, pool);  //Algorithm 4 IFCG BP
				break;
			case 21:
				printf("ALG4 IFCG V2 ILU BP\n");
				CG_IFCG_V2_ILU_BP(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, iter_fuse, orth_fac, cglog_level, pool);  //Algorithm 4 IFCG V2 BP
				break;
			case 22: 
				printf("ALG1 PCG ILU BP\n");
				CG_ALG1_ILU_BP(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, orth_fac, cglog_level, pool);  //Algorithm 1 original PCG BP
				break;
			case 23:
				printf("ALG3 Chronopoulos ILU BP\n");
				CG_ALG3_ILU_BP(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, orth_fac, cglog_level, pool);  //Algorithm 3 Chronopoulos BP
				break;
			case 24:
				printf("ALG4 Pipelined ILU BP\n");
				CG_ALG4_ILU_BP(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, orth_fac, cglog_level, pool);  //Algorithm 4 Pipelined BP
				break;
			case 25:
				printf("ALG7 Gropp ILU BP\n");
				CG_ALG7_ILU_BP(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, orth_fac, cglog_level, pool);  //Algorithm 7 Gropp BP
				break;

			case 80:
				printf("ALG4 IFCG Commutative\n"); //Algorithm 4 IFCG Commutative (commutative dot-product)
				CG_ALG4_V4(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, iter_fuse, orth_fac, cglog_level);  
				break;
			case 81:
				printf("ALG4 IFCG-AT\n"); //Algorithm 4 IFCG-AT
				CG_ALG4_AT(Acsr, Ahbh, x, rhs, cgit, bm, prec, correction, iter_fuse, orth_fac, cglog_level);  
				break;
			default:
				printf("No algorithm selected\n");
				break;
		}
	}
	return 0;
}

/* Standard PCG */
int CG_ALG1(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, double orth_fac, int cglog_level)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
	double *pool = calloc(4 * 2 * n, sizeof(double));
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha * p_i
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u = M^-1 * r_i+1
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta * p_i
	double  *s = malloc(n * sizeof(double)); // s = Ap

	double *alpha1 = calloc(2, sizeof(double));
	double *alpha2 = calloc(2, sizeof(double));

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
	bsblas_dcholsolv2(1, bm, n, S, N, r[i], u[i]);
	/* p[0] = z[0] */
	bblas_dcopy(1, bm, 1, n, 1, u[i], p[i]);
	/* alpha1[0] = <r[0], z[0]> */
	bblas_ddot(1, bm, 1, n, 1, r[i], u[i], &alpha1[i]);

	int k;
	for ( k = 0; k < cgit; k++ ) {
		start_timer();
		int iprev = i;
		i = i ^ 0x1;

		/* s = A * p[i] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, p[iprev], fp_nought, s); 
		/* alpha2[i] = <s, p[i]> */
		bblas_ddot(1, bm, 1, n, 1, s, p[iprev], &alpha2[i]);
		bblas_dcpaxpy_comb(bm, 1, n, 1, fp_mone, &alpha1[iprev], &alpha2[i], s, p[iprev], r[iprev], x[iprev], r[i], x[i]);

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		/* u[i+1] = M^-1 * r[i+1] */
		bsblas_dcholsolv2(1, bm, n, S, N, r[i], u[i]);
		/* alpha1[i+1] = <r, u> */
		bblas_ddot(1, bm, 1, n, 1, r[i], u[i], &alpha1[i]);
		/* p[i+1] = u[i+1] + transpose(beta[i]) * p[i] */
		bblas_extm_daxpy(1, bm, 1, n, 1, &alpha1[i], &alpha1[iprev], p[iprev], u[i], p[i]); 	

		#pragma omp taskwait

		stop_timer(&elapses[k]);
		alpha1[iprev] = alpha2[iprev] = (double) 0;
//		BLAS_gemm(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, 1, 1, n, FP_ONE, p[i], n, s, n, FP_NOUGHT, &orth, 1);
//		orth = FP_ABS(orth);
//		if (isgreater(orth, porth * orth_fac)){
//			fprintf(stderr, "orth fail %E %E\n", orth, porth*orth_fac);
//			break;
//		}

		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		residuals[k] = sr2norm;
		if ( isless(sr2norm, prec) ) {
			fprintf(stderr, "Precision reached\n");
			break;
		}
//		fprintf(stdout, "%d %E\n", k, residuals[k]);
	}
	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg1.log", k, residuals, elapses);
	free(pool);
	free(alpha1);
	free(alpha2);
	free(residuals);
	free(elapses);
	return 0;
}

/* Chronopoulos PCG */
int CG_ALG3(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, double orth_fac, int cglog_level)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
	double *pool = calloc(6 * 2 * n, sizeof(double));
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha * s_i
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u = M^-1 * r_i+1
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w = A * u_i+1
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta * p_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i+1 = w_i+1 + beta * s_i

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;

	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
	bsblas_dcholsolv2(1, bm, n, S, N, r[i], u[i]);
	/* w[0] = A * u[0] */
	hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]);
	/* gamma[0] = <r[0], u[0]> */
	bblas_ddot(1, bm, 1, n, 1, r[i], u[i], &gamma[i]);
	/* alpha[0] = gamma[0]/<w[0], u[0]> */
	bblas_ddot(1, bm, 1, n, 1, w[i], u[i], &alpha[i]);

	#pragma omp taskwait on (gamma[i], alpha[i])

	alpha[i] = gamma[i] / alpha[i];

	int k;
	for ( k = 0; k < cgit; k++ ) {
		start_timer();
		int iprev = i;
		i = i ^ 0x1;

		/* Grand fuse */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *pp0 = &(p[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *wp0 = &(w[iprev])[j];
			double *xp0 = &(x[iprev])[j];
			double *rp0 = &(r[iprev])[j];

			double *pp1 = &(p[i])[j];
			double *up1 = &(u[i])[j];
			double *sp1 = &(s[i])[j];
			double *wp1 = &(w[i])[j];
			double *xp1 = &(x[i])[j];
			double *rp1 = &(r[i])[j];

			#pragma omp task out([c]pp1, [c]sp1, [c]xp1, [c]rp1)  priority(1) label(alg3_fuse)
			{
				/* p_i = u_i + beta_i * p_i-1 */
				BLAS_cp(c, up0, 1, pp1, 1);
				BLAS_axpy(c, beta[iprev], pp0, 1, pp1, 1);

				/* s_i = w_i + beta_i * s_i-1 */
				BLAS_cp(c, wp0, 1, sp1, 1);
				BLAS_axpy(c, beta[iprev], sp0, 1, sp1, 1);

				/* x_i+1 = x_i + alpha_i * p_i */
				BLAS_cp(c, xp0, 1, xp1, 1);
				BLAS_axpy(c, alpha[iprev], pp1, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				BLAS_cp(c, rp0, 1, rp1, 1);
				BLAS_axpy(c, -1*alpha[iprev], sp1, 1, rp1, 1);
			}
		}

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		/* u[i+1] = M^-1 * r[i+1] */
		bsblas_dcholsolv2(1, bm, n, S, N, r[i], u[i]);
		/* w[i+1] = A * u[i+1] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]); 

		cg_ddot2(1, bm, 1, n, 1, r[i], u[i], &gamma[i], w[i], u[i], &delta);

		#pragma omp taskwait

		beta[i] = gamma[i] / gamma[iprev];
		alpha[i] = gamma[i]/(delta - beta[i] * gamma[i] / alpha[iprev]);

		stop_timer(&elapses[k]);
		gamma[iprev] = delta = (double) 0;

		//TODO Implement p-orthogonality check
#if 0
		BLAS_gemm(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, 1, 1, n, FP_ONE, p[i], n, s, n, FP_NOUGHT, &orth, 1);
		orth = FP_ABS(orth);
		if (isgreater(orth, porth * orth_fac)){
			fprintf(stderr, "orth fail %E %E\n", orth, porth*orth_fac);
			break;
		}
#endif

		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		residuals[k] = sr2norm;
		if ( isless(sr2norm, prec) ) {
			fprintf(stderr, "Precision reached\n");
			break;
		}
//		fprintf(stdout, "%d %E\n", k, residuals[k]);
	}

	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg3.log", k, residuals, elapses);
	free(pool);
	free(alpha);
	free(beta);
	free(gamma);
	free(residuals);
	free(elapses);
	return 0;
}

/* Pipelined PCG */
int CG_ALG4(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, double orth_fac, int cglog_level)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
	double *pool = calloc(10 * 2 * n, sizeof(double));
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u_i+1 = u_i - alpha_i * q_i
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w_i+1 = w_i - alpha_i * z_i
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta_i * p_i
	offs += 2 * n;
	double *m[2] = {&pool[offs], &pool[offs+n]}; // m_i = M^-1 * w_i
	offs += 2 * n;
	double *n0[2] = {&pool[offs], &pool[offs+n]}; // n_i = A * m_i
	offs += 2 * n;
	double *z[2] = {&pool[offs], &pool[offs+n]}; // z_i = n_i + beta_i * z_i-1
	offs += 2 * n;
	double *q[2] = {&pool[offs], &pool[offs+n]}; // q_i = m_i + beta_i * q_i-1
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha_i * s_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i = w_i + beta_i * s_i-1

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
	bsblas_dcholsolv2(1, bm, n, S, N, r[i], u[i]);
	/* w[0] = A * u[0] */
	hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]);

	int k;
	for ( k = 0; k < cgit; k++ ) {
		start_timer();
		int iprev = i;
		i = i ^ 0x1;

		/* m[i] = M^-1 * w[i] */
		bsblas_dcholsolv2(1, bm, n, S, N, w[iprev], m[i]);
		/* n[i] = A * m[i] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, m[i], fp_nought, n0[i]); 
		/* 
		 * gamma[i] = <r[i], u[i]>
		 * delta = <w[i], u[i]>
		 */
		cg_ddot2(1, bm, 1, n, 1, r[iprev], u[iprev], &gamma[i], w[iprev], u[iprev], &delta);

		#pragma omp taskwait on(gamma[i], delta)
		if ( k > 0 ) {
			beta[i] = gamma[i]/gamma[iprev];
			alpha[i] = gamma[i] / (delta - beta[i] * gamma[i] / alpha[iprev]);
		} else {
			beta[i] = (double) 0;
			alpha[i] = gamma[i]/delta;
		}
		gamma[iprev] = delta = 0;

		/* Grand fuse */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *zp0 = &(z[iprev])[j];
			double *qp0 = &(q[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *xp0 = &(x[iprev])[j];
			double *rp0 = &(r[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *zp1 = &(z[i])[j];
			double *qp1 = &(q[i])[j];
			double *sp1 = &(s[i])[j];
			double *pp1 = &(p[i])[j];
			double *xp1 = &(x[i])[j];
			double *rp1 = &(r[i])[j];
			double *up1 = &(u[i])[j];
			double *wp1 = &(w[i])[j];

			double *mp1 = &(m[i])[j];
			double *np1 = &(n0[i])[j];

			#pragma omp task out([c]zp1, [c]qp1, [c]sp1, [c]pp1, [c]xp1, [c]rp1, [c]up1, [c]wp1) \
				in([c]zp0, [c]qp0, [c]sp0, [c]pp0, [c]np1, [c]mp1, [c]up0, [c]xp0, [c]wp0, [c]rp0) \
				priority(1) label(alg4_fuse)
			{
				/* z_i = n_i + beta_i * z_i-1 */
				BLAS_cp(c, np1, 1, zp1, 1);
				BLAS_axpy(c, beta[i], zp0, 1, zp1, 1);

				/* q_i = m_i + beta_i * q_i-1 */
				BLAS_cp(c, mp1, 1, qp1, 1);
				BLAS_axpy(c, beta[i], qp0, 1, qp1, 1);

				/* s_i = w_i + beta_i * s_i-1 */
				BLAS_cp(c, wp0, 1, sp1, 1);
				BLAS_axpy(c, beta[i], sp0, 1, sp1, 1);

				/* p_i = u_i + beta_i * p_i-1 */
				BLAS_cp(c, up0, 1, pp1, 1);
				BLAS_axpy(c, beta[i], pp0, 1, pp1, 1);

				/* x_i+1 = x_i + alpha_i * p_i */
				BLAS_cp(c, xp0, 1, xp1, 1);
				BLAS_axpy(c, alpha[i], pp1, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				BLAS_cp(c, rp0, 1, rp1, 1);
				BLAS_axpy(c, -1*alpha[i], sp1, 1, rp1, 1);

				/* u_i+1 = u_i - alpha_i * q_i */
				BLAS_cp(c, up0, 1, up1, 1);
				BLAS_axpy(c, -1*alpha[i], qp1, 1, up1, 1);

				/* w_i+1 = w_i - alpha_i * z_i */
				BLAS_cp(c, wp0, 1, wp1, 1);
				BLAS_axpy(c, -1*alpha[i], zp1, 1, wp1, 1);
			}
		}

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		#pragma omp taskwait

//		gamma[iprev] = delta = 0;
		stop_timer(&elapses[k]);

		//TODO Implement p-orthogonality check
#if 0
		BLAS_gemm(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, 1, 1, n, FP_ONE, p[i], n, s, n, FP_NOUGHT, &orth, 1);
		orth = FP_ABS(orth);
		if (isgreater(orth, porth * orth_fac)){
			fprintf(stderr, "orth fail %E %E\n", orth, porth*orth_fac);
			break;
		}
#endif

		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		residuals[k] = sr2norm;
		if ( isless(sr2norm, prec) ) {
			fprintf(stderr, "Precision reached\n");
			break;
		}
//		fprintf(stdout, "%d %E\n", k, residuals[k]);
	}

	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg4.log", k, residuals, elapses);
	free(pool);
	free(alpha);
	free(beta);
	free(gamma);
	free(residuals);
	free(elapses);
	return 0;
}

/* Gropp PCG */
int CG_ALG7(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, double orth_fac, int cglog_level)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
	double *pool = calloc(7 * 2 * n, sizeof(double));
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u_i+1 = u_i - alpha_i * q_i
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w_i+1 = w_i - alpha_i * z_i
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta_i * p_i
	offs += 2 * n;
	double *q[2] = {&pool[offs], &pool[offs+n]}; // q_i = m_i + beta_i * q_i-1
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha_i * s_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i = w_i + beta_i * s_i-1

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
	bsblas_dcholsolv2(1, bm, n, S, N, r[i], u[i]);
	/* p[0] = u[0] */
	bblas_dcopy(1, bm, 1, n, 1, u[i], p[i]);
	/* s[0] = A * p[0] */
	hbsblas_dcsrmv(1, bm, fp_one, Ahbh, p[i], fp_nought, s[i]);
	/* gamma[0] = <r[0], u[0]> */
	bblas_ddot(1, bm, 1, n, 1, r[i], u[i], &gamma[i]);

	int k;
	for ( k = 0; k < cgit; k++ ) {
		start_timer();
		int iprev = i;
		i = i ^ 0x1;

		/* delta = <p[i], s[i]> */
		bblas_ddot(1, bm, 1, n, 1, p[iprev], s[iprev], &delta);
		/* q[i] = M^-1 * s[i] */
		bsblas_dcholsolv2(1, bm, n, S, N, s[iprev], q[i]);

		#pragma omp taskwait on(delta)
		alpha[i] = gamma[iprev]/delta;

		/* Axpy fuse x,r,u */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *xp0 = &(x[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *rp0 = &(r[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *qp0 = &(q[iprev])[j];

			double *xp1 = &(x[i])[j];
			double *pp1 = &(p[i])[j];
			double *rp1 = &(r[i])[j];
			double *sp1 = &(s[i])[j];
			double *up1 = &(u[i])[j];
			double *qp1 = &(q[i])[j];

			#pragma omp task in([c]qp1) out([c]xp1, [c]rp1, [c]up1) priority(1) label(alg7_fuse0)
			{
				/* x_i+1 = x_i + alpha_i * p_i */
				BLAS_cp(c, xp0, 1, xp1, 1);
				BLAS_axpy(c, alpha[i], pp0, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				BLAS_cp(c, rp0, 1, rp1, 1);
				BLAS_axpy(c, -1*alpha[i], sp0, 1, rp1, 1);

				/* u_i+1 = u_i - alpha_i * q_i */
				BLAS_cp(c, up0, 1, up1, 1);
				BLAS_axpy(c, -1*alpha[i], qp1, 1, up1, 1);
			}
		}

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait			
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		/* gamma[i+1] = <r[i+1], u[i+1]> */
		bblas_ddot(1, bm, 1, n, 1, r[i], u[i], &gamma[i]);
		/* w[i+1] = A * u[i+1] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]);

		#pragma omp taskwait on(gamma[i])
		beta[i] = gamma[i]/gamma[iprev];

		/* Axpy fuse p,s */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *pp0 = &(p[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *pp1 = &(p[i])[j];
			double *sp1 = &(s[i])[j];
			double *up1 = &(u[i])[j];
			double *wp1 = &(w[i])[j];

			#pragma omp task in([c]up1, [c]pp0, [c]wp1, [c]sp0) priority(1) label(alg7_fuse1)
			{
				/* p_i = u_i + beta_i * p_i-1 */
				BLAS_cp(c, up1, 1, pp1, 1);
				BLAS_axpy(c, beta[i], pp0, 1, pp1, 1);

				/* s_i = w_i + beta_i * s_i-1 */
				BLAS_cp(c, wp1, 1, sp1, 1);
				BLAS_axpy(c, beta[i], sp0, 1, sp1, 1);
			}
		}

		#pragma omp taskwait

		gamma[iprev] = delta = 0;
		stop_timer(&elapses[k]);

		//TODO Implement p-orthogonality check
#if 0
		BLAS_gemm(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, 1, 1, n, FP_ONE, p[i], n, s, n, FP_NOUGHT, &orth, 1);
		orth = FP_ABS(orth);
		if (isgreater(orth, porth * orth_fac)){
			fprintf(stderr, "orth fail %E %E\n", orth, porth*orth_fac);
			break;
		}
#endif

		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		residuals[k] = sr2norm;
		if ( isless(sr2norm, prec) ) {
			fprintf(stderr, "Precision reached\n");
			break;
		}
//		fprintf(stdout, "%d %E\n", k, residuals[k]);
	}

	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg7.log", k, residuals, elapses);
	free(pool);
	free(alpha);
	free(beta);
	free(gamma);
	free(residuals);
	free(elapses);
	return 0;
}

/* IFCG */
int CG_IFCG(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, int fuse, double orth_fac, int cglog_level)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
	double *pool = calloc(10 * 2 * n, sizeof(double));
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u_i+1 = u_i - alpha_i * q_i
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w_i+1 = w_i - alpha_i * z_i
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta_i * p_i
	offs += 2 * n;
	double *m[2] = {&pool[offs], &pool[offs+n]}; // m_i = M^-1 * w_i
	offs += 2 * n;
	double *n0[2] = {&pool[offs], &pool[offs+n]}; // n_i = A * m_i
	offs += 2 * n;
	double *z[2] = {&pool[offs], &pool[offs+n]}; // z_i = n_i + beta_i * z_i-1
	offs += 2 * n;
	double *q[2] = {&pool[offs], &pool[offs+n]}; // q_i = m_i + beta_i * q_i-1
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha_i * s_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i = w_i + beta_i * s_i-1

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
	bsblas_dcholsolv2(1, bm, n, S, N, r[i], u[i]);
	/* w[0] = A * u[0] */
	hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]);

	start_timer();
	int kk = 0;
	int k;
	for ( k = 0; k < cgit; k++) {

		int iprev = i;
		i = i ^ 0x1;

		/* m[i] = M^-1 * w[i] */
		bsblas_dcholsolv2(1, bm, n, S, N, w[iprev], m[i]);
		/* n[i] = A * m[i] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, m[i], fp_nought, n0[i]); 
		/* 
		 * gamma[i] = <r[i], u[i]>
		 * delta = <w[i], u[i]>
		 */
		cg_ddot2(1, bm, 1, n, 1, r[iprev], u[iprev], &gamma[i], w[iprev], u[iprev], &delta);

		#pragma omp taskwait on(gamma[i], delta)
		if ( k > 0 ) {
			beta[i] = gamma[i]/gamma[iprev];
			alpha[i] = gamma[i] / (delta - beta[i] * gamma[i] / alpha[iprev]);
		} else {
			beta[i] = (double) 0;
			alpha[i] = gamma[i]/delta;
		}
		gamma[iprev] = delta = 0;

		/* Grand fuse */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *zp0 = &(z[iprev])[j];
			double *qp0 = &(q[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *xp0 = &(x[iprev])[j];
			double *rp0 = &(r[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *zp1 = &(z[i])[j];
			double *qp1 = &(q[i])[j];
			double *sp1 = &(s[i])[j];
			double *pp1 = &(p[i])[j];
			double *xp1 = &(x[i])[j];
			double *rp1 = &(r[i])[j];
			double *up1 = &(u[i])[j];
			double *wp1 = &(w[i])[j];

			double *mp1 = &(m[i])[j];
			double *np1 = &(n0[i])[j];

			#pragma omp task out([c]zp1, [c]qp1, [c]sp1, [c]pp1, [c]xp1, [c]rp1, [c]up1, [c]wp1) \
				in([c]zp0, [c]qp0, [c]sp0, [c]pp0, [c]np1, [c]mp1, [c]up0, [c]xp0, [c]wp0, [c]rp0) \
				priority(1) label(alg4_fuse)
			{
				/* z_i = n_i + beta_i * z_i-1 */
				BLAS_cp(c, np1, 1, zp1, 1);
				BLAS_axpy(c, beta[i], zp0, 1, zp1, 1);

				/* q_i = m_i + beta_i * q_i-1 */
				BLAS_cp(c, mp1, 1, qp1, 1);
				BLAS_axpy(c, beta[i], qp0, 1, qp1, 1);

				/* s_i = w_i + beta_i * s_i-1 */
				BLAS_cp(c, wp0, 1, sp1, 1);
				BLAS_axpy(c, beta[i], sp0, 1, sp1, 1);

				/* p_i = u_i + beta_i * p_i-1 */
				BLAS_cp(c, up0, 1, pp1, 1);
				BLAS_axpy(c, beta[i], pp0, 1, pp1, 1);

				/* x_i+1 = x_i + alpha_i * p_i */
				BLAS_cp(c, xp0, 1, xp1, 1);
				BLAS_axpy(c, alpha[i], pp1, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				BLAS_cp(c, rp0, 1, rp1, 1);
				BLAS_axpy(c, -1*alpha[i], sp1, 1, rp1, 1);

				/* u_i+1 = u_i - alpha_i * q_i */
				BLAS_cp(c, up0, 1, up1, 1);
				BLAS_axpy(c, -1*alpha[i], qp1, 1, up1, 1);

				/* w_i+1 = w_i - alpha_i * z_i */
				BLAS_cp(c, wp0, 1, wp1, 1);
				BLAS_axpy(c, -1*alpha[i], zp1, 1, wp1, 1);
			}
		}

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		if ( k > 0 && k % fuse == 0 ) {
			#pragma omp taskwait
			stop_timer(&elapses[kk]);
			double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
			double sr2norm = norm_r/norm_b;
			residuals[kk] = sr2norm;
			if ( isless(sr2norm, prec) ) {
				fprintf(stderr, "Precision reached\n");
				break;
			}
			kk += 1;
			start_timer();
		}
	}

	#pragma omp taskwait
	if ( k == cgit ) {
		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		stop_timer(&elapses[kk]);
		residuals[kk] = sr2norm;
	}

	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg4_ifcg.log", kk, residuals, elapses);
	free(pool);
	free(alpha);
	free(beta);
	free(gamma);
	free(residuals);
	free(elapses);
	return 0;
}

/* IFCG Centinel */
int CG_IFCG_Centinel(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, int fuse, double orth_fac, int cglog_level)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
	double *pool = calloc(10 * 2 * n, sizeof(double));
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u_i+1 = u_i - alpha_i * q_i
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w_i+1 = w_i - alpha_i * z_i
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta_i * p_i
	offs += 2 * n;
	double *m[2] = {&pool[offs], &pool[offs+n]}; // m_i = M^-1 * w_i
	offs += 2 * n;
	double *n0[2] = {&pool[offs], &pool[offs+n]}; // n_i = A * m_i
	offs += 2 * n;
	double *z[2] = {&pool[offs], &pool[offs+n]}; // z_i = n_i + beta_i * z_i-1
	offs += 2 * n;
	double *q[2] = {&pool[offs], &pool[offs+n]}; // q_i = m_i + beta_i * q_i-1
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha_i * s_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i = w_i + beta_i * s_i-1

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	int bs = (n+bm-1)/bm;
	double *a_gamma = malloc(bs * sizeof(double));
	double *a_delta = malloc(bs * sizeof(double));

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
	bsblas_dcholsolv2(1, bm, n, S, N, r[i], u[i]);
	/* w[0] = A * u[0] */
	hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]);

	start_timer();
	int kk = 0;
	int k;
	for ( k = 0; k < cgit; k++) {

		int iprev = i;
		i = i ^ 0x1;

		/* m[i] = M^-1 * w[i] */
		bsblas_dcholsolv2(1, bm, n, S, N, w[iprev], m[i]);
		/* n[i] = A * m[i] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, m[i], fp_nought, n0[i]); 
		/* 
		 * gamma[i] = <r[i], u[i]>
		 * delta = <w[i], u[i]>
		 */
//		cg_ddot2_array(1, bm, 1, n, 1, r[iprev], u[iprev], a_gamma, w[iprev], u[iprev], a_delta);
		cg_ddot2(1, bm, 1, n, 1, r[iprev], u[iprev], &gamma[i], w[iprev], u[iprev], &delta);

		#pragma omp task in(gamma[i], delta) out(beta[i], alpha[i]) label(centinel)
		{
			if ( k > 0 ) {
				beta[i] = gamma[i]/gamma[iprev];
				alpha[i] = gamma[i] / (delta - beta[i] * gamma[i] / alpha[iprev]);
			} else {
				beta[i] = (double) 0;
				alpha[i] = gamma[i]/delta;
			}
			gamma[iprev] = delta = 0;
		}
#if 0
		#pragma omp task in(a_gamma[0:bs-1], a_delta[0:bs-1]) out(beta[i], alpha[i]) label(centinel)
		{
			gamma[i] = (double) 0;
			delta = (double) 0;
			for (int ii = 0; ii < bs; ii++) {
				gamma[i] += a_gamma[ii];
				delta += a_delta[ii];
			}
			if ( k > 0 ) {
				beta[i] = gamma[i]/gamma[iprev];
				alpha[i] = gamma[i] / (delta - beta[i] * gamma[i] / alpha[iprev]);
			} else {
				beta[i] = (double) 0;
				alpha[i] = gamma[i]/delta;
			}
		}
#endif

		/* Grand fuse */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *zp0 = &(z[iprev])[j];
			double *qp0 = &(q[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *xp0 = &(x[iprev])[j];
			double *rp0 = &(r[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *zp1 = &(z[i])[j];
			double *qp1 = &(q[i])[j];
			double *sp1 = &(s[i])[j];
			double *pp1 = &(p[i])[j];
			double *xp1 = &(x[i])[j];
			double *rp1 = &(r[i])[j];
			double *up1 = &(u[i])[j];
			double *wp1 = &(w[i])[j];

			double *mp1 = &(m[i])[j];
			double *np1 = &(n0[i])[j];

			#pragma omp task out([c]zp1, [c]qp1, [c]sp1, [c]pp1, [c]xp1, [c]rp1, [c]up1, [c]wp1) \
				in(alpha[i], beta[i], [c]zp0, [c]qp0, [c]sp0, [c]pp0, [c]np1, [c]mp1, [c]up0, [c]xp0, [c]wp0, [c]rp0) \
				priority(1) label(alg4_fuse)
			{
				/* z_i = n_i + beta_i * z_i-1 */
				BLAS_cp(c, np1, 1, zp1, 1);
				BLAS_axpy(c, beta[i], zp0, 1, zp1, 1);

				/* q_i = m_i + beta_i * q_i-1 */
				BLAS_cp(c, mp1, 1, qp1, 1);
				BLAS_axpy(c, beta[i], qp0, 1, qp1, 1);

				/* s_i = w_i + beta_i * s_i-1 */
				BLAS_cp(c, wp0, 1, sp1, 1);
				BLAS_axpy(c, beta[i], sp0, 1, sp1, 1);

				/* p_i = u_i + beta_i * p_i-1 */
				BLAS_cp(c, up0, 1, pp1, 1);
				BLAS_axpy(c, beta[i], pp0, 1, pp1, 1);

				/* x_i+1 = x_i + alpha_i * p_i */
				BLAS_cp(c, xp0, 1, xp1, 1);
				BLAS_axpy(c, alpha[i], pp1, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				BLAS_cp(c, rp0, 1, rp1, 1);
				BLAS_axpy(c, -1*alpha[i], sp1, 1, rp1, 1);

				/* u_i+1 = u_i - alpha_i * q_i */
				BLAS_cp(c, up0, 1, up1, 1);
				BLAS_axpy(c, -1*alpha[i], qp1, 1, up1, 1);

				/* w_i+1 = w_i - alpha_i * z_i */
				BLAS_cp(c, wp0, 1, wp1, 1);
				BLAS_axpy(c, -1*alpha[i], zp1, 1, wp1, 1);
			}
		}

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		if ( k > 0 && k % fuse == 0 ) {
			#pragma omp taskwait
			stop_timer(&elapses[kk]);
			double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
			double sr2norm = norm_r/norm_b;
			residuals[kk] = sr2norm;
			if ( isless(sr2norm, prec) ) {
				fprintf(stderr, "Precision reached\n");
				break;
			}
			kk += 1;
			start_timer();
		}
	}

	#pragma omp taskwait
	if ( k == cgit ) {
		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		stop_timer(&elapses[kk]);
		residuals[kk] = sr2norm;
	}

	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg4_ifcg_centinel.log", kk, residuals, elapses);
	free(pool);
	free(alpha);
	free(beta);
	free(gamma);
	free(a_gamma);
	free(a_delta);
	free(residuals);
	free(elapses);
	return 0;
}

/* IFCG V2 */
int CG_IFCG_V2(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, int fuse, double orth_fac, int cglog_level)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
	double *pool = calloc(10 * 2 * n, sizeof(double));
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u_i+1 = u_i - alpha_i * q_i
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w_i+1 = w_i - alpha_i * z_i
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta_i * p_i
	offs += 2 * n;
	double *m[2] = {&pool[offs], &pool[offs+n]}; // m_i = M^-1 * w_i
	offs += 2 * n;
	double *n0[2] = {&pool[offs], &pool[offs+n]}; // n_i = A * m_i
	offs += 2 * n;
	double *z[2] = {&pool[offs], &pool[offs+n]}; // z_i = n_i + beta_i * z_i-1
	offs += 2 * n;
	double *q[2] = {&pool[offs], &pool[offs+n]}; // q_i = m_i + beta_i * q_i-1
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha_i * s_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i = w_i + beta_i * s_i-1

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
	bsblas_dcholsolv2(1, bm, n, S, N, r[i], u[i]);
	/* w[0] = A * u[0] */
	hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]);

	start_timer();
	int kk = 0;
	int k;
	for ( k = 0; k < cgit; k++) {

		int iprev = i;
		i = i ^ 0x1;

		/* m[i] = M^-1 * w[i] */
		bsblas_dcholsolv2(1, bm, n, S, N, w[iprev], m[i]);
		/* 
		 * gamma[i] = <r[i], u[i]>
		 */
		bblas_ddot(1, bm, 1, n, 1, r[iprev], u[iprev], &gamma[i]);
		bblas_ddot(1, bm, 1, n, 1, w[iprev], u[iprev], &delta);

		#pragma omp taskwait on(gamma[i])
		if ( k > 0 ) {
			beta[i] = gamma[i]/gamma[iprev];
		} else {
			beta[i] = (double) 0;
		}

		/* 
		 * s_i = w_i + beta_i * s_i-1
		 * p_i = u_i + beta_i * p_i-1
		 */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *sp0 = &(s[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *sp1 = &(s[i])[j];
			double *pp1 = &(p[i])[j];

			#pragma omp task out([c]sp1, [c]pp1) in(beta[i], [c]sp0, [c]pp0, [c]up0, [c]wp0) priority(1) label(alg4_apxy)
			{
				/* s_i = w_i + beta_i * s_i-1 */
				BLAS_cp(c, wp0, 1, sp1, 1);
				BLAS_axpy(c, beta[i], sp0, 1, sp1, 1);

				/* p_i = u_i + beta_i * p_i-1 */
				BLAS_cp(c, up0, 1, pp1, 1);
				BLAS_axpy(c, beta[i], pp0, 1, pp1, 1);
			}
		}

		/* n[i] = A * m[i] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, m[i], fp_nought, n0[i]); 

		#pragma omp taskwait on(delta)
		if ( k > 0 ) {
			alpha[i] = gamma[i] / (delta - beta[i] * gamma[i] / alpha[iprev]);
		} else {
			alpha[i] = gamma[i]/delta;
		}
		gamma[iprev] = delta = 0;

		/* Grand fuse */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *zp0 = &(z[iprev])[j];
			double *qp0 = &(q[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *xp0 = &(x[iprev])[j];
			double *rp0 = &(r[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *zp1 = &(z[i])[j];
			double *qp1 = &(q[i])[j];
			double *sp1 = &(s[i])[j];
			double *pp1 = &(p[i])[j];
			double *xp1 = &(x[i])[j];
			double *rp1 = &(r[i])[j];
			double *up1 = &(u[i])[j];
			double *wp1 = &(w[i])[j];

			double *mp1 = &(m[i])[j];
			double *np1 = &(n0[i])[j];

			#pragma omp task out([c]zp1, [c]qp1, [c]xp1, [c]rp1, [c]up1, [c]wp1) \
				in(alpha[i], [c]zp0, [c]qp0, [c]sp1, [c]pp1, [c]np1, [c]mp1, [c]up0, [c]xp0, [c]wp0, [c]rp0) \
				priority(1) label(alg4_fuse)
			{
				/* z_i = n_i + beta_i * z_i-1 */
				BLAS_cp(c, np1, 1, zp1, 1);
				BLAS_axpy(c, beta[i], zp0, 1, zp1, 1);

				/* q_i = m_i + beta_i * q_i-1 */
				BLAS_cp(c, mp1, 1, qp1, 1);
				BLAS_axpy(c, beta[i], qp0, 1, qp1, 1);

				/* x_i+1 = x_i + alpha_i * p_i */
				BLAS_cp(c, xp0, 1, xp1, 1);
				BLAS_axpy(c, alpha[i], pp1, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				BLAS_cp(c, rp0, 1, rp1, 1);
				BLAS_axpy(c, -1*alpha[i], sp1, 1, rp1, 1);

				/* u_i+1 = u_i - alpha_i * q_i */
				BLAS_cp(c, up0, 1, up1, 1);
				BLAS_axpy(c, -1*alpha[i], qp1, 1, up1, 1);

				/* w_i+1 = w_i - alpha_i * z_i */
				BLAS_cp(c, wp0, 1, wp1, 1);
				BLAS_axpy(c, -1*alpha[i], zp1, 1, wp1, 1);
			}
		}

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		if ( k > 0 && k % fuse == 0 ) {
			#pragma omp taskwait
			stop_timer(&elapses[kk]);
			double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
			double sr2norm = norm_r/norm_b;
			residuals[kk] = sr2norm;
			if ( isless(sr2norm, prec) ) {
				fprintf(stderr, "Precision reached\n");
				break;
			}
			kk += 1;
			start_timer();
		}
	}

	#pragma omp taskwait
	if ( k == cgit ) {
		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		stop_timer(&elapses[kk]);
		residuals[kk] = sr2norm;
	}

	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg4_ifcg_v2.log", kk, residuals, elapses);
	free(pool);
	free(alpha);
	free(beta);
	free(gamma);
	free(residuals);
	free(elapses);
	return 0;
}

/* IFCG V2 Centinel */
int CG_IFCG_V2_Centinel(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, int fuse, double orth_fac, int cglog_level)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
	double *pool = calloc(10 * 2 * n, sizeof(double));
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u_i+1 = u_i - alpha_i * q_i
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w_i+1 = w_i - alpha_i * z_i
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta_i * p_i
	offs += 2 * n;
	double *m[2] = {&pool[offs], &pool[offs+n]}; // m_i = M^-1 * w_i
	offs += 2 * n;
	double *n0[2] = {&pool[offs], &pool[offs+n]}; // n_i = A * m_i
	offs += 2 * n;
	double *z[2] = {&pool[offs], &pool[offs+n]}; // z_i = n_i + beta_i * z_i-1
	offs += 2 * n;
	double *q[2] = {&pool[offs], &pool[offs+n]}; // q_i = m_i + beta_i * q_i-1
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha_i * s_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i = w_i + beta_i * s_i-1

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	int bs = (n+bm-1)/bm;
	double *a_gamma = malloc(bs * sizeof(double));
	double *a_delta = malloc(bs * sizeof(double));

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
	bsblas_dcholsolv2(1, bm, n, S, N, r[i], u[i]);
	/* w[0] = A * u[0] */
	hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]);

	start_timer();
	int kk = 0;
	int k;
	for ( k = 0; k < cgit; k++) {

		int iprev = i;
		i = i ^ 0x1;

		/* m[i] = M^-1 * w[i] */
		bsblas_dcholsolv2(1, bm, n, S, N, w[iprev], m[i]);

		/* 
		 * gamma[i] = <r[i], u[i]>
		 */
		bblas_ddot(1, bm, 1, n, 1, r[iprev], u[iprev], &gamma[i]);
		bblas_ddot(1, bm, 1, n, 1, w[iprev], u[iprev], &delta);

		#pragma omp task in(gamma[i]) out(beta[i]) no_copy_deps label(sentinel0)
		{
			if ( k > 0 ) {
				beta[i] = gamma[i]/gamma[iprev];
			} else {
				beta[i] = (double) 0;
			}
		}
#if 0
		bblas_ddot_array(1, bm, 1, n, 1, r[iprev], u[iprev], a_gamma);
		bblas_ddot_array(1, bm, 1, n, 1, w[iprev], u[iprev], a_delta);

		#pragma omp task in(a_gamma[0:bs-1]) out(beta[i]) no_copy_deps label(sentinel0)
		{
			gamma[i] = (double) 0;
			for (int ii = 0; ii < bs; ii++) {
				gamma[i] += a_gamma[ii];
			}
			if ( k > 0 ) {
				beta[i] = gamma[i]/gamma[iprev];
			} else {
				beta[i] = (double) 0;
			}
		}
#endif

		/* 
		 * s_i = w_i + beta_i * s_i-1
		 * p_i = u_i + beta_i * p_i-1
		 */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *sp0 = &(s[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *sp1 = &(s[i])[j];
			double *pp1 = &(p[i])[j];

			#pragma omp task out([c]sp1, [c]pp1) in(beta[i], [c]sp0, [c]pp0, [c]up0, [c]wp0) priority(1) label(alg4_apxy)
			{
				/* s_i = w_i + beta_i * s_i-1 */
				BLAS_cp(c, wp0, 1, sp1, 1);
				BLAS_axpy(c, beta[i], sp0, 1, sp1, 1);

				/* p_i = u_i + beta_i * p_i-1 */
				BLAS_cp(c, up0, 1, pp1, 1);
				BLAS_axpy(c, beta[i], pp0, 1, pp1, 1);
			}
		}

		/* n[i] = A * m[i] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, m[i], fp_nought, n0[i]); 

		#pragma omp task in(beta[i], delta) out(alpha[i]) label(sentinel1) no_copy_deps
		{
			if ( k > 0 ) {
				alpha[i] = gamma[i] / (delta - beta[i] * gamma[i] / alpha[iprev]);
			} else {
				alpha[i] = gamma[i]/delta;
			}
			gamma[iprev] = delta = 0;
		}
#if 0
		#pragma omp task in(beta[i], a_delta[0:bs-1]) out(alpha[i]) label(sentinel1) no_copy_deps
		{
			delta = (double) 0;
			for (int ii= 0; ii < bs; ii++) {
				delta += a_delta[ii];
			}
			if ( k > 0 ) {
				alpha[i] = gamma[i] / (delta - beta[i] * gamma[i] / alpha[iprev]);
			} else {
				alpha[i] = gamma[i]/delta;
			}
		}
#endif

		/* Grand fuse */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *zp0 = &(z[iprev])[j];
			double *qp0 = &(q[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *xp0 = &(x[iprev])[j];
			double *rp0 = &(r[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *zp1 = &(z[i])[j];
			double *qp1 = &(q[i])[j];
			double *sp1 = &(s[i])[j];
			double *pp1 = &(p[i])[j];
			double *xp1 = &(x[i])[j];
			double *rp1 = &(r[i])[j];
			double *up1 = &(u[i])[j];
			double *wp1 = &(w[i])[j];

			double *mp1 = &(m[i])[j];
			double *np1 = &(n0[i])[j];

			#pragma omp task out([c]zp1, [c]qp1, [c]xp1, [c]rp1, [c]up1, [c]wp1) \
				in(alpha[i], [c]zp0, [c]qp0, [c]sp1, [c]pp1, [c]np1, [c]mp1, [c]up0, [c]xp0, [c]wp0, [c]rp0) \
				priority(1) label(alg4_fuse)
			{
				/* z_i = n_i + beta_i * z_i-1 */
				BLAS_cp(c, np1, 1, zp1, 1);
				BLAS_axpy(c, beta[i], zp0, 1, zp1, 1);

				/* q_i = m_i + beta_i * q_i-1 */
				BLAS_cp(c, mp1, 1, qp1, 1);
				BLAS_axpy(c, beta[i], qp0, 1, qp1, 1);

				/* x_i+1 = x_i + alpha_i * p_i */
				BLAS_cp(c, xp0, 1, xp1, 1);
				BLAS_axpy(c, alpha[i], pp1, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				BLAS_cp(c, rp0, 1, rp1, 1);
				BLAS_axpy(c, -1*alpha[i], sp1, 1, rp1, 1);

				/* u_i+1 = u_i - alpha_i * q_i */
				BLAS_cp(c, up0, 1, up1, 1);
				BLAS_axpy(c, -1*alpha[i], qp1, 1, up1, 1);

				/* w_i+1 = w_i - alpha_i * z_i */
				BLAS_cp(c, wp0, 1, wp1, 1);
				BLAS_axpy(c, -1*alpha[i], zp1, 1, wp1, 1);
			}
		}

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		if ( k > 0 && k % fuse == 0 ) {
			#pragma omp taskwait
			stop_timer(&elapses[kk]);
			double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
			double sr2norm = norm_r/norm_b;
			residuals[kk] = sr2norm;
			if ( isless(sr2norm, prec) ) {
				fprintf(stderr, "Precision reached\n");
				break;
			}
			kk += 1;
			start_timer();
		}
	}

	#pragma omp taskwait
	if ( k == cgit ) {
		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		stop_timer(&elapses[kk]);
		residuals[kk] = sr2norm;
	}

	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg4_ifcg_v2_centinel.log", kk, residuals, elapses);
	free(pool);
	free(alpha);
	free(beta);
	free(gamma);
	free(a_gamma);
	free(a_delta);
	free(residuals);
	free(elapses);
	return 0;
}

/* Iteration-fusing Standard PCG */
int CG_IFCG_PCG(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, int fuse, double orth_fac, int cglog_level)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
	double *pool = calloc(4 * 2 * n, sizeof(double));
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha * p_i
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u = M^-1 * r_i+1
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta * p_i
	double  *s = malloc(n * sizeof(double)); // s = Ap

	double *alpha1 = calloc(2, sizeof(double));
	double *alpha2 = calloc(2, sizeof(double));

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
	bsblas_dcholsolv2(1, bm, n, S, N, r[i], u[i]);
	/* p[0] = z[0] */
	bblas_dcopy(1, bm, 1, n, 1, u[i], p[i]);
	/* alpha1[0] = <r[0], z[0]> */
	bblas_ddot(1, bm, 1, n, 1, r[i], u[i], &alpha1[i]);

	int k;
	for ( k = 0; k < cgit; k++ ) {
		start_timer();
		int iprev = i;
		i = i ^ 0x1;

		/* s = A * p[i] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, p[iprev], fp_nought, s); 
		/* alpha2[i] = <s, p[i]> */
		bblas_ddot(1, bm, 1, n, 1, s, p[iprev], &alpha2[i]);
		bblas_dcpaxpy_comb(bm, 1, n, 1, fp_mone, &alpha1[iprev], &alpha2[i], s, p[iprev], r[iprev], x[iprev], r[i], x[i]);

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		/* u[i+1] = M^-1 * r[i+1] */
		bsblas_dcholsolv2(1, bm, n, S, N, r[i], u[i]);
		/* alpha1[i+1] = <r, u> */
		bblas_ddot(1, bm, 1, n, 1, r[i], u[i], &alpha1[i]);
		/* p[i+1] = u[i+1] + transpose(beta[i]) * p[i] */
		bblas_extm_daxpy(1, bm, 1, n, 1, &alpha1[i], &alpha1[iprev], p[iprev], u[i], p[i]); 	

		#pragma omp taskwait

		stop_timer(&elapses[k]);
		alpha1[iprev] = alpha2[iprev] = (double) 0;
//		BLAS_gemm(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, 1, 1, n, FP_ONE, p[i], n, s, n, FP_NOUGHT, &orth, 1);
//		orth = FP_ABS(orth);
//		if (isgreater(orth, porth * orth_fac)){
//			fprintf(stderr, "orth fail %E %E\n", orth, porth*orth_fac);
//			break;
//		}

		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		residuals[k] = sr2norm;
		if ( isless(sr2norm, prec) ) {
			fprintf(stderr, "Precision reached\n");
			break;
		}
//		fprintf(stdout, "%d %E\n", k, residuals[k]);
	}
	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg1.log", k, residuals, elapses);
	free(pool);
	free(alpha1);
	free(alpha2);
	free(residuals);
	free(elapses);
	return 0;
}


/* TODO BPs */

int CG_ALG1_BP(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, double orth_fac, int cglog_level, double *pool)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
//	double *pool = calloc(4 * 2 * n, sizeof(double));
//	#pragma omp register ([4*2*n]pool)
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha * p_i
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u = M^-1 * r_i+1
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta * p_i
	double  *s = malloc(n * sizeof(double)); // s = Ap

	double *alpha1 = calloc(2, sizeof(double));
	double *alpha2 = calloc(2, sizeof(double));

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
	dcholsolv2_nested(1, bm, n, S, N, r[i], u[i]);
	/* p[0] = z[0] */
	bblas_dcopy(1, bm, 1, n, 1, u[i], p[i]);
	/* alpha1[0] = <r[0], z[0]> */
	bblas_ddot(1, bm, 1, n, 1, r[i], u[i], &alpha1[i]);

	int k;
	for ( k = 0; k < cgit; k++ ) {
		start_timer();
		int iprev = i;
		i = i ^ 0x1;

		/* s = A * p[i] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, p[iprev], fp_nought, s); 
		/* alpha2[i] = <s, p[i]> */
		bblas_ddot(1, bm, 1, n, 1, s, p[iprev], &alpha2[i]);
		bblas_dcpaxpy_comb(bm, 1, n, 1, fp_mone, &alpha1[iprev], &alpha2[i], s, p[iprev], r[iprev], x[iprev], r[i], x[i]);

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		/* u[i+1] = M^-1 * r[i+1] */
//		bsblas_dcholsolv2(1, bm, n, S, N, r[i], u[i]);
		dcholsolv2_nested(1, bm, n, S, N, r[i], u[i]);
		/* alpha1[i+1] = <r, u> */
		bblas_ddot(1, bm, 1, n, 1, r[i], u[i], &alpha1[i]);
		/* p[i+1] = u[i+1] + transpose(beta[i]) * p[i] */
		bblas_extm_daxpy(1, bm, 1, n, 1, &alpha1[i], &alpha1[iprev], p[iprev], u[i], p[i]); 	

		#pragma omp taskwait

		stop_timer(&elapses[k]);
		alpha1[iprev] = alpha2[iprev] = (double) 0;
//		BLAS_gemm(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, 1, 1, n, FP_ONE, p[i], n, s, n, FP_NOUGHT, &orth, 1);
//		orth = FP_ABS(orth);
//		if (isgreater(orth, porth * orth_fac)){
//			fprintf(stderr, "orth fail %E %E\n", orth, porth*orth_fac);
//			break;
//		}

		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		residuals[k] = sr2norm;
		if ( isless(sr2norm, prec) ) {
			fprintf(stderr, "Precision reached\n");
			break;
		}
//		fprintf(stdout, "%d %E\n", k, residuals[k]);
	}
	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg1_bp.log", k, residuals, elapses);
//	free(pool);
	free(alpha1);
	free(alpha2);
	free(residuals);
	free(elapses);
	return 0;
}

/* Chronopoulos PCG */
int CG_ALG3_BP(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, double orth_fac, int cglog_level, double *pool)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
//	double *pool = calloc(6 * 2 * n, sizeof(double));
//	#pragma omp register ([6*2*n]pool)	
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha * s_i
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u = M^-1 * r_i+1
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w = A * u_i+1
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta * p_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i+1 = w_i+1 + beta * s_i

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;

	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
	dcholsolv2_nested(1, bm, n, S, N, r[i], u[i]);
	/* w[0] = A * u[0] */
	hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]);
	/* gamma[0] = <r[0], u[0]> */
	bblas_ddot(1, bm, 1, n, 1, r[i], u[i], &gamma[i]);
	/* alpha[0] = gamma[0]/<w[0], u[0]> */
	bblas_ddot(1, bm, 1, n, 1, w[i], u[i], &alpha[i]);

	#pragma omp taskwait on (gamma[i], alpha[i])

	alpha[i] = gamma[i] / alpha[i];

	int k;
	for ( k = 0; k < cgit; k++ ) {
		start_timer();
		int iprev = i;
		i = i ^ 0x1;

		/* Grand fuse */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *pp0 = &(p[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *wp0 = &(w[iprev])[j];
			double *xp0 = &(x[iprev])[j];
			double *rp0 = &(r[iprev])[j];

			double *pp1 = &(p[i])[j];
			double *up1 = &(u[i])[j];
			double *sp1 = &(s[i])[j];
			double *wp1 = &(w[i])[j];
			double *xp1 = &(x[i])[j];
			double *rp1 = &(r[i])[j];

			#pragma omp task out([c]pp1, [c]sp1, [c]xp1, [c]rp1)  priority(1) label(alg3_fuse)
			{
				/* p_i = u_i + beta_i * p_i-1 */
				BLAS_cp(c, up0, 1, pp1, 1);
				BLAS_axpy(c, beta[iprev], pp0, 1, pp1, 1);

				/* s_i = w_i + beta_i * s_i-1 */
				BLAS_cp(c, wp0, 1, sp1, 1);
				BLAS_axpy(c, beta[iprev], sp0, 1, sp1, 1);

				/* x_i+1 = x_i + alpha_i * p_i */
				BLAS_cp(c, xp0, 1, xp1, 1);
				BLAS_axpy(c, alpha[iprev], pp1, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				BLAS_cp(c, rp0, 1, rp1, 1);
				BLAS_axpy(c, -1*alpha[iprev], sp1, 1, rp1, 1);
			}
		}

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		/* u[i+1] = M^-1 * r[i+1] */
		dcholsolv2_nested(1, bm, n, S, N, r[i], u[i]);
		/* w[i+1] = A * u[i+1] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]); 

		cg_ddot2(1, bm, 1, n, 1, r[i], u[i], &gamma[i], w[i], u[i], &delta);

		#pragma omp taskwait

		beta[i] = gamma[i] / gamma[iprev];
		alpha[i] = gamma[i]/(delta - beta[i] * gamma[i] / alpha[iprev]);

		stop_timer(&elapses[k]);
		gamma[iprev] = delta = (double) 0;

		//TODO Implement p-orthogonality check
#if 0
		BLAS_gemm(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, 1, 1, n, FP_ONE, p[i], n, s, n, FP_NOUGHT, &orth, 1);
		orth = FP_ABS(orth);
		if (isgreater(orth, porth * orth_fac)){
			fprintf(stderr, "orth fail %E %E\n", orth, porth*orth_fac);
			break;
		}
#endif

		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		residuals[k] = sr2norm;
		if ( isless(sr2norm, prec) ) {
			fprintf(stderr, "Precision reached\n");
			break;
		}
//		fprintf(stdout, "%d %E\n", k, residuals[k]);
	}

	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg3_bp.log", k, residuals, elapses);
//	free(pool);
	free(alpha);
	free(beta);
	free(gamma);
	free(residuals);
	free(elapses);
	return 0;
}

/* Pipelined PCG */
int CG_ALG4_BP(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, double orth_fac, int cglog_level, double *pool)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
//	double *pool = calloc(10 * 2 * n, sizeof(double));
//	#pragma omp register ([10*2*n]pool)
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u_i+1 = u_i - alpha_i * q_i
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w_i+1 = w_i - alpha_i * z_i
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta_i * p_i
	offs += 2 * n;
	double *m[2] = {&pool[offs], &pool[offs+n]}; // m_i = M^-1 * w_i
	offs += 2 * n;
	double *n0[2] = {&pool[offs], &pool[offs+n]}; // n_i = A * m_i
	offs += 2 * n;
	double *z[2] = {&pool[offs], &pool[offs+n]}; // z_i = n_i + beta_i * z_i-1
	offs += 2 * n;
	double *q[2] = {&pool[offs], &pool[offs+n]}; // q_i = m_i + beta_i * q_i-1
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha_i * s_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i = w_i + beta_i * s_i-1

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
	dcholsolv2_nested(1, bm, n, S, N, r[i], u[i]);
	/* w[0] = A * u[0] */
	hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]);

	int k;
	for ( k = 0; k < cgit; k++ ) {
		start_timer();
		int iprev = i;
		i = i ^ 0x1;

		/* m[i] = M^-1 * w[i] */
		dcholsolv2_nested(1, bm, n, S, N, w[iprev], m[i]);
		/* n[i] = A * m[i] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, m[i], fp_nought, n0[i]); 
		/* 
		 * gamma[i] = <r[i], u[i]>
		 * delta = <w[i], u[i]>
		 */
		cg_ddot2(1, bm, 1, n, 1, r[iprev], u[iprev], &gamma[i], w[iprev], u[iprev], &delta);

		#pragma omp taskwait on(gamma[i], delta)
		if ( k > 0 ) {
			beta[i] = gamma[i]/gamma[iprev];
			alpha[i] = gamma[i] / (delta - beta[i] * gamma[i] / alpha[iprev]);
		} else {
			beta[i] = (double) 0;
			alpha[i] = gamma[i]/delta;
		}
		gamma[iprev] = delta = 0;

		/* Grand fuse */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *zp0 = &(z[iprev])[j];
			double *qp0 = &(q[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *xp0 = &(x[iprev])[j];
			double *rp0 = &(r[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *zp1 = &(z[i])[j];
			double *qp1 = &(q[i])[j];
			double *sp1 = &(s[i])[j];
			double *pp1 = &(p[i])[j];
			double *xp1 = &(x[i])[j];
			double *rp1 = &(r[i])[j];
			double *up1 = &(u[i])[j];
			double *wp1 = &(w[i])[j];

			double *mp1 = &(m[i])[j];
			double *np1 = &(n0[i])[j];

			#pragma omp task out([c]zp1, [c]qp1, [c]sp1, [c]pp1, [c]xp1, [c]rp1, [c]up1, [c]wp1) \
				in([c]zp0, [c]qp0, [c]sp0, [c]pp0, [c]np1, [c]mp1, [c]up0, [c]xp0, [c]wp0, [c]rp0) \
				priority(1) label(alg4_fuse)
			{
				/* z_i = n_i + beta_i * z_i-1 */
				BLAS_cp(c, np1, 1, zp1, 1);
				BLAS_axpy(c, beta[i], zp0, 1, zp1, 1);

				/* q_i = m_i + beta_i * q_i-1 */
				BLAS_cp(c, mp1, 1, qp1, 1);
				BLAS_axpy(c, beta[i], qp0, 1, qp1, 1);

				/* s_i = w_i + beta_i * s_i-1 */
				BLAS_cp(c, wp0, 1, sp1, 1);
				BLAS_axpy(c, beta[i], sp0, 1, sp1, 1);

				/* p_i = u_i + beta_i * p_i-1 */
				BLAS_cp(c, up0, 1, pp1, 1);
				BLAS_axpy(c, beta[i], pp0, 1, pp1, 1);

				/* x_i+1 = x_i + alpha_i * p_i */
				BLAS_cp(c, xp0, 1, xp1, 1);
				BLAS_axpy(c, alpha[i], pp1, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				BLAS_cp(c, rp0, 1, rp1, 1);
				BLAS_axpy(c, -1*alpha[i], sp1, 1, rp1, 1);

				/* u_i+1 = u_i - alpha_i * q_i */
				BLAS_cp(c, up0, 1, up1, 1);
				BLAS_axpy(c, -1*alpha[i], qp1, 1, up1, 1);

				/* w_i+1 = w_i - alpha_i * z_i */
				BLAS_cp(c, wp0, 1, wp1, 1);
				BLAS_axpy(c, -1*alpha[i], zp1, 1, wp1, 1);
			}
		}

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		#pragma omp taskwait

//		gamma[iprev] = delta = 0;
		stop_timer(&elapses[k]);

		//TODO Implement p-orthogonality check
#if 0
		BLAS_gemm(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, 1, 1, n, FP_ONE, p[i], n, s, n, FP_NOUGHT, &orth, 1);
		orth = FP_ABS(orth);
		if (isgreater(orth, porth * orth_fac)){
			fprintf(stderr, "orth fail %E %E\n", orth, porth*orth_fac);
			break;
		}
#endif

		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		residuals[k] = sr2norm;
		if ( isless(sr2norm, prec) ) {
			fprintf(stderr, "Precision reached\n");
			break;
		}
//		fprintf(stdout, "%d %E\n", k, residuals[k]);
	}

	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg4_bp.log", k, residuals, elapses);
//	free(pool);
	free(alpha);
	free(beta);
	free(gamma);
	free(residuals);
	free(elapses);
	return 0;
}

/* Gropp PCG */
int CG_ALG7_BP(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, double orth_fac, int cglog_level, double *pool)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
//	double *pool = calloc(7 * 2 * n, sizeof(double));
//	#pragma omp register ([7*2*n]pool)
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u_i+1 = u_i - alpha_i * q_i
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w_i+1 = w_i - alpha_i * z_i
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta_i * p_i
	offs += 2 * n;
	double *q[2] = {&pool[offs], &pool[offs+n]}; // q_i = m_i + beta_i * q_i-1
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha_i * s_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i = w_i + beta_i * s_i-1

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
	dcholsolv2_nested(1, bm, n, S, N, r[i], u[i]);
	/* p[0] = u[0] */
	bblas_dcopy(1, bm, 1, n, 1, u[i], p[i]);
	/* s[0] = A * p[0] */
	hbsblas_dcsrmv(1, bm, fp_one, Ahbh, p[i], fp_nought, s[i]);
	/* gamma[0] = <r[0], u[0]> */
	bblas_ddot(1, bm, 1, n, 1, r[i], u[i], &gamma[i]);

	int k;
	for ( k = 0; k < cgit; k++ ) {
		start_timer();
		int iprev = i;
		i = i ^ 0x1;

		/* delta = <p[i], s[i]> */
		bblas_ddot(1, bm, 1, n, 1, p[iprev], s[iprev], &delta);
		/* q[i] = M^-1 * s[i] */
		dcholsolv2_nested(1, bm, n, S, N, s[iprev], q[i]);

		#pragma omp taskwait on(delta)
		alpha[i] = gamma[iprev]/delta;

		/* Axpy fuse x,r,u */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *xp0 = &(x[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *rp0 = &(r[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *qp0 = &(q[iprev])[j];

			double *xp1 = &(x[i])[j];
			double *pp1 = &(p[i])[j];
			double *rp1 = &(r[i])[j];
			double *sp1 = &(s[i])[j];
			double *up1 = &(u[i])[j];
			double *qp1 = &(q[i])[j];

			#pragma omp task in([c]qp1) out([c]xp1, [c]rp1, [c]up1) priority(1) label(alg7_fuse0)
			{
				/* x_i+1 = x_i + alpha_i * p_i */
				BLAS_cp(c, xp0, 1, xp1, 1);
				BLAS_axpy(c, alpha[i], pp0, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				BLAS_cp(c, rp0, 1, rp1, 1);
				BLAS_axpy(c, -1*alpha[i], sp0, 1, rp1, 1);

				/* u_i+1 = u_i - alpha_i * q_i */
				BLAS_cp(c, up0, 1, up1, 1);
				BLAS_axpy(c, -1*alpha[i], qp1, 1, up1, 1);
			}
		}

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait			
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		/* gamma[i+1] = <r[i+1], u[i+1]> */
		bblas_ddot(1, bm, 1, n, 1, r[i], u[i], &gamma[i]);
		/* w[i+1] = A * u[i+1] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]);

		#pragma omp taskwait on(gamma[i])
		beta[i] = gamma[i]/gamma[iprev];

		/* Axpy fuse p,s */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *pp0 = &(p[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *pp1 = &(p[i])[j];
			double *sp1 = &(s[i])[j];
			double *up1 = &(u[i])[j];
			double *wp1 = &(w[i])[j];

			#pragma omp task in([c]up1, [c]pp0, [c]wp1, [c]sp0) priority(1) label(alg7_fuse1)
			{
				/* p_i = u_i + beta_i * p_i-1 */
				BLAS_cp(c, up1, 1, pp1, 1);
				BLAS_axpy(c, beta[i], pp0, 1, pp1, 1);

				/* s_i = w_i + beta_i * s_i-1 */
				BLAS_cp(c, wp1, 1, sp1, 1);
				BLAS_axpy(c, beta[i], sp0, 1, sp1, 1);
			}
		}

		#pragma omp taskwait

		gamma[iprev] = delta = 0;
		stop_timer(&elapses[k]);

		//TODO Implement p-orthogonality check
#if 0
		BLAS_gemm(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, 1, 1, n, FP_ONE, p[i], n, s, n, FP_NOUGHT, &orth, 1);
		orth = FP_ABS(orth);
		if (isgreater(orth, porth * orth_fac)){
			fprintf(stderr, "orth fail %E %E\n", orth, porth*orth_fac);
			break;
		}
#endif

		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		residuals[k] = sr2norm;
		if ( isless(sr2norm, prec) ) {
			fprintf(stderr, "Precision reached\n");
			break;
		}
//		fprintf(stdout, "%d %E\n", k, residuals[k]);
	}

	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg7_bp.log", k, residuals, elapses);
//	free(pool);
	free(alpha);
	free(beta);
	free(gamma);
	free(residuals);
	free(elapses);
	return 0;
}

/* IFCG BP. */
int CG_IFCG_BP(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, int fuse, double orth_fac, int cglog_level, double *pool)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
//	double *pool = calloc(10 * 2 * n, sizeof(double));
//	#pragma omp register ([10*2*n]pool)
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u_i+1 = u_i - alpha_i * q_i
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w_i+1 = w_i - alpha_i * z_i
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta_i * p_i
	offs += 2 * n;
	double *m[2] = {&pool[offs], &pool[offs+n]}; // m_i = M^-1 * w_i
	offs += 2 * n;
	double *n0[2] = {&pool[offs], &pool[offs+n]}; // n_i = A * m_i
	offs += 2 * n;
	double *z[2] = {&pool[offs], &pool[offs+n]}; // z_i = n_i + beta_i * z_i-1
	offs += 2 * n;
	double *q[2] = {&pool[offs], &pool[offs+n]}; // q_i = m_i + beta_i * q_i-1
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha_i * s_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i = w_i + beta_i * s_i-1

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
	dcholsolv2_nested(1, bm, n, S, N, r[i], u[i]);
	/* w[0] = A * u[0] */
	hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]);

	start_timer();
	int kk = 0;
	int k;
	for ( k = 0; k < cgit; k++) {

		int iprev = i;
		i = i ^ 0x1;

		/* m[i] = M^-1 * w[i] */
//		dcholsolv2_blk(1, n, SA, NA, w[iprev], m[i]);
//		#pragma omp task in([n](w[iprev])) out([n](m[i])) label(cholsolv2_seq)
//		bsblas_dcholsolv2_seq(1, bm, n, S, N, w[iprev], m[i]);
		dcholsolv2_nested(1, bm, n, S, N, w[iprev], m[i]);
		/* n[i] = A * m[i] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, m[i], fp_nought, n0[i]); 
		/* 
		 * gamma[i] = <r[i], u[i]>
		 * delta = <w[i], u[i]>
		 */
		cg_ddot2(1, bm, 1, n, 1, r[iprev], u[iprev], &gamma[i], w[iprev], u[iprev], &delta);

		#pragma omp taskwait on(gamma[i], delta)
		if ( k > 0 ) {
			beta[i] = gamma[i]/gamma[iprev];
			alpha[i] = gamma[i] / (delta - beta[i] * gamma[i] / alpha[iprev]);
		} else {
			beta[i] = (double) 0;
			alpha[i] = gamma[i]/delta;
		}
		gamma[iprev] = delta = 0;

		/* Grand fuse */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *zp0 = &(z[iprev])[j];
			double *qp0 = &(q[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *xp0 = &(x[iprev])[j];
			double *rp0 = &(r[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *zp1 = &(z[i])[j];
			double *qp1 = &(q[i])[j];
			double *sp1 = &(s[i])[j];
			double *pp1 = &(p[i])[j];
			double *xp1 = &(x[i])[j];
			double *rp1 = &(r[i])[j];
			double *up1 = &(u[i])[j];
			double *wp1 = &(w[i])[j];

			double *mp1 = &(m[i])[j];
			double *np1 = &(n0[i])[j];

			#pragma omp task out([c]zp1, [c]qp1, [c]sp1, [c]pp1, [c]xp1, [c]rp1, [c]up1, [c]wp1) \
				in([c]zp0, [c]qp0, [c]sp0, [c]pp0, [c]np1, [c]mp1, [c]up0, [c]xp0, [c]wp0, [c]rp0) \
				priority(1) label(alg4_fuse)
			{
				/* z_i = n_i + beta_i * z_i-1 */
				BLAS_cp(c, np1, 1, zp1, 1);
				BLAS_axpy(c, beta[i], zp0, 1, zp1, 1);

				/* q_i = m_i + beta_i * q_i-1 */
				BLAS_cp(c, mp1, 1, qp1, 1);
				BLAS_axpy(c, beta[i], qp0, 1, qp1, 1);

				/* s_i = w_i + beta_i * s_i-1 */
				BLAS_cp(c, wp0, 1, sp1, 1);
				BLAS_axpy(c, beta[i], sp0, 1, sp1, 1);

				/* p_i = u_i + beta_i * p_i-1 */
				BLAS_cp(c, up0, 1, pp1, 1);
				BLAS_axpy(c, beta[i], pp0, 1, pp1, 1);

				/* x_i+1 = x_i + alpha_i * p_i */
				BLAS_cp(c, xp0, 1, xp1, 1);
				BLAS_axpy(c, alpha[i], pp1, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				BLAS_cp(c, rp0, 1, rp1, 1);
				BLAS_axpy(c, -1*alpha[i], sp1, 1, rp1, 1);

				/* u_i+1 = u_i - alpha_i * q_i */
				BLAS_cp(c, up0, 1, up1, 1);
				BLAS_axpy(c, -1*alpha[i], qp1, 1, up1, 1);

				/* w_i+1 = w_i - alpha_i * z_i */
				BLAS_cp(c, wp0, 1, wp1, 1);
				BLAS_axpy(c, -1*alpha[i], zp1, 1, wp1, 1);
			}
		}

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		if ( k > 0 && k % fuse == 0 ) {
			#pragma omp taskwait
			stop_timer(&elapses[kk]);
			double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
			double sr2norm = norm_r/norm_b;
			residuals[kk] = sr2norm;
			if ( isless(sr2norm, prec) ) {
				fprintf(stderr, "Precision reached\n");
				break;
			}
			kk += 1;
			start_timer();
		}
	}

	#pragma omp taskwait
	if ( k == cgit ) {
		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		stop_timer(&elapses[kk]);
		residuals[kk] = sr2norm;
	}

	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg4_ifcg_bp.log", kk, residuals, elapses);
//	free(pool);
	free(alpha);
	free(beta);
	free(gamma);
	free(residuals);
	free(elapses);
	return 0;
}

/* IFCG V2 BP*/
int CG_IFCG_V2_BP(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, int fuse, double orth_fac, int cglog_level, double *pool)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
//	double *pool = calloc(10 * 2 * n, sizeof(double));
//	#pragma omp register ([10*2*n]pool)
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u_i+1 = u_i - alpha_i * q_i
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w_i+1 = w_i - alpha_i * z_i
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta_i * p_i
	offs += 2 * n;
	double *m[2] = {&pool[offs], &pool[offs+n]}; // m_i = M^-1 * w_i
	offs += 2 * n;
	double *n0[2] = {&pool[offs], &pool[offs+n]}; // n_i = A * m_i
	offs += 2 * n;
	double *z[2] = {&pool[offs], &pool[offs+n]}; // z_i = n_i + beta_i * z_i-1
	offs += 2 * n;
	double *q[2] = {&pool[offs], &pool[offs+n]}; // q_i = m_i + beta_i * q_i-1
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha_i * s_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i = w_i + beta_i * s_i-1

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
//	dcholsolv2_blk(1, n, SA, NA, r[i], u[i]);
	dcholsolv2_nested(1, bm, n, S, N, r[i], u[i]);
//	#pragma omp task in([n](r[i])) out([n](u[i])) label(cholsolv2_seq)
//	bsblas_dcholsolv2_seq(1, bm, n, S, N, r[i], u[i]);
//	bsblas_dcholsolv2(1, bm, n, S, N, r[i], u[i]);
	/* w[0] = A * u[0] */
	hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]);

	start_timer();
	int kk = 0;
	int k;
	for ( k = 0; k < cgit; k++) {

		int iprev = i;
		i = i ^ 0x1;

		/* m[i] = M^-1 * w[i] */
		dcholsolv2_nested(1, bm, n, S, N, w[iprev], m[i]);
//		dcholsolv2_blk(1, n, SA, NA, w[iprev], m[i]);
//		bsblas_dcholsolv2(1, bm, n, S, N, w[iprev], m[i]);
//		#pragma omp task in([n](w[iprev])) out([n](m[i])) label(cholsolv2_seq)
//		bsblas_dcholsolv2_seq(1, bm, n, S, N, w[iprev], m[i]);

		/* 
		 * gamma[i] = <r[i], u[i]>
		 */
		bblas_ddot(1, bm, 1, n, 1, r[iprev], u[iprev], &gamma[i]);
		bblas_ddot(1, bm, 1, n, 1, w[iprev], u[iprev], &delta);

		#pragma omp taskwait on(gamma[i])
		if ( k > 0 ) {
			beta[i] = gamma[i]/gamma[iprev];
//			alpha[i] = gamma[i] / (delta - beta[i] * gamma[i] / alpha[iprev]);
		} else {
			beta[i] = (double) 0;
//			alpha[i] = gamma[i]/delta;
		}
//		gamma[iprev] = delta = 0;

		/* 
		 * s_i = w_i + beta_i * s_i-1
		 * p_i = u_i + beta_i * p_i-1
		 */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *sp0 = &(s[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *sp1 = &(s[i])[j];
			double *pp1 = &(p[i])[j];

			#pragma omp task out([c]sp1, [c]pp1) in([c]sp0, [c]pp0, [c]up0, [c]wp0) priority(1) label(alg4_apxy)
			{
				/* s_i = w_i + beta_i * s_i-1 */
				BLAS_cp(c, wp0, 1, sp1, 1);
				BLAS_axpy(c, beta[i], sp0, 1, sp1, 1);

				/* p_i = u_i + beta_i * p_i-1 */
				BLAS_cp(c, up0, 1, pp1, 1);
				BLAS_axpy(c, beta[i], pp0, 1, pp1, 1);
			}
		}

		/* n[i] = A * m[i] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, m[i], fp_nought, n0[i]); 

//		bblas_ddot(1, bm, 1, n, 1, w[iprev], u[iprev], &delta);
		#pragma omp taskwait on(delta)
		if ( k > 0 ) {
			alpha[i] = gamma[i] / (delta - beta[i] * gamma[i] / alpha[iprev]);
		} else {
			alpha[i] = gamma[i]/delta;
		}
		gamma[iprev] = delta = 0;

		/* Grand fuse */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *zp0 = &(z[iprev])[j];
			double *qp0 = &(q[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *xp0 = &(x[iprev])[j];
			double *rp0 = &(r[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *zp1 = &(z[i])[j];
			double *qp1 = &(q[i])[j];
			double *sp1 = &(s[i])[j];
			double *pp1 = &(p[i])[j];
			double *xp1 = &(x[i])[j];
			double *rp1 = &(r[i])[j];
			double *up1 = &(u[i])[j];
			double *wp1 = &(w[i])[j];

			double *mp1 = &(m[i])[j];
			double *np1 = &(n0[i])[j];

			#pragma omp task out([c]zp1, [c]qp1, [c]xp1, [c]rp1, [c]up1, [c]wp1) \
				in([c]zp0, [c]qp0, [c]sp0, [c]pp0, [c]np1, [c]mp1, [c]up0, [c]xp0, [c]wp0, [c]rp0) \
				priority(1) label(alg4_fuse)
			{
				/* z_i = n_i + beta_i * z_i-1 */
				BLAS_cp(c, np1, 1, zp1, 1);
				BLAS_axpy(c, beta[i], zp0, 1, zp1, 1);

				/* q_i = m_i + beta_i * q_i-1 */
				BLAS_cp(c, mp1, 1, qp1, 1);
				BLAS_axpy(c, beta[i], qp0, 1, qp1, 1);

				/* x_i+1 = x_i + alpha_i * p_i */
				BLAS_cp(c, xp0, 1, xp1, 1);
				BLAS_axpy(c, alpha[i], pp1, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				BLAS_cp(c, rp0, 1, rp1, 1);
				BLAS_axpy(c, -1*alpha[i], sp1, 1, rp1, 1);

				/* u_i+1 = u_i - alpha_i * q_i */
				BLAS_cp(c, up0, 1, up1, 1);
				BLAS_axpy(c, -1*alpha[i], qp1, 1, up1, 1);

				/* w_i+1 = w_i - alpha_i * z_i */
				BLAS_cp(c, wp0, 1, wp1, 1);
				BLAS_axpy(c, -1*alpha[i], zp1, 1, wp1, 1);
			}
		}

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		if ( k > 0 && k % fuse == 0 ) {
			#pragma omp taskwait
			stop_timer(&elapses[kk]);
			double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
			double sr2norm = norm_r/norm_b;
			residuals[kk] = sr2norm;
			if ( isless(sr2norm, prec) ) {
				fprintf(stderr, "Precision reached\n");
				break;
			}
			kk += 1;
			start_timer();
		}
	}

	#pragma omp taskwait
	if ( k == cgit ) {
		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		stop_timer(&elapses[kk]);
		residuals[kk] = sr2norm;
	}

	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg4_ifcg_v2_bp.log", kk, residuals, elapses);
//	free(pool);
	free(alpha);
	free(beta);
	free(gamma);
	free(residuals);
	free(elapses);
	return 0;
}


/* TODO ILU */

int CG_ALG1_ILU_BP(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, double orth_fac, int cglog_level, double *pool)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha * p_i
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u = M^-1 * r_i+1
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta * p_i
	double  *s = malloc(n * sizeof(double)); // s = Ap

	double *alpha1 = calloc(2, sizeof(double));
	double *alpha2 = calloc(2, sizeof(double));

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
//	dcholsolv2_nested(1, bm, n, S, N, r[i], u[i]);
	dcholsolv2_blk(1, n, SA, NA, r[i], u[i]);
	/* p[0] = z[0] */
	bblas_dcopy(1, bm, 1, n, 1, u[i], p[i]);
	/* alpha1[0] = <r[0], z[0]> */
	bblas_ddot(1, bm, 1, n, 1, r[i], u[i], &alpha1[i]);

	int k;
	for ( k = 0; k < cgit; k++ ) {
		start_timer();
		int iprev = i;
		i = i ^ 0x1;

		/* s = A * p[i] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, p[iprev], fp_nought, s); 
		/* alpha2[i] = <s, p[i]> */
		bblas_ddot(1, bm, 1, n, 1, s, p[iprev], &alpha2[i]);
		bblas_dcpaxpy_comb(bm, 1, n, 1, fp_mone, &alpha1[iprev], &alpha2[i], s, p[iprev], r[iprev], x[iprev], r[i], x[i]);

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		/* u[i+1] = M^-1 * r[i+1] */
//		bsblas_dcholsolv2(1, bm, n, S, N, r[i], u[i]);
//		dcholsolv2_nested(1, bm, n, S, N, r[i], u[i]);
		dcholsolv2_blk(1, n, SA, NA, r[i], u[i]);
		/* alpha1[i+1] = <r, u> */
		bblas_ddot(1, bm, 1, n, 1, r[i], u[i], &alpha1[i]);
		/* p[i+1] = u[i+1] + transpose(beta[i]) * p[i] */
		bblas_extm_daxpy(1, bm, 1, n, 1, &alpha1[i], &alpha1[iprev], p[iprev], u[i], p[i]); 	

		#pragma omp taskwait

		stop_timer(&elapses[k]);
		alpha1[iprev] = alpha2[iprev] = (double) 0;
//		BLAS_gemm(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, 1, 1, n, FP_ONE, p[i], n, s, n, FP_NOUGHT, &orth, 1);
//		orth = FP_ABS(orth);
//		if (isgreater(orth, porth * orth_fac)){
//			fprintf(stderr, "orth fail %E %E\n", orth, porth*orth_fac);
//			break;
//		}

		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		residuals[k] = sr2norm;
		if ( isless(sr2norm, prec) ) {
			fprintf(stderr, "Precision reached\n");
			break;
		}
//		fprintf(stdout, "%d %E\n", k, residuals[k]);
	}
	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg1_ilu_bp.log", k, residuals, elapses);
	free(alpha1);
	free(alpha2);
	free(residuals);
	free(elapses);
	return 0;
}

/* Chronopoulos PCG */
int CG_ALG3_ILU_BP(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, double orth_fac, int cglog_level, double *pool)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha * s_i
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u = M^-1 * r_i+1
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w = A * u_i+1
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta * p_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i+1 = w_i+1 + beta * s_i

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;

	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
//	dcholsolv2_nested(1, bm, n, S, N, r[i], u[i]);
	dcholsolv2_blk(1, n, SA, NA, r[i], u[i]);
	/* w[0] = A * u[0] */
	hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]);
	/* gamma[0] = <r[0], u[0]> */
	bblas_ddot(1, bm, 1, n, 1, r[i], u[i], &gamma[i]);
	/* alpha[0] = gamma[0]/<w[0], u[0]> */
	bblas_ddot(1, bm, 1, n, 1, w[i], u[i], &alpha[i]);

	#pragma omp taskwait on (gamma[i], alpha[i])

	alpha[i] = gamma[i] / alpha[i];

	int k;
	for ( k = 0; k < cgit; k++ ) {
		start_timer();
		int iprev = i;
		i = i ^ 0x1;

		/* Grand fuse */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *pp0 = &(p[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *wp0 = &(w[iprev])[j];
			double *xp0 = &(x[iprev])[j];
			double *rp0 = &(r[iprev])[j];

			double *pp1 = &(p[i])[j];
			double *up1 = &(u[i])[j];
			double *sp1 = &(s[i])[j];
			double *wp1 = &(w[i])[j];
			double *xp1 = &(x[i])[j];
			double *rp1 = &(r[i])[j];

			#pragma omp task out([c]pp1, [c]sp1, [c]xp1, [c]rp1)  priority(1) label(alg3_fuse)
			{
				/* p_i = u_i + beta_i * p_i-1 */
				BLAS_cp(c, up0, 1, pp1, 1);
				BLAS_axpy(c, beta[iprev], pp0, 1, pp1, 1);

				/* s_i = w_i + beta_i * s_i-1 */
				BLAS_cp(c, wp0, 1, sp1, 1);
				BLAS_axpy(c, beta[iprev], sp0, 1, sp1, 1);

				/* x_i+1 = x_i + alpha_i * p_i */
				BLAS_cp(c, xp0, 1, xp1, 1);
				BLAS_axpy(c, alpha[iprev], pp1, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				BLAS_cp(c, rp0, 1, rp1, 1);
				BLAS_axpy(c, -1*alpha[iprev], sp1, 1, rp1, 1);
			}
		}

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		/* u[i+1] = M^-1 * r[i+1] */
//		dcholsolv2_nested(1, bm, n, S, N, r[i], u[i]);
		dcholsolv2_blk(1, n, SA, NA, r[i], u[i]);
		/* w[i+1] = A * u[i+1] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]); 

		cg_ddot2(1, bm, 1, n, 1, r[i], u[i], &gamma[i], w[i], u[i], &delta);

		#pragma omp taskwait

		beta[i] = gamma[i] / gamma[iprev];
		alpha[i] = gamma[i]/(delta - beta[i] * gamma[i] / alpha[iprev]);

		stop_timer(&elapses[k]);
		gamma[iprev] = delta = (double) 0;

		//TODO Implement p-orthogonality check
#if 0
		BLAS_gemm(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, 1, 1, n, FP_ONE, p[i], n, s, n, FP_NOUGHT, &orth, 1);
		orth = FP_ABS(orth);
		if (isgreater(orth, porth * orth_fac)){
			fprintf(stderr, "orth fail %E %E\n", orth, porth*orth_fac);
			break;
		}
#endif

		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		residuals[k] = sr2norm;
		if ( isless(sr2norm, prec) ) {
			fprintf(stderr, "Precision reached\n");
			break;
		}
//		fprintf(stdout, "%d %E\n", k, residuals[k]);
	}

	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg3_ilu_bp.log", k, residuals, elapses);
	free(alpha);
	free(beta);
	free(gamma);
	free(residuals);
	free(elapses);
	return 0;
}

/* Pipelined PCG */
int CG_ALG4_ILU_BP(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, double orth_fac, int cglog_level, double *pool)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u_i+1 = u_i - alpha_i * q_i
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w_i+1 = w_i - alpha_i * z_i
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta_i * p_i
	offs += 2 * n;
	double *m[2] = {&pool[offs], &pool[offs+n]}; // m_i = M^-1 * w_i
	offs += 2 * n;
	double *n0[2] = {&pool[offs], &pool[offs+n]}; // n_i = A * m_i
	offs += 2 * n;
	double *z[2] = {&pool[offs], &pool[offs+n]}; // z_i = n_i + beta_i * z_i-1
	offs += 2 * n;
	double *q[2] = {&pool[offs], &pool[offs+n]}; // q_i = m_i + beta_i * q_i-1
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha_i * s_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i = w_i + beta_i * s_i-1

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
//	dcholsolv2_nested(1, bm, n, S, N, r[i], u[i]);
	dcholsolv2_blk(1, n, SA, NA, r[i], u[i]);
	/* w[0] = A * u[0] */
	hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]);

	int k;
	for ( k = 0; k < cgit; k++ ) {
		start_timer();
		int iprev = i;
		i = i ^ 0x1;

		/* m[i] = M^-1 * w[i] */
//		dcholsolv2_nested(1, bm, n, S, N, w[iprev], m[i]);
		dcholsolv2_blk(1, n, SA, NA, w[iprev], m[i]);
		/* n[i] = A * m[i] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, m[i], fp_nought, n0[i]); 
		/* 
		 * gamma[i] = <r[i], u[i]>
		 * delta = <w[i], u[i]>
		 */
		cg_ddot2(1, bm, 1, n, 1, r[iprev], u[iprev], &gamma[i], w[iprev], u[iprev], &delta);

		#pragma omp taskwait on(gamma[i], delta)
		if ( k > 0 ) {
			beta[i] = gamma[i]/gamma[iprev];
			alpha[i] = gamma[i] / (delta - beta[i] * gamma[i] / alpha[iprev]);
		} else {
			beta[i] = (double) 0;
			alpha[i] = gamma[i]/delta;
		}
		gamma[iprev] = delta = 0;

		/* Grand fuse */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *zp0 = &(z[iprev])[j];
			double *qp0 = &(q[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *xp0 = &(x[iprev])[j];
			double *rp0 = &(r[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *zp1 = &(z[i])[j];
			double *qp1 = &(q[i])[j];
			double *sp1 = &(s[i])[j];
			double *pp1 = &(p[i])[j];
			double *xp1 = &(x[i])[j];
			double *rp1 = &(r[i])[j];
			double *up1 = &(u[i])[j];
			double *wp1 = &(w[i])[j];

			double *mp1 = &(m[i])[j];
			double *np1 = &(n0[i])[j];

			#pragma omp task out([c]zp1, [c]qp1, [c]sp1, [c]pp1, [c]xp1, [c]rp1, [c]up1, [c]wp1) \
				in([c]zp0, [c]qp0, [c]sp0, [c]pp0, [c]np1, [c]mp1, [c]up0, [c]xp0, [c]wp0, [c]rp0) \
				priority(1) label(alg4_fuse)
			{
				/* z_i = n_i + beta_i * z_i-1 */
				BLAS_cp(c, np1, 1, zp1, 1);
				BLAS_axpy(c, beta[i], zp0, 1, zp1, 1);

				/* q_i = m_i + beta_i * q_i-1 */
				BLAS_cp(c, mp1, 1, qp1, 1);
				BLAS_axpy(c, beta[i], qp0, 1, qp1, 1);

				/* s_i = w_i + beta_i * s_i-1 */
				BLAS_cp(c, wp0, 1, sp1, 1);
				BLAS_axpy(c, beta[i], sp0, 1, sp1, 1);

				/* p_i = u_i + beta_i * p_i-1 */
				BLAS_cp(c, up0, 1, pp1, 1);
				BLAS_axpy(c, beta[i], pp0, 1, pp1, 1);

				/* x_i+1 = x_i + alpha_i * p_i */
				BLAS_cp(c, xp0, 1, xp1, 1);
				BLAS_axpy(c, alpha[i], pp1, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				BLAS_cp(c, rp0, 1, rp1, 1);
				BLAS_axpy(c, -1*alpha[i], sp1, 1, rp1, 1);

				/* u_i+1 = u_i - alpha_i * q_i */
				BLAS_cp(c, up0, 1, up1, 1);
				BLAS_axpy(c, -1*alpha[i], qp1, 1, up1, 1);

				/* w_i+1 = w_i - alpha_i * z_i */
				BLAS_cp(c, wp0, 1, wp1, 1);
				BLAS_axpy(c, -1*alpha[i], zp1, 1, wp1, 1);
			}
		}

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		#pragma omp taskwait

//		gamma[iprev] = delta = 0;
		stop_timer(&elapses[k]);

		//TODO Implement p-orthogonality check
#if 0
		BLAS_gemm(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, 1, 1, n, FP_ONE, p[i], n, s, n, FP_NOUGHT, &orth, 1);
		orth = FP_ABS(orth);
		if (isgreater(orth, porth * orth_fac)){
			fprintf(stderr, "orth fail %E %E\n", orth, porth*orth_fac);
			break;
		}
#endif

		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		residuals[k] = sr2norm;
		if ( isless(sr2norm, prec) ) {
			fprintf(stderr, "Precision reached\n");
			break;
		}
//		fprintf(stdout, "%d %E\n", k, residuals[k]);
	}

	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg4_ilu_bp.log", k, residuals, elapses);
	free(alpha);
	free(beta);
	free(gamma);
	free(residuals);
	free(elapses);
	return 0;
}

/* Gropp PCG */
int CG_ALG7_ILU_BP(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, double orth_fac, int cglog_level, double *pool)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u_i+1 = u_i - alpha_i * q_i
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w_i+1 = w_i - alpha_i * z_i
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta_i * p_i
	offs += 2 * n;
	double *q[2] = {&pool[offs], &pool[offs+n]}; // q_i = m_i + beta_i * q_i-1
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha_i * s_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i = w_i + beta_i * s_i-1

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
//	dcholsolv2_nested(1, bm, n, S, N, r[i], u[i]);
	dcholsolv2_blk(1, n, SA, NA, r[i], u[i]);
	/* p[0] = u[0] */
	bblas_dcopy(1, bm, 1, n, 1, u[i], p[i]);
	/* s[0] = A * p[0] */
	hbsblas_dcsrmv(1, bm, fp_one, Ahbh, p[i], fp_nought, s[i]);
	/* gamma[0] = <r[0], u[0]> */
	bblas_ddot(1, bm, 1, n, 1, r[i], u[i], &gamma[i]);

	int k;
	for ( k = 0; k < cgit; k++ ) {
		start_timer();
		int iprev = i;
		i = i ^ 0x1;

		/* delta = <p[i], s[i]> */
		bblas_ddot(1, bm, 1, n, 1, p[iprev], s[iprev], &delta);
		/* q[i] = M^-1 * s[i] */
//		dcholsolv2_nested(1, bm, n, S, N, s[iprev], q[i]);
		dcholsolv2_blk(1, n, SA, NA, s[iprev], q[i]);

		#pragma omp taskwait on(delta)
		alpha[i] = gamma[iprev]/delta;

		/* Axpy fuse x,r,u */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *xp0 = &(x[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *rp0 = &(r[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *qp0 = &(q[iprev])[j];

			double *xp1 = &(x[i])[j];
			double *pp1 = &(p[i])[j];
			double *rp1 = &(r[i])[j];
			double *sp1 = &(s[i])[j];
			double *up1 = &(u[i])[j];
			double *qp1 = &(q[i])[j];

			#pragma omp task in([c]qp1) out([c]xp1, [c]rp1, [c]up1) priority(1) label(alg7_fuse0)
			{
				/* x_i+1 = x_i + alpha_i * p_i */
				BLAS_cp(c, xp0, 1, xp1, 1);
				BLAS_axpy(c, alpha[i], pp0, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				BLAS_cp(c, rp0, 1, rp1, 1);
				BLAS_axpy(c, -1*alpha[i], sp0, 1, rp1, 1);

				/* u_i+1 = u_i - alpha_i * q_i */
				BLAS_cp(c, up0, 1, up1, 1);
				BLAS_axpy(c, -1*alpha[i], qp1, 1, up1, 1);
			}
		}

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait			
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		/* gamma[i+1] = <r[i+1], u[i+1]> */
		bblas_ddot(1, bm, 1, n, 1, r[i], u[i], &gamma[i]);
		/* w[i+1] = A * u[i+1] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]);

		#pragma omp taskwait on(gamma[i])
		beta[i] = gamma[i]/gamma[iprev];

		/* Axpy fuse p,s */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *pp0 = &(p[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *pp1 = &(p[i])[j];
			double *sp1 = &(s[i])[j];
			double *up1 = &(u[i])[j];
			double *wp1 = &(w[i])[j];

			#pragma omp task in([c]up1, [c]pp0, [c]wp1, [c]sp0) priority(1) label(alg7_fuse1)
			{
				/* p_i = u_i + beta_i * p_i-1 */
				BLAS_cp(c, up1, 1, pp1, 1);
				BLAS_axpy(c, beta[i], pp0, 1, pp1, 1);

				/* s_i = w_i + beta_i * s_i-1 */
				BLAS_cp(c, wp1, 1, sp1, 1);
				BLAS_axpy(c, beta[i], sp0, 1, sp1, 1);
			}
		}

		#pragma omp taskwait

		gamma[iprev] = delta = 0;
		stop_timer(&elapses[k]);

		//TODO Implement p-orthogonality check
#if 0
		BLAS_gemm(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, 1, 1, n, FP_ONE, p[i], n, s, n, FP_NOUGHT, &orth, 1);
		orth = FP_ABS(orth);
		if (isgreater(orth, porth * orth_fac)){
			fprintf(stderr, "orth fail %E %E\n", orth, porth*orth_fac);
			break;
		}
#endif

		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		residuals[k] = sr2norm;
		if ( isless(sr2norm, prec) ) {
			fprintf(stderr, "Precision reached\n");
			break;
		}
//		fprintf(stdout, "%d %E\n", k, residuals[k]);
	}

	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg7_ilu_bp.log", k, residuals, elapses);
	free(alpha);
	free(beta);
	free(gamma);
	free(residuals);
	free(elapses);
	return 0;
}

/* IFCG ILU */
int CG_IFCG_ILU_BP(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, int fuse, double orth_fac, int cglog_level, double *pool)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u_i+1 = u_i - alpha_i * q_i
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w_i+1 = w_i - alpha_i * z_i
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta_i * p_i
	offs += 2 * n;
	double *m[2] = {&pool[offs], &pool[offs+n]}; // m_i = M^-1 * w_i
	offs += 2 * n;
	double *n0[2] = {&pool[offs], &pool[offs+n]}; // n_i = A * m_i
	offs += 2 * n;
	double *z[2] = {&pool[offs], &pool[offs+n]}; // z_i = n_i + beta_i * z_i-1
	offs += 2 * n;
	double *q[2] = {&pool[offs], &pool[offs+n]}; // q_i = m_i + beta_i * q_i-1
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha_i * s_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i = w_i + beta_i * s_i-1

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
//	dcholsolv2_nested(1, bm, n, S, N, r[i], u[i]);
	dcholsolv2_blk(1, n, SA, NA, r[i], u[i]);
	/* w[0] = A * u[0] */
	hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]);

	start_timer();
	int kk = 0;
	int k;
	for ( k = 0; k < cgit; k++) {

		int iprev = i;
		i = i ^ 0x1;

		/* m[i] = M^-1 * w[i] */
		dcholsolv2_blk(1, n, SA, NA, w[iprev], m[i]);
//		dcholsolv2_nested(1, bm, n, S, N, w[iprev], m[i]);

		/* n[i] = A * m[i] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, m[i], fp_nought, n0[i]); 
		/* 
		 * gamma[i] = <r[i], u[i]>
		 * delta = <w[i], u[i]>
		 */
		cg_ddot2(1, bm, 1, n, 1, r[iprev], u[iprev], &gamma[i], w[iprev], u[iprev], &delta);

		#pragma omp taskwait on(gamma[i], delta)
		if ( k > 0 ) {
			beta[i] = gamma[i]/gamma[iprev];
			alpha[i] = gamma[i] / (delta - beta[i] * gamma[i] / alpha[iprev]);
		} else {
			beta[i] = (double) 0;
			alpha[i] = gamma[i]/delta;
		}
		gamma[iprev] = delta = 0;

		/* Grand fuse */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *zp0 = &(z[iprev])[j];
			double *qp0 = &(q[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *xp0 = &(x[iprev])[j];
			double *rp0 = &(r[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *zp1 = &(z[i])[j];
			double *qp1 = &(q[i])[j];
			double *sp1 = &(s[i])[j];
			double *pp1 = &(p[i])[j];
			double *xp1 = &(x[i])[j];
			double *rp1 = &(r[i])[j];
			double *up1 = &(u[i])[j];
			double *wp1 = &(w[i])[j];

			double *mp1 = &(m[i])[j];
			double *np1 = &(n0[i])[j];

			#pragma omp task out([c]zp1, [c]qp1, [c]sp1, [c]pp1, [c]xp1, [c]rp1, [c]up1, [c]wp1) \
				in([c]zp0, [c]qp0, [c]sp0, [c]pp0, [c]np1, [c]mp1, [c]up0, [c]xp0, [c]wp0, [c]rp0) \
				priority(1) label(alg4_fuse)
			{
				/* z_i = n_i + beta_i * z_i-1 */
				BLAS_cp(c, np1, 1, zp1, 1);
				BLAS_axpy(c, beta[i], zp0, 1, zp1, 1);

				/* q_i = m_i + beta_i * q_i-1 */
				BLAS_cp(c, mp1, 1, qp1, 1);
				BLAS_axpy(c, beta[i], qp0, 1, qp1, 1);

				/* s_i = w_i + beta_i * s_i-1 */
				BLAS_cp(c, wp0, 1, sp1, 1);
				BLAS_axpy(c, beta[i], sp0, 1, sp1, 1);

				/* p_i = u_i + beta_i * p_i-1 */
				BLAS_cp(c, up0, 1, pp1, 1);
				BLAS_axpy(c, beta[i], pp0, 1, pp1, 1);

				/* x_i+1 = x_i + alpha_i * p_i */
				BLAS_cp(c, xp0, 1, xp1, 1);
				BLAS_axpy(c, alpha[i], pp1, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				BLAS_cp(c, rp0, 1, rp1, 1);
				BLAS_axpy(c, -1*alpha[i], sp1, 1, rp1, 1);

				/* u_i+1 = u_i - alpha_i * q_i */
				BLAS_cp(c, up0, 1, up1, 1);
				BLAS_axpy(c, -1*alpha[i], qp1, 1, up1, 1);

				/* w_i+1 = w_i - alpha_i * z_i */
				BLAS_cp(c, wp0, 1, wp1, 1);
				BLAS_axpy(c, -1*alpha[i], zp1, 1, wp1, 1);
			}
		}

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		if ( k > 0 && k % fuse == 0 ) {
			#pragma omp taskwait
			stop_timer(&elapses[kk]);
			double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
			double sr2norm = norm_r/norm_b;
			residuals[kk] = sr2norm;
			if ( isless(sr2norm, prec) ) {
				fprintf(stderr, "Precision reached\n");
				break;
			}
			kk += 1;
			start_timer();
		}
	}

	#pragma omp taskwait
	if ( k == cgit ) {
		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		stop_timer(&elapses[kk]);
		residuals[kk] = sr2norm;
	}

	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg4_ifcg_ilu_bp.log", kk, residuals, elapses);
	free(alpha);
	free(beta);
	free(gamma);
	free(residuals);
	free(elapses);
	return 0;
}

/* IFCG V2 BP*/
int CG_IFCG_V2_ILU_BP(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, int fuse, double orth_fac, int cglog_level, double *pool)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u_i+1 = u_i - alpha_i * q_i
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w_i+1 = w_i - alpha_i * z_i
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta_i * p_i
	offs += 2 * n;
	double *m[2] = {&pool[offs], &pool[offs+n]}; // m_i = M^-1 * w_i
	offs += 2 * n;
	double *n0[2] = {&pool[offs], &pool[offs+n]}; // n_i = A * m_i
	offs += 2 * n;
	double *z[2] = {&pool[offs], &pool[offs+n]}; // z_i = n_i + beta_i * z_i-1
	offs += 2 * n;
	double *q[2] = {&pool[offs], &pool[offs+n]}; // q_i = m_i + beta_i * q_i-1
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha_i * s_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i = w_i + beta_i * s_i-1

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
	dcholsolv2_blk(1, n, SA, NA, r[i], u[i]);
//	dcholsolv2_nested(1, bm, n, S, N, r[i], u[i]);
//	bsblas_dcholsolv2(1, bm, n, S, N, r[i], u[i]);
	/* w[0] = A * u[0] */
	hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]);

	start_timer();
	int kk = 0;
	int k;
	for ( k = 0; k < cgit; k++) {

		int iprev = i;
		i = i ^ 0x1;

		/* m[i] = M^-1 * w[i] */
		dcholsolv2_blk(1, n, SA, NA, w[iprev], m[i]);
//		dcholsolv2_nested(1, bm, n, S, N, w[iprev], m[i]);
//		bsblas_dcholsolv2(1, bm, n, S, N, w[iprev], m[i]);

		/* 
		 * gamma[i] = <r[i], u[i]>
		 */
		bblas_ddot(1, bm, 1, n, 1, r[iprev], u[iprev], &gamma[i]);

		#pragma omp taskwait on(gamma[i])
		if ( k > 0 ) {
			beta[i] = gamma[i]/gamma[iprev];
//			alpha[i] = gamma[i] / (delta - beta[i] * gamma[i] / alpha[iprev]);
		} else {
			beta[i] = (double) 0;
//			alpha[i] = gamma[i]/delta;
		}
//		gamma[iprev] = delta = 0;

		/* 
		 * s_i = w_i + beta_i * s_i-1
		 * p_i = u_i + beta_i * p_i-1
		 */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *sp0 = &(s[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *sp1 = &(s[i])[j];
			double *pp1 = &(p[i])[j];

			#pragma omp task out([c]sp1, [c]pp1) in([c]sp0, [c]pp0, [c]up0, [c]wp0) priority(1) label(alg4_apxy)
			{
				/* s_i = w_i + beta_i * s_i-1 */
				BLAS_cp(c, wp0, 1, sp1, 1);
				BLAS_axpy(c, beta[i], sp0, 1, sp1, 1);

				/* p_i = u_i + beta_i * p_i-1 */
				BLAS_cp(c, up0, 1, pp1, 1);
				BLAS_axpy(c, beta[i], pp0, 1, pp1, 1);
			}
		}

		/* n[i] = A * m[i] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, m[i], fp_nought, n0[i]); 

		bblas_ddot(1, bm, 1, n, 1, w[iprev], u[iprev], &delta);
//		cg_ddot2(1, bm, 1, n, 1, r[iprev], u[iprev], &gamma[i], w[iprev], u[iprev], &delta);
		#pragma omp taskwait on(delta)
		if ( k > 0 ) {
			alpha[i] = gamma[i] / (delta - beta[i] * gamma[i] / alpha[iprev]);
		} else {
			alpha[i] = gamma[i]/delta;
		}
		gamma[iprev] = delta = 0;

		/* Grand fuse */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *zp0 = &(z[iprev])[j];
			double *qp0 = &(q[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *xp0 = &(x[iprev])[j];
			double *rp0 = &(r[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *zp1 = &(z[i])[j];
			double *qp1 = &(q[i])[j];
			double *sp1 = &(s[i])[j];
			double *pp1 = &(p[i])[j];
			double *xp1 = &(x[i])[j];
			double *rp1 = &(r[i])[j];
			double *up1 = &(u[i])[j];
			double *wp1 = &(w[i])[j];

			double *mp1 = &(m[i])[j];
			double *np1 = &(n0[i])[j];

			#pragma omp task out([c]zp1, [c]qp1, [c]xp1, [c]rp1, [c]up1, [c]wp1) \
				in([c]zp0, [c]qp0, [c]sp0, [c]pp0, [c]np1, [c]mp1, [c]up0, [c]xp0, [c]wp0, [c]rp0) \
				priority(1) label(alg4_fuse)
			{
				/* z_i = n_i + beta_i * z_i-1 */
				BLAS_cp(c, np1, 1, zp1, 1);
				BLAS_axpy(c, beta[i], zp0, 1, zp1, 1);

				/* q_i = m_i + beta_i * q_i-1 */
				BLAS_cp(c, mp1, 1, qp1, 1);
				BLAS_axpy(c, beta[i], qp0, 1, qp1, 1);

				/* x_i+1 = x_i + alpha_i * p_i */
				BLAS_cp(c, xp0, 1, xp1, 1);
				BLAS_axpy(c, alpha[i], pp1, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				BLAS_cp(c, rp0, 1, rp1, 1);
				BLAS_axpy(c, -1*alpha[i], sp1, 1, rp1, 1);

				/* u_i+1 = u_i - alpha_i * q_i */
				BLAS_cp(c, up0, 1, up1, 1);
				BLAS_axpy(c, -1*alpha[i], qp1, 1, up1, 1);

				/* w_i+1 = w_i - alpha_i * z_i */
				BLAS_cp(c, wp0, 1, wp1, 1);
				BLAS_axpy(c, -1*alpha[i], zp1, 1, wp1, 1);
			}
		}

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		if ( k > 0 && k % fuse == 0 ) {
			#pragma omp taskwait
			stop_timer(&elapses[kk]);
			double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
			double sr2norm = norm_r/norm_b;
			residuals[kk] = sr2norm;
			if ( isless(sr2norm, prec) ) {
				fprintf(stderr, "Precision reached\n");
				break;
			}
			kk += 1;
			start_timer();
		}
	}

	#pragma omp taskwait
	if ( k == cgit ) {
		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		stop_timer(&elapses[kk]);
		residuals[kk] = sr2norm;
	}

	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg4_ifcg_v2_ilu_bp.log", kk, residuals, elapses);
	free(alpha);
	free(beta);
	free(gamma);
	free(residuals);
	free(elapses);
	return 0;
}


/* Pipelined PCG V3*/
int CG_ALG4_V4(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, int fuse, double orth_fac, int cglog_level)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;

	int offs = 0;
	double *pool = calloc(10 * 2 * n, sizeof(double));
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u_i+1 = u_i - alpha_i * q_i
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w_i+1 = w_i - alpha_i * z_i
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta_i * p_i
	offs += 2 * n;
	double *m[2] = {&pool[offs], &pool[offs+n]}; // m_i = M^-1 * w_i
	offs += 2 * n;
	double *n0[2] = {&pool[offs], &pool[offs+n]}; // n_i = A * m_i
	offs += 2 * n;
	double *z[2] = {&pool[offs], &pool[offs+n]}; // z_i = n_i + beta_i * z_i-1
	offs += 2 * n;
	double *q[2] = {&pool[offs], &pool[offs+n]}; // q_i = m_i + beta_i * q_i-1
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha_i * s_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i = w_i + beta_i * s_i-1

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(double));

	int i = 0;
	/* r[0] = b */
	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
	bsblas_dcholsolv2(1, bm, n, S, N, r[i], u[i]);
	/* w[0] = A * u[0] */
	hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]);

	start_timer();
	int kk = 0;
	int k;
	for ( k = 0; k < cgit; k++) {

		int iprev = i;
		i = i ^ 0x1;

		/* m[i] = M^-1 * w[i] */
		bsblas_dcholsolv2(1, bm, n, S, N, w[iprev], m[i]);
		/* n[i] = A * m[i] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, m[i], fp_nought, n0[i]); 
		/* 
		 * gamma[i] = <r[i], u[i]>
		 * delta = <w[i], u[i]>
		 */
		cg_ddot2_commutative(1, bm, 1, n, 1, r[iprev], u[iprev], &gamma[i], w[iprev], u[iprev], &delta);

		#pragma omp taskwait on(gamma[i], delta)
		if ( k > 0 ) {
			beta[i] = gamma[i]/gamma[iprev];
			alpha[i] = gamma[i] / (delta - beta[i] * gamma[i] / alpha[iprev]);
		} else {
			beta[i] = (double) 0;
			alpha[i] = gamma[i]/delta;
		}
		gamma[iprev] = delta = 0;

		/* Grand fuse */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *zp0 = &(z[iprev])[j];
			double *qp0 = &(q[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *xp0 = &(x[iprev])[j];
			double *rp0 = &(r[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *zp1 = &(z[i])[j];
			double *qp1 = &(q[i])[j];
			double *sp1 = &(s[i])[j];
			double *pp1 = &(p[i])[j];
			double *xp1 = &(x[i])[j];
			double *rp1 = &(r[i])[j];
			double *up1 = &(u[i])[j];
			double *wp1 = &(w[i])[j];

			double *mp1 = &(m[i])[j];
			double *np1 = &(n0[i])[j];

			#pragma omp task out([c]zp1, [c]qp1, [c]sp1, [c]pp1, [c]xp1, [c]rp1, [c]up1, [c]wp1) \
				in([c]zp0, [c]qp0, [c]sp0, [c]pp0, [c]np1, [c]mp1, [c]up0, [c]xp0, [c]wp0, [c]rp0) \
				priority(1) label(alg4_fuse)
			{
				/* z_i = n_i + beta_i * z_i-1 */
				BLAS_cp(c, np1, 1, zp1, 1);
				BLAS_axpy(c, beta[i], zp0, 1, zp1, 1);

				/* q_i = m_i + beta_i * q_i-1 */
				BLAS_cp(c, mp1, 1, qp1, 1);
				BLAS_axpy(c, beta[i], qp0, 1, qp1, 1);

				/* s_i = w_i + beta_i * s_i-1 */
				BLAS_cp(c, wp0, 1, sp1, 1);
				BLAS_axpy(c, beta[i], sp0, 1, sp1, 1);

				/* p_i = u_i + beta_i * p_i-1 */
				BLAS_cp(c, up0, 1, pp1, 1);
				BLAS_axpy(c, beta[i], pp0, 1, pp1, 1);

				/* x_i+1 = x_i + alpha_i * p_i */
				BLAS_cp(c, xp0, 1, xp1, 1);
				BLAS_axpy(c, alpha[i], pp1, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				BLAS_cp(c, rp0, 1, rp1, 1);
				BLAS_axpy(c, -1*alpha[i], sp1, 1, rp1, 1);

				/* u_i+1 = u_i - alpha_i * q_i */
				BLAS_cp(c, up0, 1, up1, 1);
				BLAS_axpy(c, -1*alpha[i], qp1, 1, up1, 1);

				/* w_i+1 = w_i - alpha_i * z_i */
				BLAS_cp(c, wp0, 1, wp1, 1);
				BLAS_axpy(c, -1*alpha[i], zp1, 1, wp1, 1);
			}
		}

		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}

		if ( k > 0 && k % fuse == 0 ) {
			#pragma omp taskwait
			stop_timer(&elapses[kk]);
			double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
			double sr2norm = norm_r/norm_b;
			residuals[kk] = sr2norm;
			if ( isless(sr2norm, prec) ) {
				fprintf(stderr, "Precision reached\n");
				break;
			}
			kk += 1;
			start_timer();
		}
	}

	#pragma omp taskwait
	if ( k == cgit ) {
		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		stop_timer(&elapses[kk]);
		residuals[kk] = sr2norm;
	}

	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info("cg_alg4_v3.log", kk, residuals, elapses);
	free(pool);
	free(alpha);
	free(beta);
	free(gamma);
	free(residuals);
	free(elapses);
	return 0;
}

#define STRIDE 10
#define SFUSE 10
/* Pipelined PCG Auto-tuning */
int CG_ALG4_AT(void *A, void *Ahbh, double *solution, double *b, int cgit, int bm, double prec, int ACCIMP, int fuse, double orth_fac, int cglog_level)
{
	hbmat_t *Ahb = (hbmat_t*) A;
	int n = Ahb->m;
	int FUSE = SFUSE;
	int HALT = FUSE;
	int trial = 1;

	int offs = 0;
	double *pool = calloc(10 * 2 * n, sizeof(double));
	double *x[2] = {&pool[offs], &pool[offs+n]};
	offs += 2 * n;
	double *u[2] = {&pool[offs], &pool[offs+n]};   // u_i+1 = u_i - alpha_i * q_i
	offs += 2 * n;
	double *w[2] = {&pool[offs], &pool[offs+n]};   // w_i+1 = w_i - alpha_i * z_i
	offs += 2 * n;
	double *p[2] = {&pool[offs], &pool[offs+n]}; // p_i+1 = u_i+1 + beta_i * p_i
	offs += 2 * n;
	double *m[2] = {&pool[offs], &pool[offs+n]}; // m_i = M^-1 * w_i
	offs += 2 * n;
	double *n0[2] = {&pool[offs], &pool[offs+n]}; // n_i = A * m_i
	offs += 2 * n;
	double *z[2] = {&pool[offs], &pool[offs+n]}; // z_i = n_i + beta_i * z_i-1
	offs += 2 * n;
	double *q[2] = {&pool[offs], &pool[offs+n]}; // q_i = m_i + beta_i * q_i-1
	offs += 2 * n;
	double *r[2] = {&pool[offs], &pool[offs+n]}; // r_i+1 = r_i - alpha_i * s_i
	offs += 2 * n;
	double *s[2] = {&pool[offs], &pool[offs+n]}; // s_i = w_i + beta_i * s_i-1

	double *alpha = calloc(2, sizeof(double));
	double *beta = calloc(2, sizeof(double));
	double *gamma = calloc(2, sizeof(double));
	double delta = (double) 0;

	double orth;
	double porth = DBL_MAX;
	double norm_b = cblas_ddot(n, b, 1, b, 1);

	double *residuals = malloc(cgit * sizeof(double));
	unsigned int *elapses = malloc(cgit * sizeof(unsigned int));
	int *fuses = calloc(cgit, sizeof(int));
	int *iters = calloc(cgit, sizeof(int));

	int i = 0;
	/* r[0] = b */
	bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
	/* r[0] = b - A * x[0] */
	hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
	/* u[0] = M^-1 * r[0] */
	bsblas_dcholsolv2(1, bm, n, S, N, r[i], u[i]);
	/* w[0] = A * u[0] */
	hbsblas_dcsrmv(1, bm, fp_one, Ahbh, u[i], fp_nought, w[i]);

	start_timer();
	int kk = 0;
	int k;
	for ( k = 0; k < cgit; k++) {

		int iprev = i;
		i = i ^ 0x1;

		/* m[i] = M^-1 * w[i] */
		bsblas_dcholsolv2(1, bm, n, S, N, w[iprev], m[i]);
		/* n[i] = A * m[i] */
		hbsblas_dcsrmv(1, bm, fp_one, Ahbh, m[i], fp_nought, n0[i]); 
		/* 
		 * gamma[i] = <r[i], u[i]>
		 * delta = <w[i], u[i]>
		 */
		cg_ddot2(1, bm, 1, n, 1, r[iprev], u[iprev], &gamma[i], w[iprev], u[iprev], &delta);

		#pragma omp taskwait on(gamma[i], delta)
		if ( k > 0 ) {
			beta[i] = gamma[i]/gamma[iprev];
			alpha[i] = gamma[i] / (delta - beta[i] * gamma[i] / alpha[iprev]);
		} else {
			beta[i] = (double) 0;
			alpha[i] = gamma[i]/delta;
		}
		gamma[iprev] = delta = 0;

		/* Grand fuse */
		for (int j = 0; j < n; j += bm ) {
			int cs = n - j;
			int c = cs < bm ? cs : bm;
			double *zp0 = &(z[iprev])[j];
			double *qp0 = &(q[iprev])[j];
			double *sp0 = &(s[iprev])[j];
			double *pp0 = &(p[iprev])[j];
			double *xp0 = &(x[iprev])[j];
			double *rp0 = &(r[iprev])[j];
			double *up0 = &(u[iprev])[j];
			double *wp0 = &(w[iprev])[j];

			double *zp1 = &(z[i])[j];
			double *qp1 = &(q[i])[j];
			double *sp1 = &(s[i])[j];
			double *pp1 = &(p[i])[j];
			double *xp1 = &(x[i])[j];
			double *rp1 = &(r[i])[j];
			double *up1 = &(u[i])[j];
			double *wp1 = &(w[i])[j];

			double *mp1 = &(m[i])[j];
			double *np1 = &(n0[i])[j];

			#pragma omp task out([c]zp1, [c]qp1, [c]sp1, [c]pp1, [c]xp1, [c]rp1, [c]up1, [c]wp1) \
				in([c]zp0, [c]qp0, [c]sp0, [c]pp0, [c]np1, [c]mp1, [c]up0, [c]xp0, [c]wp0, [c]rp0) \
				priority(1) label(alg4_fuse)
			{
				/* z_i = n_i + beta_i * z_i-1 */
				BLAS_cp(c, np1, 1, zp1, 1);
				BLAS_axpy(c, beta[i], zp0, 1, zp1, 1);

				/* q_i = m_i + beta_i * q_i-1 */
				BLAS_cp(c, mp1, 1, qp1, 1);
				BLAS_axpy(c, beta[i], qp0, 1, qp1, 1);

				/* s_i = w_i + beta_i * s_i-1 */
				BLAS_cp(c, wp0, 1, sp1, 1);
				BLAS_axpy(c, beta[i], sp0, 1, sp1, 1);

				/* p_i = u_i + beta_i * p_i-1 */
				BLAS_cp(c, up0, 1, pp1, 1);
				BLAS_axpy(c, beta[i], pp0, 1, pp1, 1);

				/* x_i+1 = x_i + alpha_i * p_i */
				BLAS_cp(c, xp0, 1, xp1, 1);
				BLAS_axpy(c, alpha[i], pp1, 1, xp1, 1);

				/* r_i+1 = r_i - alpha_i * s_i */
				BLAS_cp(c, rp0, 1, rp1, 1);
				BLAS_axpy(c, -1*alpha[i], sp1, 1, rp1, 1);

				/* u_i+1 = u_i - alpha_i * q_i */
				BLAS_cp(c, up0, 1, up1, 1);
				BLAS_axpy(c, -1*alpha[i], qp1, 1, up1, 1);

				/* w_i+1 = w_i - alpha_i * z_i */
				BLAS_cp(c, wp0, 1, wp1, 1);
				BLAS_axpy(c, -1*alpha[i], zp1, 1, wp1, 1);
			}
		}

#if 0
		/* Accuracy improvement */
		if ( k > 0 && k % ACCIMP == 0 ) {
			#pragma omp taskwait
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
		}
#endif

		if ( k == HALT ) {
			#pragma omp taskwait
			/* Accuracy improvement */
			bblas_dcopy(1, bm, 1, n, 1, b, r[i]);
			hbsblas_dcsrmv(1, bm, fp_mone, Ahbh, x[i], fp_one, r[i]);
			#pragma omp taskwait

			stop_timer(&elapses[kk]);
			double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
			double sr2norm = norm_r/norm_b;
			residuals[kk] = sr2norm;
			fuses[kk] = FUSE;
			iters[kk] = k;
			if ( isless(sr2norm, prec) ) {
				fprintf(stderr, "Precision reached\n");
				break;
			}

			/* Fuse auto-tuning */
			if ( trial ) {
				if ( kk == 0 ) {
					FUSE += STRIDE;
				} else {
					int elp_inc = ( elapses[kk] <= elapses[kk-1] );
					int res_inc = ( isless(residuals[kk], residuals[kk-1]) );
					if (elp_inc && res_inc) {
						FUSE += STRIDE;
					} else {
						FUSE -= STRIDE;
						trial = 0;
						fprintf(stderr, "Trial turned off\n");
					}
				}
			}
			HALT += FUSE;
//			printf("k: %d kk: %d FUSE: %d\n", k, kk, FUSE);
//			printf("elp[%d] %d elp[%d] %d\n", kk, elapses[kk], kk-1, elapses[kk-1]);

			kk += 1;
			start_timer();
		}
	}

	#pragma omp taskwait
	if ( k == cgit ) {
		double norm_r = sqrt(cblas_ddot(n, r[i], 1, r[i], 1));
		double sr2norm = norm_r/norm_b;
		stop_timer(&elapses[kk]);
		residuals[kk] = sr2norm;
		fuses[kk] = FUSE;
	}

	memcpy(solution, x[i], n * sizeof(double));
	if ( cglog_level )
		dump_info4("cg_alg4_at.log", kk, residuals, elapses, fuses, iters);
	free(pool);
	free(alpha);
	free(beta);
	free(gamma);
	free(residuals);
	free(elapses);
	return 0;
}

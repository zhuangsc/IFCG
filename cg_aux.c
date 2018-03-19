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


#include "cg_aux.h"

const char *scan_dconspec = "%lf";
const char *scan_sconspec = "%f";

void hb_read_double(char *input_file, int *m, int *n, int *elemc, int **vptr, int **vpos, double **vval)
{
  double *exact = NULL;
  double *guess = NULL;
  int i;
  int indcrd;
  char *indfmt = NULL;
  FILE *input;
  int j;
  char *key = NULL;
  int khi;
  int klo;
  char *mxtype = NULL;
  int neltvl;
  int nrhs;
  int nrhsix;
  int ptrcrd;
  char *ptrfmt = NULL;
  int rhscrd;
  char *rhsfmt = NULL;
  int *rhsind = NULL;
  int *rhsptr = NULL;
  char *rhstyp = NULL;
  double *rhsval = NULL;
  double *rhsvec = NULL;
  char *title = NULL;
  int totcrd;
  int valcrd;
  char *valfmt = NULL;
  int nrow;
  int ncol;
  int nnzero;
  int *colptr = NULL;
  int *rowind = NULL;
  double *values = NULL;

  input = fopen ( input_file, "rt" );

  if ( !input )
  {
    printf ( "  Error opening the file.\n" );
    return;
  }

  hb_file_read ( input, &title, &key, &totcrd, &ptrcrd, &indcrd, 
    &valcrd, &rhscrd, &mxtype, &nrow, &ncol, &nnzero, &neltvl, 
    &ptrfmt, &indfmt, &valfmt, &rhsfmt, &rhstyp, &nrhs, &nrhsix, 
    &colptr, &rowind, &values, &rhsval, &rhsptr, &rhsind, &rhsvec, 
    &guess, &exact );

  fclose ( input );

  if ( exact )
  {
    free ( exact );
  }
  if ( guess )
  {
    free ( guess );
  }
  if ( rhsind )
  {
    free ( rhsind );
  }
  if ( rhsptr )
  {
    free ( rhsptr );
  }
  if ( rhsval )
  {
    free ( rhsval );
  }
  if ( rhsvec )
  {
    free ( rhsvec );
  }

  *m = nrow;
  *n = ncol;
  *elemc = nnzero;
  *vptr = colptr;
  *vpos = rowind;
  *vval = values;

  return;
}

void hb_reset(hbmat_t *A)
{
	A->m = 0; A->n =0; A->elemc = 0;
	A->vptr = 0; A->vpos = 0; A->vval = 0;
	A->vdiag = NULL; 
	A->b = 0; 
	A->trans = 0; A->orig = 0; A->hyper = 0;
	A->orig_row = 0; A->orig_col = 0; A->e_tree = 0;
	A->type = 0;
	A->FACT = 0;
}

void one2zero(hbmat_t* in_matrix)
{
	int m = in_matrix->m;
	int elemc = in_matrix->elemc;
	int *vptr = in_matrix->vptr;
	int *vpos = in_matrix->vpos;

	int i;
	for ( i = 0; i <= m; i++ ) {
		vptr[i]--;
	}
	for ( i=0; i<elemc; i++ ) {
		vpos[i]--;
	}
}

typedef struct _sparse_nodes {
	int row;
	int col;
	double val;
	struct sparse_node *next;
	struct sparse_node *current;
} _sn_t;

/*	Expand an symmetric half matrix B to its full counterpart A */
void hb_sym_expand(hbmat_t *A, hbmat_t *B)
{
	hb_init_basic(A, B);
	int m = A->m;
	A->elemc = B->elemc * 2 - m;
	int nnz = A->elemc;
	A->vptr = malloc((m+1) * sizeof(int));
	A->vpos = malloc(nnz * sizeof(int));
	A->vval = malloc(nnz * sizeof(double));

	int *vptra = A->vptr; int *vposa = A->vpos; double *vvala = A->vval;
	int *vptrb = B->vptr; int *vposb = B->vpos; double *vvalb = B->vval;

	_sn_t *ll_mat = malloc(m * sizeof(_sn_t));

	int i;
	for ( i = 0; i < m; i++ ) {
		ll_mat[i].row = m;
		ll_mat[i].col = -1;
		ll_mat[i].val = 0.0;
		ll_mat[i].next = NULL;
		ll_mat[i].current = &ll_mat[i];
	}

	int vptr_c = 0;
	int elemc_c = 0;

	for ( i = 0; i < m; i++ ) {
		vptra[i] = vptr_c;
		int bptr = vptrb[i]; 
		int eptr = vptrb[i+1];

		/* Fill the lower csr info */
		_sn_t *c = &ll_mat[i];
		while ( c != NULL && c->col != -1 ) {
			vposa[elemc_c] = c->col;
			vvala[elemc_c] = c->val;
			elemc_c++;
			vptr_c++;
			c = c->next;
		}
					
		/* Copy the upper csr info */
		int j;
		for ( j = bptr; j < eptr; j++ ) {
			int col = vposb[j];
			double val = vvalb[j];
			vposa[elemc_c] = col;
			vvala[elemc_c] = val;
			elemc_c++;
			/*-------------------------------------------------*
			 * linked list insert
			 *-------------------------------------------------*/
			_sn_t *head = &ll_mat[col];
			_sn_t *current = ll_mat[col].current;
			if ( current->col != -1 ) {
				current->next = malloc(sizeof(_sn_t));
				current = current->next;
				head->current = current;
			}
			current->row = col;
			current->col = i;
			current->val = val;
			current->next = NULL;
			current->current = current;
			/*-------------------------------------------------*
			 * linked list insert end
			 *-------------------------------------------------*/
		}
		vptr_c += eptr - bptr;
	}

	vptra[m] = vptr_c;

	for ( i = 0; i < m; i++ ) {
		_sn_t *c= ll_mat[i].next;
		while ( c != NULL ) {
			_sn_t *n = c->next;
			free(c);
			c = n;
		}
	}
	free(ll_mat);
}

/* Copy basic info from B to A */
void hb_init_basic(hbmat_t *A, hbmat_t *B)
{
	hb_reset(A);
	int M = B->m;
	int elemc = B->elemc;
	A->m = A->n = M;
	A->elemc = elemc;
}

void hb_free(hbmat_t *A)
{
	free(A->vptr); 
	free(A->vpos);
	free(A->vval);
	free(A->vdiag);
	free(A->e_tree);

	free(A);
}

void* __hb2hbh_block(int I, int J, hbmat_t *A, int b, hbmat_t *Bp) 
{
	int alloc = Bp == NULL;
	
	if ( b < 0 ) {
		fprintf(stderr, "err: b must be positive\n");
		return NULL;
	}

	int m = A->m; int n = A->n;
	int* vptr = A->vptr; 
	int* vpos = A->vpos; 
	double* vval = A->vval;
	int offs = vptr[0] == 0 ? 0 : 1;
	int csr = 1;//hb_CSR(A);

	int brow = I*b;  
	int bcol = J*b;
	int rleft = m - brow;
	int cleft = n - bcol;
	int rows = b > rleft ? rleft : b;
	int cols = b > cleft ? cleft : b;
	int erow = brow + rows;
	int ecol = bcol + cols;
	int dimb = csr ? brow : bcol;
	int dime = csr ? erow : ecol;
	int rngb = csr ? bcol : brow;
	int rnge = csr ? ecol : erow;

	vector_t* ab_vptr = vector_create();
	vector_t* ab_vpos = vector_create();
	vector_t* ab_vval = vector_create();
	vel_t vel;

	int L;
	for ( L = dimb; L < dime; ++L ) {
		vel.i = ab_vpos->elemc + offs; 
		vector_insert(ab_vptr, vel);
	
		int k;
		for ( k = vptr[L]; k < vptr[L+1]; ++k ) {
			int lk = k - offs;
			int c = vpos[lk] - offs; 

			if ( c >= rngb && c < rnge ) {
				vel.i = c - rngb + offs;
				vector_insert(ab_vpos, vel);
				vel.d = vval[lk];
				vector_insert(ab_vval, vel);
			}
		}
	}

	vel.i = ab_vpos->elemc + offs;
	vector_insert(ab_vptr, vel);

	if ( alloc ) {
		Bp = malloc(sizeof(hbmat_t));
		hb_reset(Bp);
	}

	if ( ab_vpos->elemc ) {
		Bp->m = rows; 
		Bp->n = cols; 
		Bp->elemc = ab_vpos->elemc;
		Bp->vdiag = NULL;
		Bp->vptr = vector2int(ab_vptr); 
		Bp->vpos = vector2int(ab_vpos);
		Bp->vval = vector2double(ab_vval);
	} else {
		vector_free(ab_vptr);
		vector_free(ab_vpos);
		vector_free(ab_vval);

		return NULL;
	}

	return Bp;
}

hbmat_t* hb2hbh(hbmat_t *A, int b, int is_csr) 
{
	int m = A->m; 
	int n = A->n; 
	int elemc = A->elemc;
	int *vptr = A->vptr; 
	int *vpos = A->vpos; 
	double* vval = A->vval;
	int M = (m+b-1) / b;
	int N = (n+b-1) / b;
	int num = M * N;
	int offs = vptr[0] == 0 ? 0 : 1;

	hbmat_t* hyper = malloc(sizeof(hbmat_t));
	hb_reset(hyper);
	hyper->m = M; hyper->n = N; hyper->vdiag = NULL;
	hyper->orig = A;
	hyper->vval = malloc(num * sizeof(hbmat_t*));
	hbmat_t** hbmat_array = malloc(num * sizeof(hbmat_t*));

	vector_t* ab_vptr = vector_create(); 
	vector_t* ab_vpos = vector_create();
	vel_t pos_val;

	int acc0 = 0;
	int acc = 0;
	int I, J;
	if ( is_csr ) {
		for ( I = 0; I < M; ++I ) {
			pos_val.i = ab_vpos->elemc + offs;
			vector_insert(ab_vptr, pos_val);
			for ( J = 0; J < N; ++J ) {
				hbmat_t *B = __hb2hbh_block(I, J, A, b, NULL);  
				if ( B != NULL ) {
					pos_val.i = J + offs;
					vector_insert(ab_vpos, pos_val);
					((hbmat_t**)hyper->vval)[acc0] = B;
					++acc0;
				}
				++acc;
			}
		}
	} else {
		printf("warn: hb2hbh for csc not yet implemented\n");
	}

	pos_val.i = ab_vpos->elemc + offs;
	vector_insert(ab_vptr, pos_val);
	hyper->elemc = ab_vpos->elemc;
	hyper->vptr = vector2int(ab_vptr);
	hyper->vpos = vector2int(ab_vpos);

//	hb_setdiag(hyper);

	return hyper;
}

/* Construct an array of block diagonal submatrices */
void hb_sym_diag_block(hbmat_t *src_mat, int bsze, hbmat_t *diagb)
{
	/* Assuming CSR */
	int m = src_mat->m;
	/* Number of subblocks */
	int bs = (m+bsze-1)/bsze;
	int *svptr = src_mat->vptr; int *svpos = src_mat->vpos;
	double *svval = src_mat->vval;
	int i;
	/* Loop for generating all the diagonal blocks*/
	for ( i = 0; i < bs; i++ ) {
		hbmat_t *d = &diagb[i];
		int elemc = 0;
		int brow = i*bsze; int erow = brow+bsze;
		erow = erow > m ? m : erow;
		int dim = erow - brow;
		d->m = d->n = dim;
		// Allocate individual HB structures
		// Note that vpos and vval size are over-estimated
		int *vptr = malloc((dim+1) * sizeof(int));
		int esze = (svptr[erow] - svptr[brow]);
		int *vpos = malloc(esze * sizeof(int));
		double *vval = malloc(esze * sizeof(double));
		int idx;
		int row;
		/* Traverse through rows */
		for ( row = brow, idx = 0; row < erow; row++ ,idx++) {
			vptr[idx] = elemc;
			int pos = svptr[row]; int epos = svptr[row+1];
			while ( pos < epos ) {
				int col = svpos[pos];
				/* Only take the lower triangular part of the matrix */
//				if ( col >= row && col < erow ) { //Upper
//				if ( col >= brow && col < row ) { //Lower
				if ( col >= brow && col < erow ) { //Complete
					vpos[elemc] = col - brow;
					vval[elemc] = svval[pos];
					elemc++;
				}
				pos++; 
			}
		}
		vptr[idx] = elemc;
		d->elemc = elemc;
		d->vptr = vptr;
		//FIXME using realloc to reduce memory consumption
		d->vpos = vpos;
		d->vval = vval;
		//TODO Remove verifications
//		hb_sanity_check("A_hb", d, 0);
//		assert(idx == dim);
//		assert(d->vptr != NULL && d->vpos != NULL && d->vval != NULL);
	}
}

/* Block diagonal (non-split) */
void hb_sym_diag(hbmat_t *src_mat, int bsze, hbmat_t *d)
{
	/* Assuming CSR */
	int m = src_mat->m;
	int bs = (m+bsze-1)/bsze;
	int *svptr = src_mat->vptr; 
	int *svpos = src_mat->vpos;
	double *svval = src_mat->vval;
	d->m = d->n = m;
	int elemc = d->elemc = 0;
	d->vptr = malloc((m+1) * sizeof(int));
	d->vpos = malloc(src_mat->elemc * sizeof(int));
	d->vval = malloc(src_mat->elemc * sizeof(double));
	int *dvptr = d->vptr;
	int *dvpos = d->vpos;
	double *dvval = d->vval;
	int idx = 0;

	int i;
	for ( i = 0; i < bs; i++ ) {
		int brow = i*bsze; 
		int erow = brow+bsze;
		erow = erow > m ? m : erow;
		int dim = erow - brow;
//		int esze = (svptr[erow] - svptr[brow]);
//		int idx;
		int row;
		/* Traverse through rows */
		for ( row = brow; row < erow; row++) {
			dvptr[idx] = elemc;
			idx += 1;
			int pos = svptr[row]; 
			int epos = svptr[row+1];
			while ( pos < epos ) {
				int col = svpos[pos];
				/* Only take the lower triangular part of the matrix */
//				if ( col >= row && col < erow ) { //Upper
//				if ( col >= brow && col < row ) { //Lower
				if ( col >= brow && col < erow ) { //Complete
					//TODO Verify
//					vpos[elemc] = col - brow;
					dvpos[elemc] = col;
					dvval[elemc] = svval[pos];
					elemc++;
				}
				pos++; 
			}
		}
		dvptr[idx] = elemc;
		d->elemc = elemc;
		d->vptr = dvptr;
		d->vpos = dvpos;
		d->vval = dvval;
	}
//	printf("m %d n %d elemc : %d\n", d->m, d->n, d->elemc);
//	for(int i = 0; i < m; i++ ){
//		printf("[%d]: %d ", i, dvptr[i]);
//	}
//	printf("\n\n");
//	for(int i = 0; i < elemc; i++ ){
//		printf("[%d]: %d ", i, dvpos[i]);
//	}
//	printf("\n\n");
//	for(int i = 0; i < elemc; i++ ){
//		printf("[%d]: %E ", i, dvval[i]);
//	}
//	printf("\n\n");
}


int read_mm2dense(FILE *f, int m, int n, double *A) 
{
	char buf[1024];
	double el;

	int i = 0;
	while ( fgets(buf, sizeof(buf), f) != NULL && i < m) {
    	if ( buf[0] != '#' ) {
			sscanf(buf, FP_SCANSPEC, &el);
			*A++ = el;
			++i;
		}
	}

	return 0;
}

// column-major
void print_dense2mm(FILE *f, const char *name, int m, int n, const double *A, int lda) 
{
	printf("warning: writing obj %s\n", name);

	fprintf(f, "# name: %s\n", name);
	fprintf(f, "# type: matrix\n");
	fprintf(f, "# rows: %i\n", m);
	fprintf(f, "# columns: %i\n", n);

	int i;
	for ( i=0; i<m; ++i ) {
		int j;
		for ( j=0; j<n; ++j ) {
				fprintf(f, "%.16e \n", A[j*lda+i]);
//				fprintf(f, "\n");
		}
//		fprintf(f, "\n");
	}
}

void fprint_dense2mm(const char *fname, const char *name, int m, int n, const double *A, int lda) 
{
	FILE *f = fopen(fname, "w");
	if ( f == NULL ) {
		fprintf(stderr, "err: cannot open %s for writing\n", fname);
	}

	print_dense2mm(f, name, m, n, A, lda);

	fclose(f);
}

/* 
 * BLAS/LAPACK task wrappers
 * */
void __t_copy(int p, int bm, int bn, int m, int n, double *x, double *y, int initx, int inity) 
{
	double *X = &x[initx];
	double *Y = &y[inity];
	int i_one = 1;
	int j;
	for ( j=0; j<bn; ++j ) {
		BLAS_cp(bm, &X[j*m], i_one, &Y[j*m], i_one);
	}
}

void __t_dot(int p, int bm, int bn, int m, int n, double *x, double *y, int initx, int inity, double *result) 
{
	double *X = &x[initx];
	double *Y = &y[inity];
	int i_one = 1;
	double local_result[bn];
	double fp_one = 1.0;
	int j;
	for ( j=0; j<bn; ++j ) {
		local_result[j] = BLAS_dot(bm, X, i_one, Y, i_one);
		X += m;
		Y += m;
	}

	#pragma omp critical
	{
		BLAS_axpy(bn, fp_one, local_result, i_one, result, i_one);
	}
}

void __t_dot_array(int p, int bm, int bn, int m, int n, double *x, double *y, int initx, int inity, double *result, int initr)
{
	double *X = &x[initx];
	double *Y = &y[inity];
	int i_one = 1;
	double local_result[bn];
	double fp_one = 1.0;
	int j;
	for ( j=0; j<bn; ++j ) {
		result[initr+j] = BLAS_dot(bm, X, i_one, Y, i_one);
		X += m;
		Y += m;
	}
}

void _cg_dot2_array(int p, int bm, int bn, int m, int n, double *x, double *y, int initx, int inity, double *result, int initr, double *a, double *b, int inita, int initb, double *result2, int initr2) 
{
	double *X = &x[initx];
	double *Y = &y[inity];
	double *A = &a[inita];
	double *B = &b[initb];
	//double fp_one = 1.0;
	int i_one = 1;

	for ( int j=0; j<bn; ++j ) {
		result[initr+j] = BLAS_dot(bm, X, i_one, Y, i_one);
		X += m;
		Y += m;
	}

	int j;
	for ( int j=0; j<bn; ++j ) {
		result2[initr2+j] = BLAS_dot(bm, A, i_one, B, i_one);
		A += m;
		B += m;
	}
}


void _cg_dot2(int p, int bm, int bn, int m, int n, double *x, double *y, int initx, int inity, double *result, double *a, double *b, int inita, int initb, double *result2) 
{
	double *X = &x[initx];
	double *Y = &y[inity];
	double *A = &a[inita];
	double *B = &b[initb];
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

	#pragma omp critical
	{
		BLAS_axpy(bn, fp_one, local_result, i_one, result, i_one);
		BLAS_axpy(bn, fp_one, local_result2, i_one, result2, i_one);
	}
}

void __t_cpaxpy_comb(int bm, int bn, int m, int n, double alpha, double *Anum, double *Aden, double *X1, double *X2, double *Y1, double *Y2, double *Z1, double *Z2) 
{
	int i_one = 1;
	int j;
	for ( j=0; j<bn; ++j) {
		/* update of x */
		double factor = Anum[j] / Aden[j];
		BLAS_cp(bm, Y2, i_one, Z2, i_one);
		BLAS_axpy(bm, factor, X2, i_one, Z2, i_one);
		X2 += m;
		Y2 += m;
		Z2 += m;

		
		/* update of r */
		factor = alpha * factor;
		BLAS_cp(bm, Y1, i_one, Z1, i_one);
		BLAS_axpy(bm, factor, X1, i_one, Z1, i_one);
		X1 += m;
		Y1 += m;
		Z1 += m;
	}
}

void __t_extm_axpy(int bm, int bn, int m, int n, double *SAnum, double *SAden, double *X, double *Y, double *Z, int p) 
{
	int i_one = 1;
	int j;
	for ( j=0; j<bn; ++j) {
		double f = SAnum[j] / SAden[j];
		BLAS_cp(bm, &Y[j*m], i_one, &Z[j*m], i_one);
		BLAS_axpy(bm, f, &X[j*m], i_one, &Z[j*m], i_one);
	}
}

/* Non-optimal implementation of the mkl_csrmv when it is not present */
void manual_csrmv(char *trans, int m, int n, double alpha, double *avval, int *avpos, int *avptr, double *Bptr, double beta, double *Cptr)
{
	if ( strcmp(trans, "N") || strcmp(trans, "n") ) {
		for ( int i = 0; i < m; i++ ) {
			double c = Cptr[i];
			c = beta * c;
			for ( int v = avptr[i]; v < avptr[i+1]; v++ ) {
				int pos = avpos[v];
				double val = avval[v];
				c += alpha * val * Bptr[pos];
			}
			Cptr[i] = c;
		}
	} else if ( strcmp(trans, "T") || strcmp(trans, "t") ) {
		int count = 0;
		for ( int i = 0; i < n; i++, count++ ) {
			for ( int v = avptr[i]; v < avptr[i+1]; v++ ) {
				int pos = avpos[v];
				double val = avval[v];
				if ( ! count ) {
					Cptr[pos] = alpha * val * Bptr[i] + beta * Cptr[pos];
				} else {
					Cptr[pos] += alpha * val * Bptr[i];
				}
			}
		}
	}
}

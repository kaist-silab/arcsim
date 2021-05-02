/*
  Copyright ©2013 The Regents of the University of California
  (Regents). All Rights Reserved. Permission to use, copy, modify, and
  distribute this software and its documentation for educational,
  research, and not-for-profit purposes, without fee and without a
  signed licensing agreement, is hereby granted, provided that the
  above copyright notice, this paragraph and the following two
  paragraphs appear in all copies, modifications, and
  distributions. Contact The Office of Technology Licensing, UC
  Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620,
  (510) 643-7201, for commercial licensing opportunities.

  IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT,
  INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
  LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS
  DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE POSSIBILITY
  OF SUCH DAMAGE.

  REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING
  DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS
  IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
*/

#include "taucs.hpp"
#include "timer.hpp"
#include "display.hpp"
#include "alglib/solvers.h"
#include <cstdlib>
#include <iostream>
using namespace std;
using namespace alglib;

vector<Node*>* debug_nodes = 0;

extern "C" {
#include "taucs.h"
int taucs_linsolve (taucs_ccs_matrix* A, // input matrix
                    void** factorization, // an approximate inverse
                    int nrhs, // number of right-hand sides
                    void* X, // unknowns
                    void* B, // right-hand sides
                    char* options[], // options (what to do and how)
                    void* arguments[]); // option arguments
}

ostream &operator<< (ostream &out, taucs_ccs_matrix *A) {
    out << "n: " << A->n << endl;
    out << "m: " << A->m << endl;
    out << "flags: " << A->flags << endl;
    out << "colptr: ";
    for (int i = 0; i <= A->n; i++)
        out << (i==0?"":", ") << A->colptr[i];
    out << endl;
    out << "rowind: ";
    for (int j = 0; j <= A->colptr[A->n]; j++)
        out << (j==0?"":", ") << A->rowind[j];
    out << endl;
    out << "values.d: ";
    for (int j = 0; j <= A->colptr[A->n]; j++)
        out << (j==0?"":", ") << A->values.d[j];
    out << endl;
    return out;
}

taucs_ccs_matrix *sparse_to_taucs (const SpMat<double> &As) {
    // assumption: A is square and symmetric
    int n = As.n;
    int nnz = 0;
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < (int)As.rows[i].indices.size(); k++) {
            int j = As.rows[i].indices[k];
            if (j < i)
                continue;
            nnz++;
        }
    }
    taucs_ccs_matrix *At = taucs_ccs_create(n,n, nnz, TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER);
    int pos = 0;
    for (int i = 0; i < n; i++) {
        At->colptr[i] = pos;
        for (int k = 0; k < (int)As.rows[i].indices.size(); k++) {
            int j = As.rows[i].indices[k];
            if (j < i)
                continue;
            At->rowind[pos] = j;
            At->values.d[pos] = As.rows[i].entries[k];
            pos++;
        }
    }
    At->colptr[n] = pos;
    return At;
}

template <int m> taucs_ccs_matrix *sparse_to_taucs (const SpMat< Mat<m,m> > &As) {
    // assumption: A is square and symmetric
    int n = As.n;
    int nnz = 0;
    for (int i = 0; i < n; i++) {
        for (int jj = 0; jj < (int)As.rows[i].indices.size(); jj++) {
            int j = As.rows[i].indices[jj];
            if (j < i)
                continue;
            nnz += (j==i) ? m*(m+1)/2 : m*m;
        }
    }
    taucs_ccs_matrix *At = taucs_ccs_create
        (n*m,n*m, nnz, TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER);
    int pos = 0;
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < m; k++) {
            At->colptr[i*m+k] = pos;
            for (int jj = 0; jj < (int)As.rows[i].indices.size(); jj++) {
                int j = As.rows[i].indices[jj];
                if (j < i)
                    continue;
                const Mat<m,m> &Aij = As.rows[i].entries[jj];
                for (int l = (i==j) ? k : 0; l < m; l++) {
                    At->rowind[pos] = j*m+l;
                    At->values.d[pos] = Aij(k,l);
                    pos++;
                }
            }
        }
    }
    At->colptr[n*m] = pos;
    return At;
}

vector<double> alglib_linear_solve(const SpMat<double>& A, const vector<double>& b) {
    const int n = b.size();
    real_2d_array M;
    real_1d_array x,c;
    M.setlength(n,n);
    x.setlength(n);
    c.setcontent(n,&b[0]);
    for (int i = 0; i < A.m; i++) {
        const SpVec<double> &row = A.rows[i];
        for (int j = 0; j < A.n; j++)
            M(i,j) = 0;
        for (size_t jj = 0; jj < row.indices.size(); jj++) {
            int j = row.indices[jj];
            M(i,j) = row.entries[jj];
        }
    }
    
    ae_int_t info;
    densesolverreport rep;    
    rmatrixsolve(M,n,c,info,rep,x);
    
    vector<double> ret(n);
    for (int i=0; i<n; i++)
        ret[i] = x[i];
    
    return ret;
}

template<int C>
vector<Vec<C> > alglib_linear_solve_vec(const SpMat<Mat<C,C> >& A, const vector<Vec<C> >& b) {
    const int n = b.size()*C;
    real_2d_array M;
    real_1d_array x,c;
    M.setlength(n,n);
    x.setlength(n);
    c.setlength(n);
    for (int i=0; i<n; i++)
        c[i] = b[i/C][i%C];
    for (int i = 0; i < A.m; i++) {
        const SpVec<Mat<C,C> > &row = A.rows[i];
        for (size_t jj = 0; jj < row.indices.size(); jj++) {
            int j = row.indices[jj];
            for (int si=0; si<C; si++)
                for (int sj=0; sj<C; sj++)
                    M(i*C+si,j*C+sj) = row.entries[jj](si,sj);
        }
    }

    ae_int_t info;
    densesolverreport rep;    
    rmatrixsolve(M,n,c,info,rep,x);
    
    vector<Vec<C> > ret(n);
    for (int i=0; i<n; i++)
        ret[i/C][i%C] = x[i];
    
    return ret;
}

vector<double> taucs_linear_solve (const SpMat<double> &A, const vector<double> &b) {
    if (b.size() < 20) 
        return alglib_linear_solve(A,b);
    
    // taucs_logfile("stdout");
    taucs_ccs_matrix *Ataucs = sparse_to_taucs(A);
    vector<double> x(b.size());
    char *options[] = {(char*)"taucs.factor.LLT=true", NULL};
    int retval = taucs_linsolve(Ataucs, NULL, 1, &x[0], (double*)&b[0], options, NULL);
    if (retval != TAUCS_SUCCESS) {
        cerr << "Error: TAUCS(scalar) failed with return value " << retval << endl;
        segfault();
        exit(EXIT_FAILURE);
    }
    taucs_ccs_free(Ataucs);
    return x;
}

template <int m> vector< Vec<m> > taucs_linear_solve
    (const SpMat< Mat<m,m> > &A, const vector< Vec<m> > &b) {
    if (b.size() < 6) 
        return alglib_linear_solve_vec(A,b);
    
    // taucs_logfile("stdout");
    taucs_ccs_matrix *Ataucs = sparse_to_taucs(A);
    vector< Vec<m> > x(b.size());
    char *options[] = {(char*)"taucs.factor.LLT=true", NULL};
    int retval = taucs_linsolve(Ataucs, NULL, 1, &x[0], (double*)&b[0], options, NULL);
    if (retval != TAUCS_SUCCESS) {
        cerr << "Error: TAUCS(vector) failed with return value " << retval << endl;
        segfault();
        exit(EXIT_FAILURE);
    }
    taucs_ccs_free(Ataucs);
    return x;
}

template vector<Vec3> taucs_linear_solve (const SpMat<Mat3x3> &A,
                                          const vector<Vec3> &b);

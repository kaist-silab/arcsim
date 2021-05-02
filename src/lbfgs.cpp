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

#include "optimization.hpp"

#include "alglib/optimization.h"
#include <omp.h>
#include <vector>
using namespace std;
using namespace alglib;

static const NLOpt *problem;

static void lbfgs_value_and_grad (const real_1d_array &x, double &value,
                                  real_1d_array &grad, void *ptr=NULL);

void l_bfgs_method (const NLOpt &problem, OptOptions opt, bool verbose) {
    ::problem = &problem;
    real_1d_array x;
    x.setlength(::problem->nvar);
    ::problem->initialize(&x[0]);
    minlbfgsstate state;
    minlbfgsreport rep;
    int m = 10;
    minlbfgscreate(m, x, state);
    minlbfgssetcond(state, opt.eps_g(), opt.eps_f(), opt.eps_x(),
                           opt.max_iter());
    minlbfgsoptimize(state, lbfgs_value_and_grad);
    minlbfgsresults(state, x, rep);
    if (verbose)
        cout << rep.iterationscount << " iterations" << endl;
    ::problem->finalize(&x[0]);
}

/*static void add (real_1d_array &x, const vector<double> &y) {
    for (int i = 0; i < (int)y.size(); i++)
        x[i] += y[i];
}*/

static void lbfgs_value_and_grad (const real_1d_array &x, double &value,
                                 real_1d_array &grad, void *ptr) {
    ::problem->precompute(&x[0]);
    value = ::problem->objective(&x[0]);
    ::problem->gradient(&x[0], &grad[0]);
}

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

#include "physics.hpp"

#include "blockvectors.hpp"
#include "collisionutil.hpp"
#include "sparse.hpp"
#include "taucs.hpp"
#include "io.hpp"
#include "display.hpp"

using namespace std;

static const bool verbose = false;
extern bool consistency_check;
extern vector<Node*>* debug_nodes;

static void consistency(vector<Vec3>& b, const string& name) {
	if (consistency_check) {
		cout << "> physics: " << name << " ";
		test_state(b, "/tmp/b");
	}
}

typedef Mat<9,9> Mat9x9;
typedef Mat<9,6> Mat9x6;
typedef Mat<6,6> Mat6x6;
typedef Mat<4,6> Mat4x6;
typedef Mat<3,4> Mat3x4;
typedef Mat<4,9> Mat4x9;
typedef Vec<9> Vec9;

// A kronecker B = [a11 B, a12 B, ..., a1n B;
//                  a21 B, a22 B, ..., a2n B;
//                   ... ,  ... , ...,  ... ;
//                  am1 B, am2 B, ..., amn B]
template <int m, int n, int p, int q>
Mat<m*p,n*q> kronecker (const Mat<m,n> &A, const Mat<p,q> &B) {
    Mat<m*p,n*q> C;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < p; k++)
                for (int l = 0; l < q; l++)
                    C(i*p+k,j*q+l) = A(i,j)*B(k,l);
    return C;
}

template <int m> Mat<m,1> colmat (const Vec<m> &v) {
    Mat<1,m> A; for (int i = 0; i < m; i++) A(i,0) = v[i]; return A;}
template <int n> Mat<1,n> rowmat (const Vec<n> &v) {
    Mat<1,n> A; for (int i = 0; i < n; i++) A(0,i) = v[i]; return A;}

template <Space s>
Mat3x3 deformation_gradient (const Face *face) {
    return derivative(pos<s>(face->v[0]->node), pos<s>(face->v[1]->node), 
                      pos<s>(face->v[2]->node), normal<s>(face), face) * face->Sp_str;
}

Mat3x3 material_model (const Face *face, const Mat3x3& G) {
    const Material* mat = face->material;
    double weakening_mult = 1/(1 + mat->weakening * face->damage);
    if (mat->use_dde) {
        Vec4 k = stretching_stiffness(reduce_xy(G), mat->dde_stretching) * weakening_mult;
        Mat3x3 sigma (Vec3(k[0]*G(0,0)+k[1]*G(1,1), 0.5*k[3]*G(0,1), 0),
                      Vec3(0.5*k[3]*G(1,0), k[2]*G(1,1)+k[1]*G(0,0), 0),
                      Vec3(0, 0, 0));
        return sigma;
    } else {
    	double A = mat->alt_stretching * weakening_mult;
        Mat3x3 sigma = A*(1.0-mat->alt_poisson)*G + Mat3x3(A*mat->alt_poisson*trace(G));
        return sigma;
    }
}

template <Space s>
double stretching_energy (const Face *face) {
    Mat3x3 F = deformation_gradient<WS>(face);
    Mat3x3 G = (F.t()*F - Mat3x3(1)) * 0.5;
    Mat3x3 sigma = material_model(face, G);

    return face->a * 0.5 * inner(sigma, G);
}

template <Space s>
pair<Mat9x9,Vec9> stretching_force (const Face *face) {
    const Material* mat = face->material;
    // compute stress, strain
    const Vec3 x[3] = { pos<s>(face->v[0]->node), pos<s>(face->v[1]->node), pos<s>(face->v[2]->node) };
    Mat3x3 F = deformation_gradient<s>(face);
    Mat3x3 G = (F.t()*F - Mat3x3(1)) * 0.5;
    double weakening_mult = 1/(1 + mat->weakening * face->damage);

    Mat3x3 Y = face->invDm * face->Sp_str;
    Mat3x3 D = Mat3x3::rows(-Y.row(0)-Y.row(1),Y.row(0),Y.row(1));
    Mat<3,9> DD[3] = { kronecker(rowmat(D.col(0)), Mat3x3(1)),
                       kronecker(rowmat(D.col(1)), Mat3x3(1)),
                       kronecker(rowmat(D.col(2)), Mat3x3(1)) };
    Vec9 X;
    for (int i=0; i<9; i++) X[i] = x[i/3][i%3];
    Vec3 f[3] = {DD[0]*X,DD[1]*X,DD[2]*X};
    
    Vec9 grad_f(0);
    Mat9x9 hess_f(0);
    if (mat->use_dde) {
        Vec4 k = stretching_stiffness(reduce_xy(G), mat->dde_stretching) * weakening_mult;

        const Mat<3,9>& Du = DD[0];
        const Mat<3,9>& Dv = DD[1];
        Vec9 fuu = Du.t()*f[0], fvv = Dv.t()*f[1], fuv = (Du.t()*f[1] + Dv.t()*f[0])/2.;
        grad_f = k[0]*G(0,0)*fuu + k[2]*G(1,1)*fvv
               + k[1]*(G(0,0)*fvv + G(1,1)*fuu) + 2*k[3]*G(0,1)*fuv;
        hess_f = k[0]*(outer(fuu,fuu) + max(G(0,0),0.)*Du.t()*Du)
               + k[2]*(outer(fvv,fvv) + max(G(1,1),0.)*Dv.t()*Dv)
               + k[1]*(outer(fuu,fvv) + max(G(0,0),0.)*Dv.t()*Dv
                     + outer(fvv,fuu) + max(G(1,1),0.)*Du.t()*Du)
               + 2.*k[3]*(outer(fuv,fuv));
        // ignoring G(0,1)*(Du.t()*Dv+Dv.t()*Du)/2. term
        // because may not be positive definite
        return make_pair(-face->a*hess_f, -face->a*grad_f);
    } else {
        double gf = -face->a * mat->alt_stretching;
        //Vec9 d_trace = 0.5 * (DD[0].t()*DD[0]+DD[1].t()*DD[1]+DD[2].t()*DD[2]) * X;
        Mat3x3 Gc = max(G,Mat3x3(0)); // posdef
        for (int i=0; i<3; i++)
            for(int j=0; j<=i; j++) {
                Vec9 dG = 0.5 * (DD[i].t()*f[j] + DD[j].t()*f[i]);
                Mat9x9 dG2 = 0.5 * (DD[i].t()*DD[j] + DD[j].t()*DD[i]);
                Vec9 d_trace_c = 0.5*dG; // posdef

                if (i==j)
                    grad_f += gf * (1.0 - mat->alt_poisson) * G(i,j) * dG +
                              gf * mat->alt_poisson * trace(G) * dG;
                else
                    grad_f += 2.0 * gf * (1.0 - mat->alt_poisson) * G(i,j) * dG;
                
                if (i==j)
                    hess_f += gf * (1.0 - mat->alt_poisson) * (outer(dG,dG) + Gc(i,j) * dG2) +
                              gf * mat->alt_poisson * (trace(Gc) * dG2 + outer(d_trace_c, dG));
                else
                    hess_f += 2.0 * gf * (1.0 - mat->alt_poisson) * (outer(dG,dG) +
                              Mat9x9(0)); // posdef
            }
        return make_pair(hess_f, grad_f);
    }
}

typedef Mat<12,12> Mat12x12;
typedef Vec<12> Vec12;

double bending_coeff(const Edge* edge, double theta) {
    const Face *face0 = edge->adjf[0], *face1 = edge->adjf[1];
    double a = face0->a + face1->a;
    double l = norm(edge->n[1]->x - edge->n[0]->x);
    
    double ke0 = (face0->material->use_dde) ? 
    	bending_stiffness(edge, 0, face0->material->dde_bending, l, theta) : face0->material->alt_bending;
    double ke1 = (face1->material->use_dde) ? 
    	bending_stiffness(edge, 1, face1->material->dde_bending, l, theta) : face1->material->alt_bending;
    
    double ke = min(ke0, ke1);
    double weakening = max(face0->material->weakening, face1->material->weakening);
    ke *= 1/(1 + weakening*edge->damage);
    double shape = sq(l) / (2*a);

    return ke * shape;
}

template <Space s>
double bending_energy (const Edge *edge) {
    if (!edge->adjf[0] || !edge->adjf[1])
        return 0;

    double theta = dihedral_angle<s>(edge);
    return bending_coeff(edge, theta) * sq(theta - edge->theta_ideal)/4;
}

double distance (const Vec3 &x, const Vec3 &a, const Vec3 &b) {
    Vec3 e = b-a;
    Vec3 xp = e*dot(e, x-a)/dot(e,e);
    // return norm((x-a)-xp);
    return max(norm((x-a)-xp), 1e-3*norm(e));
}

Vec2 barycentric_weights (const Vec3 &x, const Vec3 &a, const Vec3 &b) {
    Vec3 e = b-a;
    double t = dot(e, x-a)/dot(e,e);
    return Vec2(1-t, t);
}

template <Space s>
pair<Mat12x12,Vec12> bending_force (const Edge *edge) {
    const Face *face0 = edge->adjf[0], *face1 = edge->adjf[1];
    if (!face0 || !face1)
        return make_pair(Mat12x12(0), Vec12(0));
    double theta = dihedral_angle<s>(edge);
    Vec3 x0 = pos<s>(edge->n[0]),
         x1 = pos<s>(edge->n[1]),
         x2 = pos<s>(edge_opp_vert(edge, 0)->node),
         x3 = pos<s>(edge_opp_vert(edge, 1)->node);
    double h0 = distance(x2, x0, x1), h1 = distance(x3, x0, x1);
    Vec3 n0 = normal<s>(face0), n1 = normal<s>(face1);
    Vec2 w_f0 = barycentric_weights(x2, x0, x1),
         w_f1 = barycentric_weights(x3, x0, x1);
    Vec12 dtheta = mat_to_vec(Mat3x4(-(w_f0[0]*n0/h0 + w_f1[0]*n1/h1),
                                     -(w_f0[1]*n0/h0 + w_f1[1]*n1/h1),
                                     n0/h0,
                                     n1/h1));

    double coeff = bending_coeff(edge, theta);
    return make_pair(-coeff * outer(dtheta, dtheta)/2.,
                     -coeff * (theta - edge->theta_ideal)*dtheta/2.);
}

template <int m, int n> Mat<3,3> submat3 (const Mat<m,n> &A, int i, int j) {
    Mat3x3 Asub;
    for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
            Asub(k,l) = A(i*3+k, j*3+l);
    return Asub;
}

template <int n> Vec<3> subvec3 (const Vec<n> &b, int i) {
    Vec3 bsub;
    for (int k = 0; k < 3; k++)
        bsub[k] = b[i*3+k];
    return bsub;
}

template <int m> void add_submat (const Mat<m*3,m*3> &Asub, const Vec<m,int> &ix, SpMat<Mat3x3> &A) {
    for (int i = 0; i < m; i++) {
        if (ix[i] < 0) continue;
        for (int j = 0; j < m; j++) {
            if (ix[j] < 0) continue;
            A(ix[i],ix[j]) += submat3(Asub, i,j);
        }
    }
}

template <int m> void add_subvec (const Vec<m*3> &bsub, const Vec<m,int> &ix, vector<Vec3> &b) {
    for (int i = 0; i < m; i++) {
        if (ix[i] < 0) continue;
        b[ix[i]] += subvec3(bsub, i);
    }
}

Vec<3,int> indices (const Node *n0, const Node *n1, const Node *n2) {
    Vec<3,int> ix;
    ix[0] = n0->active() ? n0->index : -1;
    ix[1] = n1->active() ? n1->index : -1;
    ix[2] = n2->active() ? n2->index : -1;
    return ix;
}

Vec<4,int> indices (const Node *n0, const Node *n1,
                    const Node *n2, const Node *n3) {
    Vec<4,int> ix;
    ix[0] = n0->active() ? n0->index : -1;
    ix[1] = n1->active() ? n1->index : -1;
    ix[2] = n2->active() ? n2->index : -1;
    ix[3] = n3->active() ? n3->index : -1;
    return ix;
}

template <Space s>
double internal_energy (const vector<Face*>& faces, const vector<Edge*>& edges) {
    double E = 0;
    for (size_t f = 0; f < faces.size(); f++)
        E += stretching_energy<s>(faces[f]);
    for (size_t e = 0; e < edges.size(); e++) {
        E += bending_energy<s>(edges[e]);
    }
    return E;
}
template double internal_energy<PS> (const vector<Face*>&, const vector<Edge*>&);
template double internal_energy<WS> (const vector<Face*>&, const vector<Edge*>&);

// A = dt^2 J + dt damp J
// b = dt f + dt^2 J v + dt damp J v

template <Space s>
void add_internal_forces (const vector<Face*>& faces, const vector<Edge*>& edges,
						  SpMat<Mat3x3> &A, vector<Vec3> &b, double dt) {
    
    for (size_t f = 0; f < faces.size(); f++) {
        const Face* face = faces[f];
        const Node *n0 = face->v[0]->node, *n1 = face->v[1]->node,
                   *n2 = face->v[2]->node;
        Vec9 vs = mat_to_vec(Mat3x3(n0->v, n1->v, n2->v));
        pair<Mat9x9,Vec9> membF = stretching_force<s>(face);
        Mat9x9 J = membF.first;
        Vec9 F = membF.second;
        if (dt == 0) {
            add_submat(-J, indices(n0,n1,n2), A);
            add_subvec(F, indices(n0,n1,n2), b);
        } else {
            double damping = face->material->damping;
            add_submat(-dt*(dt+damping)*J, indices(n0,n1,n2), A);
            add_subvec(dt*(F + (dt+damping)*J*vs), indices(n0,n1,n2), b);
        }
    }
    for (size_t e = 0; e < edges.size(); e++) {
        const Edge *edge = edges[e];
        if (!edge->adjf[0] || !edge->adjf[1])
            continue;
        pair<Mat12x12,Vec12> bendF = bending_force<s>(edge);
        const Node *n0 = edge->n[0],
                   *n1 = edge->n[1],
                   *n2 = edge_opp_vert(edge, 0)->node,
                   *n3 = edge_opp_vert(edge, 1)->node;
        Vec12 vs = mat_to_vec(Mat3x4(n0->v, n1->v, n2->v, n3->v));
        Mat12x12 J = bendF.first;
        Vec12 F = bendF.second;
        
        if (dt == 0) {
            add_submat(-J, indices(n0,n1,n2,n3), A);
            add_subvec(F, indices(n0,n1,n2,n3), b);
        } else {
            double damping = (edge->adjf[0]->material->damping +
							  edge->adjf[1]->material->damping) * 0.5;
            add_submat(-dt*(dt+damping)*J, indices(n0,n1,n2,n3), A);
            add_subvec(dt*(F + (dt+damping)*J*vs), indices(n0,n1,n2,n3), b);
        }
    }
}
template void add_internal_forces<PS> (const vector<Face*>&, const vector<Edge*>&, 
                                       SpMat<Mat3x3> &, vector<Vec3>&, double);
template void add_internal_forces<WS> (const vector<Face*>&, const vector<Edge*>&, 
                                       SpMat<Mat3x3> &, vector<Vec3>&, double);

double constraint_energy (const vector<Constraint*> &cons) {
    double E = 0;
    for (int c = 0; c < (int)cons.size(); c++) {
        double value = cons[c]->value();
        double e = cons[c]->energy(value);
        E += e;
    }
    return E;
}

void add_constraint_forces (const vector<Constraint*> &cons,
                            SpMat<Mat3x3> &A, vector<Vec3> &b, double dt) {
    for (int c = 0; c < (int)cons.size(); c++) {
        double value = cons[c]->value();
        double g = cons[c]->energy_grad(value);
        double h = cons[c]->energy_hess(value);
        MeshGrad grad = cons[c]->gradient();
        // f = -g*grad
        // J = -h*outer(grad,grad)
        double v_dot_grad = 0;
        for (MeshGrad::iterator it = grad.begin(); it != grad.end(); it++) {
            v_dot_grad += dot(it->f, it->node->v);
        }
        for (MeshGrad::iterator it = grad.begin(); it != grad.end(); it++) {
            const Node *nodei = it->node;
            if (!nodei->active())
                continue;
            int ni = nodei->index;
            for (MeshGrad::iterator jt = grad.begin(); jt != grad.end(); jt++) {
                const Node *nodej = jt->node;
                if (!nodej->active())
                    continue;
                int nj = nodej->index;
                if (dt == 0)
                    A(ni,nj) += h*outer(it->f, jt->f);
                else
                    A(ni,nj) += dt*dt*h*outer(it->f, jt->f);
            }
            if (dt == 0)
                b[ni] -= g*it->f;
            else
                b[ni] -= dt*(g + dt*h*v_dot_grad)*it->f;
        }
    }
}

void add_friction_forces (const vector<Constraint*> cons,
                          SpMat<Mat3x3> &A, vector<Vec3> &b, double dt) {
    for (int c = 0; c < (int)cons.size(); c++) {
        MeshHess jac;
        MeshGrad force = cons[c]->friction(dt, jac);
        for (MeshGrad::iterator it = force.begin(); it != force.end(); it++) {
            const Node *node = it->node;
            if (!node->active())
                continue;
            b[node->index] += dt*it->f;
        }
        for (MeshHess::iterator it = jac.begin(); it != jac.end(); it++) {
            const Node *nodei = it->i, *nodej = it->j;
            if (!nodei->active() || !nodej->active())
                continue;
            A(nodei->index, nodej->index) -= dt*it->J;
        }
    }
}

vector<Vec3> implicit_update (vector<Node*>& nodes, const vector<Edge*>& edges, const vector<Face*>& faces,
					  const vector<Vec3> &fext, const vector<Mat3x3> &Jext,
                      const vector<Constraint*> &cons, double dt) {
    vector<Vert*>::iterator vert_it;
    vector<Face*>::iterator face_it;
    int nn = nodes.size();

    // M Dv/Dt = F (x + Dx) = F (x + Dt (v + Dv))
    // Dv = Dt (M - Dt2 F)i F (x + Dt v)
    // A = M - Dt2 F
    // b = Dt F (x + Dt v)
    SpMat<Mat3x3> A(nn,nn);
    vector<Vec3> b(nn, Vec3(0));
    for (size_t n = 0; n < nodes.size(); n++) {
        A(n,n) += Mat3x3(nodes[n]->m) - dt*dt*Jext[n];
        b[n] += dt*fext[n];        
    }
    consistency((vector<Vec3>&)fext, "fext");
    consistency(b, "init");
    add_internal_forces<WS>(faces, edges, A, b, dt);
    consistency(b, "internal forces");
    add_constraint_forces(cons, A, b, dt);
    consistency(b, "constraints");
    add_friction_forces(cons, A, b, dt);
    consistency(b, "friction");
    ::debug_nodes = &nodes;
    vector<Vec3> dv = taucs_linear_solve(A, b);
    ::debug_nodes = 0;
    consistency(dv, "taucs");
    return dv;    
}

Vec3 wind_force (const Face *face, const Wind &wind) {
    Vec3 vface = (face->v[0]->node->v + face->v[1]->node->v
                  + face->v[2]->node->v)/3.;
    Vec3 vrel = wind.velocity - vface;
    double vn = dot(face->n, vrel);
    Vec3 vt = vrel - vn*face->n;
    return wind.density*face->a*abs(vn)*vn*face->n + wind.drag*face->a*vt;
}

void add_external_forces (const vector<Node*>& nodes, const vector<Face*>& faces, const Vec3 &gravity,
                          const Wind &wind, vector<Vec3> &fext,
                          vector<Mat3x3> &Jext) {
    for (size_t n = 0; n < nodes.size(); n++) {
        fext[n] += nodes[n]->m*gravity;
    }
    for (size_t f = 0; f < faces.size(); f++) {
        const Face *face = faces[f];
        Vec3 fw = wind_force(face, wind);
        if (norm2(fw)>0)
            for (int v = 0; v < 3; v++)
                fext[face->v[v]->node->index] += fw/3.;
    }
}

void add_morph_forces (const Cloth &cloth, const Morph &morph, double t,
                       double dt, vector<Vec3> &fext, vector<Mat3x3> &Jext) {
    const Mesh &mesh = cloth.mesh;
    for (int v = 0; v < (int)mesh.verts.size(); v++) {
        const Vert *vert = mesh.verts[v];
        Vec3 x = morph.pos(t, vert->u);
        double stiffness = exp(morph.log_stiffness.pos(t));
        //Vec3 n = vert->node->n;
        double s = stiffness*vert->node->a;
        // // lower stiffness in tangential direction
        // Mat3x3 k = s*outer(n,n) + (s/10)*(Mat3x3(1) - outer(n,n));
        Mat3x3 k = Mat3x3(s);
        double c = sqrt(s*vert->node->m); // subcritical damping
        Mat3x3 d = c/s * k;
        fext[vert->node->index] -= k*(vert->node->x - x);
        fext[vert->node->index] -= d*vert->node->v;
        Jext[vert->node->index] -= k + d/dt;
    }
}

void project_outside (vector<Node*>& nodes, const vector<Constraint*> &cons) {
    int nn = nodes.size();
    vector<double> w(nn, 0);
    vector<Vec3> dx(nn, Vec3(0));
    for (int c = 0; c < (int)cons.size(); c++) {
        MeshGrad dxc = cons[c]->project();
        for (MeshGrad::iterator it = dxc.begin(); it != dxc.end(); it++) {
            const Node *node = it->node;
            if (!node->active())
                continue;
            double wn = norm2(it->f);
            int n = node->index;
            w[n] += wn;
            dx[n] += wn*it->f;
        }
    }
    for (int n = 0; n < nn; n++) {
        if (w[n] == 0)
            continue;
        nodes[n]->x += dx[n]/w[n];
    }
}

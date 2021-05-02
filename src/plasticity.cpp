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

#include "plasticity.hpp"

#include "geometry.hpp"
#include "optimization.hpp"
#include "physics.hpp"
#include <omp.h>

using namespace std;

static const double mu = 1e-6;

Vec3 face_to_edges (const Mat3x3 &S, const Face *face);

void reset_plasticity (Cloth &cloth) {
    Mesh &mesh = cloth.mesh;
    for (int n = 0; n < (int)mesh.nodes.size(); n++)
        mesh.nodes[n]->y = mesh.nodes[n]->x;
    for (int e = 0; e < (int)mesh.edges.size(); e++) {
        Edge *edge = mesh.edges[e];
        edge->theta_ideal = dihedral_angle<WS>(edge);
        edge->damage = 0;
    }
    for (int f = 0; f < (int)mesh.faces.size(); f++) {
        recompute_Sp_bend(mesh.faces[f]);
        mesh.faces[f]->damage = 0;
    }
}

void recompute_edge_plasticity (Mesh &mesh);

void optimize_plastic_embedding (Cloth &cloth);

void plastic_update (Cloth &cloth) {
    Mesh &mesh = cloth.mesh;
    for (size_t f = 0; f < mesh.faces.size(); f++) {
        Face *face = mesh.faces[f];
        double S_yield = face->material->yield_curv;
        Mat3x3 S_total = curvature<WS>(face);
        Mat3x3 S_elastic = S_total - face->Sp_bend;
        double dS = norm_F(S_elastic);
        if (dS > S_yield) {
            face->Sp_bend += S_elastic/dS*(dS - S_yield);
            face->damage += dS/S_yield - 1;
        }
    }
    recompute_edge_plasticity(cloth.mesh);
}

// ------------------------------------------------------------------ //

struct EmbedOpt: public NLOpt {
    Cloth &cloth;
    Mesh &mesh;
    vector<Vec3> y0;
    mutable vector<Vec3> f;
    mutable SpMat<Mat3x3> J;
    EmbedOpt (Cloth &cloth): cloth(cloth), mesh(cloth.mesh) {
        int nn = mesh.nodes.size();
        nvar = nn*3;
        y0.resize(nn);
        for (int n = 0; n < nn; n++)
            y0[n] = mesh.nodes[n]->y;
        // f.resize(nv);
        // J = SpMat<Mat3x3>(nv,nv);
    }
    void initialize (double *x) const;
    void precompute (const double *x) const;
    double objective (const double *x) const;
    void gradient (const double *x, double *g) const;
    bool hessian (const double *x, SpMat<double> &H) const;
    void finalize (const double *x) const;
};

void reduce_stretching_stiffnesses (vector<Material*> &materials);
void restore_stretching_stiffnesses (vector<Material*> &materials);

void optimize_plastic_embedding (Cloth &cloth) {
    // vector<Material> materials = cloth.materials;
    reduce_stretching_stiffnesses(cloth.materials);
    line_search_newtons_method(EmbedOpt(cloth), OptOptions().max_iter(1));
    restore_stretching_stiffnesses(cloth.materials);
    // cloth.materials = materials;
}

void EmbedOpt::initialize (double *x) const {
    for (int n = 0; n < (int)mesh.nodes.size(); n++)
        set_subvec(x, n, Vec3(0));
}

void EmbedOpt::precompute (const double *x) const {
    int nn = mesh.nodes.size();
    f.assign(nn, Vec3(0));
    J = SpMat<Mat3x3>(nn,nn);
    for (int n = 0; n < (int)mesh.nodes.size(); n++) {
        mesh.nodes[n]->y = y0[n] + get_subvec(x, n);
        mesh.nodes[n]->index = n;
    }
    add_internal_forces<PS>(cloth.mesh.faces, cloth.mesh.edges, J, f, 0);
}

double EmbedOpt::objective (const double *x) const {
    for (int n = 0; n < (int)mesh.nodes.size(); n++)
        mesh.nodes[n]->y = y0[n] + get_subvec(x, n);
    return internal_energy<PS>(cloth.mesh.faces, cloth.mesh.edges);
}

void EmbedOpt::gradient (const double *x, double *g) const {
    for (int n = 0; n < (int)mesh.nodes.size(); n++) {
        //const Node *node = mesh.nodes[n];
        set_subvec(g, n, -f[n]);
    }
}

bool EmbedOpt::hessian (const double *x, SpMat<double> &H) const {
    for (int i = 0; i < (int)mesh.nodes.size(); i++) {
        const SpVec<Mat3x3> &Ji = J.rows[i];
        for (int jj = 0; jj < (int)Ji.indices.size(); jj++) {
            int j = Ji.indices[jj];
            const Mat3x3 &Jij = Ji.entries[jj];
            set_submat(H, i, j, Jij);
        }
        add_submat(H, i, i, Mat3x3(::mu));
    }
    return true;
}

void EmbedOpt::finalize (const double *x) const {
    for (int n = 0; n < (int)mesh.nodes.size(); n++)
        mesh.nodes[n]->y = y0[n] + get_subvec(x, n);
}

void reduce_stretching_stiffnesses (vector<Material*> &materials) {
    for (int m = 0; m < (int)materials.size(); m++)
        for (int i = 0; i < 40; i++)
            for (int j = 0; j < 40; j++)
                for (int k = 0; k < 40; k++)
                    materials[m]->dde_stretching.s[i][j][k] *= 1e-2;
}

void restore_stretching_stiffnesses (vector<Material*> &materials) {
    for (int m = 0; m < (int)materials.size(); m++)
        for (int i = 0; i < 40; i++)
            for (int j = 0; j < 40; j++)
                for (int k = 0; k < 40; k++)
                    materials[m]->dde_stretching.s[i][j][k] *= 1e2;
}

// ------------------------------------------------------------------ //

Mat3x3 edges_to_face (const Vec3 &theta, const Face *face) {
    Mat3x3 S;
    Vec3 n = normal<MS>(face);    
    for (int e = 0; e < 3; e++) {
        //const Edge *edge = face->adje[e];
        Vec3 e_mat = face->v[PREV(e)]->u - face->v[NEXT(e)]->u,
             t_mat = cross(normalize(e_mat),n);
        S -= 1/2.*theta[e]*norm(e_mat)*outer(t_mat, t_mat);
    }
    S /= face->a;
    return S;
}

Vec3 face_to_edges (const Mat3x3 &S, const Face *face) {
	const Vec3 n = normal<MS>(face);
    Mat<6,3> A;
    
    for (int e = 0; e < 3; e++) {
        Vec3 e_mat = face->v[PREV(e)]->u - face->v[NEXT(e)]->u,
             t_mat = cross(normalize(e_mat),n);
        Mat3x3 Se = -1/2.*norm(e_mat)*outer(t_mat, t_mat);
        A.col(e) = Vec<6>(Se(0,0),Se(1,1),Se(2,2),Se(0,1),Se(0,2),Se(1,2));
    }

    Vec<6> y = face->a * Vec<6>(S(0,0),S(1,1),S(2,2),S(0,1),S(0,2),S(1,2));
    return solve_llsq(A,y);
}

void recompute_edge_plasticity (Mesh &mesh) {
    for (int e = 0; e < (int)mesh.edges.size(); e++) {
        mesh.edges[e]->theta_ideal = 0;
        mesh.edges[e]->damage = 0;
    }
    for (int f = 0; f < (int)mesh.faces.size(); f++) {
        const Face *face = mesh.faces[f];
        Vec3 theta = face_to_edges(face->Sp_bend, face);
        for (int e = 0; e < 3; e++) {
            face->adje[e]->theta_ideal += theta[e];
            face->adje[e]->damage += face->damage;
        }
    }
    for (int e = 0; e < (int)mesh.edges.size(); e++) {
        Edge *edge = mesh.edges[e];
        if (edge->adjf[0] && edge->adjf[1]) {// edge has two adjacent faces
            edge->theta_ideal /= 2;
            edge->damage /= 2;
        }
    }
}

// ------------------------------------------------------------------ //

void recompute_Sp_bend (Face *face) {
    Vec3 theta;
    for (int e = 0; e < 3; e++)
        theta[e] = face->adje[e]->theta_ideal;
    face->Sp_bend = edges_to_face(theta, face); // add residual and PS reconstruction
}

// ------------------------------------------------------------------ //

Mat3x3 stretch_plasticity_from_embedding(const Face *face) {
    double limit = face->material->plastic_limit;
    Mat3x3 F_ps = derivative(pos<PS>(face->v[0]->node),
                             pos<PS>(face->v[1]->node),
                             pos<PS>(face->v[2]->node), normal<PS>(face), face);
    SVD<3,3> svd = singular_value_decomposition(F_ps);
    for (int j=0;j<3; j++)
        svd.s[j] = clamp(1.0/svd.s[j],1./limit,limit);
    Mat3x3 M = svd.Vt.t() * diag(svd.s) * svd.Vt;
    return M;
}

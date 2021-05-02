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

#include "separate.hpp"

#include <omp.h>
#include <set>

#include "collisionutil.hpp"
#include "geometry.hpp"
#include "io.hpp"
#include "magic.hpp"
#include "optimization.hpp"
#include "simulation.hpp"
#include "util.hpp"
using namespace std;

static const double &thickness = ::magic.projection_thickness;

static double obs_mass;
static bool deform_obstacles;
static double get_mass (const Node *node) {
    return is_free(node) ? node->m : obs_mass;}

typedef Vec3 Bary; // barycentric coordinates
typedef std::pair<Vec3, Vec3> Line3;

struct Ixn {  // intersection
    Face *f0, *f1;
    double l;
    Vec3 g0[3], g1[3]; // dl/dx for each of the faces' nodes
    Ixn () {}
    Ixn (const Face *f0, const Face *f1)
        : f0((Face*)f0), f1((Face*)f1) {}
};

void update_active (const vector<AccelStruct*> &accs,
                    const vector<AccelStruct*> &obs_accs,
                    const vector<Ixn> &ixns);

vector<Ixn> find_intersections (const vector<AccelStruct*> &accs,
        const vector<AccelStruct*> &obs_accs);
vector<Ixn> find_overlappings (const vector<AccelStruct*> &accs,
        const vector<AccelStruct*> &obs_accs);

bool edge_face_intersection (const Edge* edge, const Face *face, Vec3& pt);
bool edge_face_intersection (const Vec3& e0, const Vec3& e1,
                const Vec3& f0, const Vec3& f1, const Vec3& f2, const Vec3& fn, Vec3& pt);
bool get_line_of_intersection (const Face *face0, const Face *face1, Line3& line);

void compute_length_and_gradient (Ixn &ixn);

vector< vector<Node*> > connected_components (const vector<Ixn> &ixns);

// #include "display.hpp"

void separate (vector<Mesh*> &meshes, const vector<Mesh*> &obs_meshes) {
    if (false) {
        // Unit tests.
        Vec3 e0, e1;
        Vec3 f0, f1, f2, fn;
        Vec3 pt;

        f0 = Vec3(1, 1, 0); f1 = Vec3(-1, 1, 0); f2 = Vec3(0, -1, 0);
        fn = Vec3(0, 0, 1);

        e0 = Vec3(0, 0.5, -0.5); e1 = Vec3(0, 0.5, 0.5);
        assert(edge_face_intersection(e0, e1, f0, f1, f2, fn, pt));
        e0 = Vec3(1, 1, -1.0); e1 = Vec3(1, 1, 1.0);
        assert(edge_face_intersection(e0, e1, f0, f1, f2, fn, pt));
        e0 = Vec3(0, 0.5, 0.1); e1 = Vec3(0, 0.5, 1.0);
        assert(!edge_face_intersection(e0, e1, f0, f1, f2, fn, pt));
        e0 = Vec3(1, 0.5, -0.5); e1 = Vec3(1, 0.5, 0.5);
        assert(!edge_face_intersection(e0, e1, f0, f1, f2, fn, pt));
        // End of unit tests.
    }

    vector<AccelStruct*> accs = create_accel_structs(meshes, false),
                         obs_accs = create_accel_structs(obs_meshes, false);

    // double D = ::magic.separation_step_size;  // step size
    // //std::cout << "Separation with step_size: " << D << std::endl;

    std::map<Node*, Vec3> g;

    ::obs_mass = 1e3;
    for (int deform = 0; deform <= 1; deform++) {
        ::deform_obstacles = deform;

        for (size_t iter = 0; iter < 100; ++iter) {
            const bool print_logs = true; // (iter % 16 == 0);

            vector<Ixn> ixns = find_intersections(accs, obs_accs);
            if (ixns.empty())
                break;

            // Aggregate the gradient vectors.
            double l_total = 0;
            g.clear();
            for (size_t i = 0; i < ixns.size(); i++) {
                Ixn &ixn = ixns[i];
                l_total += ixn.l;
                for (int v = 0; v < 3; v++) {
                    g[ixn.f0->v[v]->node] += ixn.g0[v];
                    g[ixn.f1->v[v]->node] += ixn.g1[v];
                }
            }

            if (l_total == 0)
                break;

            // for (size_t i = 0; i < ixns.size(); i++) {
            //     if (ixns[i].l == 0)
            //         continue;
            //     Annotation::add(ixns[i].f0);
            //     Annotation::add(ixns[i].f1);
            // }
            // wait_key();

            if (print_logs)
                std::cout << "iter: " << iter
                          << " n_ixns: " << ixns.size()
                          << " l_total: " << l_total << std::endl;

            vector< vector<Node*> > components = connected_components(ixns);

            // for (size_t c = 0; c < components.size(); c++) {
            //     double t = 2*M_PI*c/components.size();
            //     Vec3 color = (Vec3(cos(t),cos(t+2*M_PI/3),cos(t+4*M_PI/3))
            //                   + Vec3(1))/2.;
            //     for (size_t i = 0; i < components[c].size(); i++)
            //         Annotation::add(components[c][i], color);
            // }
            // wait_key();

            // Find per-component displacement x = C y
            // which minimizes 1/2 x' M x = 1/2 y' Mc y
            // such that -l >= g' x = gc' y
            // So y = -l (Mc^-1 gc) / (gc' Mc^-1 gc)

            // Compute inverse-mass-weighted gradient
            vector<double> mc(components.size()); // mass of component
            vector<Vec3> gc(components.size()); // per-component gradient
            double denom = 0;
            for (size_t c = 0; c < components.size(); c++) {
                const vector<Node*> &component = components[c];
                mc[c] = 0;
                gc[c] = Vec3(0);
                for (size_t n = 0; n < component.size(); n++) {
                    Node *node = component[n];
                    if (!is_free(node) && !deform_obstacles)
                        continue;
                    mc[c] += get_mass(node);
                    gc[c] += g[node];
                }
                denom += norm2(gc[c])/mc[c];
            }
            double step_length = (l_total + ::thickness)/denom;

            // Apply the displacements.
            for (size_t c = 0; c < components.size(); c++) {
                const vector<Node*> &component = components[c];
                Vec3 displacement = -step_length*gc[c]/mc[c];
                for (size_t n = 0; n < component.size(); n++) {
                    Node *node = component[n];
                    if (!is_free(node) && !deform_obstacles)
                        continue;
                    node->x += displacement;
                }
            }

            // Update the world-space information.
            for (size_t m = 0; m < meshes.size(); m++) {
                compute_ws_data(*meshes[m]);
                update_accel_struct(*accs[m]);
            }
            for (size_t o = 0; o < obs_meshes.size(); o++)
                compute_ws_data(*obs_meshes[o]);
            for (size_t o = 0; o < obs_accs.size(); o++)
                update_accel_struct(*obs_accs[o]);            
        }

        if (deform_obstacles)
            ::obs_mass /= 2;
    }

    // Wrap-up.
    for (int m = 0; m < (int)meshes.size(); m++) {
        // Save collision-free configurations.
        update_x0(*meshes[m]);
    }
    for (int o = 0; o < (int)obs_meshes.size(); o++) {
        // Save collision-free configurations.
        update_x0(*obs_meshes[o]);
    }
    destroy_accel_structs(accs);
    destroy_accel_structs(obs_accs);
    // std::cout << "Separation done." << std::endl;
}

void add_gradient (const Face *face0, const Edge *edge0, const Face *face1,
                   Vec3 *g0, Vec3 *g1);

double compute_coplanar (const Face *face0, const Edge *edge0,
                         const Face *face1, Vec3 *g0, Vec3 *g1);

void compute_length_and_gradient (Ixn &ixn) {
    ixn.l = 0;
    for (int i = 0; i < 3; i++) {
        ixn.g0[i] = Vec3(0);
        ixn.g1[i] = Vec3(0);
    }
    if (norm2(cross(ixn.f0->n, ixn.f1->n)) < 1e-12) {
        // coplanar case
        Vec3 x0 = ixn.f0->v[0]->node->x,
             x1 = ixn.f1->v[0]->node->x,
             n = ixn.f0->n; // which is equal to ixn.f1->n
        if (abs(dot(x1 - x0, n)) > ::thickness)
            return;
        for (int f = 0; f < 2; f++) {
            const Face *face0 = (f==0) ? ixn.f0 : ixn.f1,
                       *face1 = (f==0) ? ixn.f1 : ixn.f0;
            Vec3 *g0 = (f==0) ? ixn.g0 : ixn.g1,
                 *g1 = (f==0) ? ixn.g1 : ixn.g0;
            for (int i = 0; i < 3; i++)
                ixn.l += compute_coplanar(face0, face0->adje[i], face1, g0, g1);
        }
    } else {
        // general case
        Line3 line;
        bool is_intersect = get_line_of_intersection(ixn.f0, ixn.f1, line);
        if (!is_intersect)
            return;
        ixn.l = norm(line.second - line.first);
        for (int f = 0; f < 2; f++) {
            const Face *face0 = (f==0) ? ixn.f0 : ixn.f1,
                       *face1 = (f==0) ? ixn.f1 : ixn.f0;
            Vec3 *g0 = (f==0) ? ixn.g0 : ixn.g1,
                 *g1 = (f==0) ? ixn.g1 : ixn.g0;
            for (int i = 0; i < 3; i++)
                add_gradient(face0, face0->adje[i], face1, g0, g1);
        }
    }
}

Bary barycentric_coords (const Vec3 &x, const Face *face);

void add_gradient (const Face *face0, const Edge *edge0, const Face *face1,
                   Vec3 *g0, Vec3 *g1) {
    Vec3 pt_ixn;
    bool is_ixn = edge_face_intersection(edge0, face1, pt_ixn);
    if (!is_ixn)
        return;
    Vec3 e = normalize(edge0->n[1]->x - edge0->n[0]->x);
    if (face0 != edge0->adjf[0])
        e = -e; // make sure edge vector is CCW w.r.t. face0
    // Get vector parallel to intersection line in correct direction.
    const Vec3 &n0 = face0->n, &n1 = face1->n;
    Vec3 r = normalize(cross(n0, n1));
    if (!right_handed(e, r, n0))
        r = -r;
    if (abs(dot(e, n1)) < 1e-12)
        return; // edge is parallel to face, ignore
    Vec3 face_grad = -dot(e,r)/dot(e,n1) * n1,
         edge_grad = -r - face_grad;
    Bary b0 = barycentric_coords(pt_ixn, face0),
         b1 = barycentric_coords(pt_ixn, face1);
    for (int v = 0; v < 3; v++) {
        g0[v] += b0[v]*edge_grad;
        g1[v] += b1[v]*face_grad;
    }
}

struct EdgeClipping {
    double t[2];
    Edge *edge[2]; // edges causing clip
    EdgeClipping () {t[0] = 0; t[1] = 1; edge[0] = edge[1] = NULL;}
    bool empty () {return t[0] > t[1];}
};

EdgeClipping clip_edge_to_face (const Edge *edge, const Face *face);

double compute_coplanar (const Face *face0, const Edge *edge0,
                         const Face *face1, Vec3 *g0, Vec3 *g1) {
    Vec3 x0 = edge0->n[0]->x, x1 = edge0->n[1]->x;
    Vec3 n = face1->n;
    EdgeClipping clip = clip_edge_to_face(edge0, face1);
    if (clip.empty()) // edge0 entirely outside face1
        return 0;
    for (int c = 0; c < 2; c++) {
        Vec3 r = normalize(edge0->n[1-c]->x - edge0->n[c]->x);
        const Edge *edge1 = clip.edge[c];
        if (!edge1) { // endpoint inside face1
            Bary b = barycentric_coords(edge0->n[c]->x, face0);
            for (int v = 0; v < 3; v++)
                g0[v] += b[v]*(-r);
            continue;
        } else {
            Vec3 e = normalize(edge1->n[1]->x - edge1->n[0]->x);
            Vec3 rp = cross(n,r), ep = cross(n,e);
            if (abs(dot(r,ep)) < 1e-12)
                continue; // edges are parallel, ignore
            Vec3 edge0_grad = dot(rp,ep)/dot(r,ep) * rp,
                 edge1_grad = -1/dot(r,ep) * ep;
            Vec3 x = x0 + clip.t[c]*(x1 - x0);
            Bary b0 = barycentric_coords(x, face0),
                 b1 = barycentric_coords(x, face1);
            for (int v = 0; v < 3; v++) {
                g0[v] += b0[v]*edge0_grad;
                g1[v] += b1[v]*edge1_grad;
            }
        }
    }
    return (clip.t[1] - clip.t[0])*norm(edge0->n[1]->x - edge0->n[0]->x);
}

EdgeClipping clip_edge_to_face (const Edge *edge, const Face *face) {
    EdgeClipping clip;
    Vec3 x0 = edge->n[0]->x, x1 = edge->n[1]->x;
    Vec3 n = face->n;
    for (int i = 0; i < 3; i++) {
        const Edge *edge1 = face->adje[i];
        Vec3 y0 = edge1->n[0]->x, y1 = edge1->n[1]->x;
        Vec3 e = normalize(y1 - y0);
        double a0 = stp(e, x0 - y0, n),
               a1 = stp(e, x1 - y0, n);
        if (stp(e, face->v[i]->node->x - y0, n) < 0) {
            a0 = -a0;
            a1 = -a1;
        } // a >= 0 if inside face
        if (a0 < 0 && a1 < 0) { // edge is completely outside
            clip.t[0] = 1;
            clip.t[1] = 0;
            return clip;
        } else if (a0 > 0 && a1 > 0) { // edge is completely inside
            continue;
        } else {
            for (int c = 0; c < 2; c++) {
                double ac = a0 + clip.t[c]*(a1 - a0);
                if (ac < 0) {
                    clip.t[c] = -a0/(a1 - a0);
                    clip.edge[c] = (Edge*)edge1;
                }
            }
        }
    }
    return clip;
}

Bary barycentric_coords (const Vec3 &x, const Face *face) {
    Bary b;
    double sum = 0;
    for (int i = 0; i < 3; i++) {
        b[i] = dot(face->n, cross(face->v[NEXT(i)]->node->x - x,
                                  face->v[PREV(i)]->node->x - x));
        b[i] = max(b[i], 0.);
        sum += b[i];
    }
    return b/sum;
}

Vec3 pos (const Face *face, const Bary &b) {
    return b[0]*face->v[0]->node->x
         + b[1]*face->v[1]->node->x
         + b[2]*face->v[2]->node->x;
}

void update_active (const vector<AccelStruct*> &accs,
                    const vector<AccelStruct*> &obs_accs,
                    const vector<Ixn> &ixns) {
    assert(false);
    // for (int a = 0; a < accs.size(); a++)
    //     mark_all_inactive(*accs[a]);
    // for (int i = 0; i < ixns.size(); i++)
    //     for (int f = 0; f < 2; f++) {
    //         const Face *face = f==0 ? ixns[i].f0 : ixns[i].f1;
    //         int m = find(face, *::meshes);
    //         if (m == -1)
    //             continue;
    //         mark_active(*accs[m], face);
    //     }
}

static int nthreads = 0;
static vector<Ixn> *ixns = NULL;

void find_face_intersection (const Face *face0, const Face *face1);

vector<Ixn> find_intersections (const vector<AccelStruct*> &accs,
                                const vector<AccelStruct*> &obs_accs) {
    if (!::ixns) {
        ::nthreads = omp_get_max_threads();
        ::ixns = new vector<Ixn>[::nthreads];
    }
    for (int t = 0; t < ::nthreads; t++)
        ::ixns[t].clear();
    for_overlapping_faces(accs, obs_accs, ::thickness, find_face_intersection);
    vector<Ixn> ixns;
    for (int t = 0; t < ::nthreads; t++)
        append(ixns, ::ixns[t]);
    return ixns;
}

void find_face_overlappings (const Face *face0, const Face *face1);

// Find all the pairs of faces that might be intersections between them.
vector<Ixn> find_overlappings (const vector<AccelStruct*> &accs,
        const vector<AccelStruct*> &obs_accs) {
    if (!::ixns) {
        ::nthreads = omp_get_max_threads();
        ::ixns = new vector<Ixn>[::nthreads];
    }
    for (int t = 0; t < ::nthreads; t++)
        ::ixns[t].clear();
    for_overlapping_faces(accs, obs_accs, ::thickness, find_face_overlappings);
    vector<Ixn> ixns;
    for (int t = 0; t < ::nthreads; t++)
        append(ixns, ::ixns[t]);
    return ixns;
}

bool adjacent (const Face *face0, const Face *face1);

bool intersection_midpoint (const Face *face0, const Face *face1,
                            Bary &b0, Bary &b1);

void find_face_intersection (const Face *face0, const Face *face1) {
    if (adjacent(face0, face1))
        return;
    int t = omp_get_thread_num();
    Ixn ixn(face0, face1);
    compute_length_and_gradient(ixn);
    if (ixn.l == 0)
        return;
    ::ixns[t].push_back(ixn);
}

void find_face_overlappings (const Face *face0, const Face *face1) {
    if (adjacent(face0, face1))
        return;
    int t = omp_get_thread_num();
    ::ixns[t].push_back(Ixn(face0, face1));
}

bool adjacent (const Face *face0, const Face *face1) {
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            if (face0->v[i]->node == face1->v[j]->node)
                return true;
    return false;
}

bool face_plane_intersection (const Face *face, const Face *plane,
                              Bary &b0, Bary &b1);
int major_axis (const Vec3 &v);

bool intersection_midpoint (const Face *face0, const Face *face1,
                            Bary &b0, Bary &b1) {
    if (norm2(cross(face0->n, face1->n)) < 1e-12)
        return false;
    Bary b00, b01, b10, b11;
    bool ix0 = face_plane_intersection(face0, face1, b00, b01),
         ix1 = face_plane_intersection(face1, face0, b10, b11);
    if (!ix0 || !ix1)
        return false;
    int axis = major_axis(cross(face0->n, face1->n));
    double a00 = pos(face0, b00)[axis], a01 = pos(face0, b01)[axis],
           a10 = pos(face1, b10)[axis], a11 = pos(face1, b11)[axis];
    double amin = max(min(a00, a01), min(a10, a11)),
           amax = min(max(a00, a01), max(a10, a11)),
           amid = (amin + amax)/2;
    if (amin > amax)
        return false;
    b0 = (a01==a00) ? b00 : b00 + (amid-a00)/(a01-a00)*(b01-b00);
    b1 = (a11==a10) ? b10 : b10 + (amid-a10)/(a11-a10)*(b11-b10);
    return true;
}

struct cmpOneAxis {
    size_t axis;

    cmpOneAxis(size_t axis)
        : axis(axis) {}
    bool operator() (const Vec3& lhs, const Vec3& rhs) const {
        return lhs[axis] < rhs[axis];
    }
};

bool get_line_of_intersection (const Face *face0, const Face *face1, Line3& line) {
    if (norm2(cross(face0->n, face1->n)) < 1e-12)
        return false;

    Bary b00, b01, b10, b11;
    bool ix0 = face_plane_intersection(face0, face1, b00, b01),
         ix1 = face_plane_intersection(face1, face0, b10, b11);
    if (!ix0 || !ix1)
        return false;

    const Vec3& p00 = pos(face0, b00);
    const Vec3& p01 = pos(face0, b01);
    const Vec3& p10 = pos(face1, b10);
    const Vec3& p11 = pos(face1, b11);

    int axis = major_axis(cross(face0->n, face1->n));
    double a00 = p00[axis], a01 = p01[axis],
           a10 = p10[axis], a11 = p11[axis];

    // If there is no overlapping between the two lines of intersection,
    if (max(a00, a01) < min(a10, a11) || max(a10, a11) < min(a00, a01))
        return false;

    // Sort to pick the overlapping region.
    std::vector<Vec3> pts;
    pts.push_back(p00);
    pts.push_back(p01);
    pts.push_back(p10);
    pts.push_back(p11);
    std::sort(pts.begin(), pts.end(), cmpOneAxis(axis));

    line = std::make_pair(pts[2], pts[1]);
    return true;
}

bool face_plane_intersection (const Face *face, const Face *plane,
                              Bary &b0, Bary &b1) {
    const Vec3 &x0 = plane->v[0]->node->x, &n = plane->n;
    double h[3];
    int sign_sum = 0;
    for (int v = 0; v < 3; v++) {
        h[v] = dot(face->v[v]->node->x - x0, n);
        sign_sum += sgn(h[v]);
    }
    if (sign_sum == -3 || sign_sum == 3)
        return false;
    int v0 = -1;
    for (int v = 0; v < 3; v++)
        if (sgn(h[v]) == -sign_sum)
            v0 = v;
    double t0 = h[v0]/(h[v0] - h[NEXT(v0)]), t1 = h[v0]/(h[v0] - h[PREV(v0)]);
    b0[v0] = 1 - t0;
    b0[NEXT(v0)] = t0;
    b0[PREV(v0)] = 0;
    b1[v0] = 1 - t1;
    b1[PREV(v0)] = t1;
    b1[NEXT(v0)] = 0;
    return true;
}

bool edge_face_intersection (const Vec3& e0, const Vec3& e1,
        const Vec3& x0, const Vec3& x1, const Vec3& x2, const Vec3& n, Vec3& pt) {
    const double numer = dot(x0 - e0, n);
    const double denom = dot(e1 - e0, n);
    if (std::abs(numer) < 1e-12 || std::abs(denom) < 1e-12)
        return false;

    // Find a parameter t in [0, 1].
    const double t = numer / denom;
    if (t < 0.0 || t > 1.0)
        return false;

    // Find a point of intersection.
    pt = e0 + t * (e1 - e0);

    // Check CCW (insdie) or CW (outside).
    const Vec3& c0 = cross(x1 - x0, pt - x0);
    const Vec3& c1 = cross(x2 - x1, pt - x1);
    const Vec3& c2 = cross(x0 - x2, pt - x2);

    bool is_inside = (dot(n, c0) >= -1e-12
            && dot(n, c1) >= -1e-12
            && dot(n, c2) >= -1e-12);
    return is_inside;
}

bool edge_face_intersection (const Edge* edge, const Face *face, Vec3& pt) {
    // edge
    const Vec3& e0 = edge->n[0]->x;
    const Vec3& e1 = edge->n[1]->x;

    // face
    const Vec3& f0 = face->v[0]->node->x;
    const Vec3& f1 = face->v[1]->node->x;
    const Vec3& f2 = face->v[2]->node->x;
    const Vec3& fn = normal<WS>(face);

    return edge_face_intersection(e0, e1, f0, f1, f2, fn, pt);
}

int major_axis (const Vec3 &v) {
    return (abs(v[0]) > abs(v[1]) && abs(v[0]) > abs(v[2])) ? 0
         : (abs(v[1]) > abs(v[2])) ? 1 : 2;
}

struct UnionFind {
    vector<size_t> parent, rank;
    UnionFind (size_t n): parent(n), rank(n, 0) {
        for (size_t i = 0; i < n; i++) parent[i] = i;
    }
    size_t find (size_t i) {
        if (parent[i] != i)
            parent[i] = find(parent[i]);
        return parent[i];
    }
    void unify (size_t x, size_t y) {
        size_t x0 = find(x), y0 = find(y);
        if (x0 == y0)
            return;
        if (rank[x0] < rank[y0])
            parent[x0] = y0;
        else if (rank[x0] > rank[y0])
            parent[y0] = x0;
        else {
            parent[y0] = x0;
            rank[x0]++;
        }
    }
};

vector< vector<Node*> > connected_components (const vector<Ixn> &ixns) {
    vector<Node*> nodes;
    for (size_t i = 0; i < ixns.size(); i++) {
        for (int v = 0; v < 3; v++) {
            include(ixns[i].f0->v[v]->node, nodes);
            include(ixns[i].f1->v[v]->node, nodes);
        }
    }
    for (size_t n = 0; n < nodes.size(); n++)
        nodes[n]->index = n;
    UnionFind uf(nodes.size());
    for (size_t i = 0; i < ixns.size(); i++) {
        const Face *face0 = ixns[i].f0, *face1 = ixns[i].f1;
        for (int v = 0; v < 3; v++) {
            uf.unify(face0->v[v]->node->index, face0->v[NEXT(v)]->node->index);
            uf.unify(face1->v[v]->node->index, face1->v[NEXT(v)]->node->index);
        }
    }
    int c = 0;
    for (size_t i = 0; i < nodes.size(); i++)
        if (uf.find(i) == i)
            nodes[i]->index = c++;
    vector< vector<Node*> > components(c);
    for (size_t i = 0; i < nodes.size(); i++)
        components[nodes[uf.find(i)]->index].push_back(nodes[i]);
    return components;
}

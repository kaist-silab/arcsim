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

#include "proximity.hpp"

#include "collisionutil.hpp"
#include "geometry.hpp"
#include "magic.hpp"
#include "io.hpp"
#include "proxy.hpp"
#include "simulation.hpp"
#include <vector>
using namespace std;

template <typename T> struct Min {
    double key;
    T val;
    Min (): key(infinity), val() {}
    void add (double key, T val) {
        if (key < this->key) {
            this->key = key;
            this->val = val;
        }
    }
};

// useful for debugging
extern bool consistency_check;

template<class T> void serialize_minvec(vector<Min<T*> >& v, Serialize& s, const string& name) {
	serializer_array(v, s, name);
	for(int i=0; i<v.size(); i++) {
		serializer(v[i].key, s, name);
		int idx = v[i].val ? v[i].val->index : -1;
		serializer(idx, s, name);		
	}
}
template<> void serializer<vector<Min<Face*> > >(vector<Min<Face*> >& x, Serialize& s, const string& n) { serialize_minvec(x,s,n); }
template<> void serializer<vector<Min<Edge*> > >(vector<Min<Edge*> >& x, Serialize& s, const string& n) { serialize_minvec(x,s,n); }
template<> void serializer<vector<Min<Node*> > >(vector<Min<Node*> >& x, Serialize& s, const string& n) { serialize_minvec(x,s,n); }

static vector< Min<Face*> > node_prox[2];
static vector< Min<Edge*> > edge_prox[2];
static vector< Min<Node*> > face_prox[2];
static vector< Min<Node*> > edge_node_prox;
static vector< Min<Edge*> > node_edge_prox;

void find_proximities (const Face *face0, const Face *face1);
Constraint *make_constraint (const Node *node, const Face *face,
                             double mu, double mu_obs);
Constraint *make_constraint (const Edge *edge0, const Edge *edge1,
                             double mu, double mu_obs);
Constraint *make_constraint (const Edge *edge, const Node *node,
                             double mu, double mu_obs);
void make_proxy_constraints (Mesh& mesh, CollisionProxy& proxy, vector<Constraint*>& cons);

vector<Constraint*> proximity_constraints (vector<Mesh*> &meshes,
                                           const vector<Mesh*> &obs_meshes,
                                           double mu, double mu_obs, bool proxy_only, bool obs_only) {
    ::meshes = &meshes;
    const double dmin = 2*::magic.repulsion_thickness;
    vector<Constraint*> cons;
    
    if (proxy_only || obs_only) {
        for (size_t m = 0; m<meshes.size(); m++) {
            for (size_t i=0; i<obs_meshes.size(); i++) {
                if (obs_meshes[i]->proxy)
                    make_proxy_constraints(*meshes[m], * (CollisionProxy*)(obs_meshes[i]->proxy), cons);
            }
        }
        if (proxy_only)
            return cons;
    }

    vector<AccelStruct*> accs = create_accel_structs(meshes, false),
                         obs_accs = create_accel_structs(obs_meshes, false);

    set_indices(meshes);
    int nn=0, ne=0, nf=0;
    for (size_t m=0; m<meshes.size(); m++) {
    	nn += meshes[m]->nodes.size();
    	ne += meshes[m]->edges.size();
    	nf += meshes[m]->faces.size();
    }
    for (int i = 0; i < 2; i++) {
        ::node_prox[i].assign(nn, Min<Face*>());
        ::edge_prox[i].assign(ne, Min<Edge*>());
        ::face_prox[i].assign(nf, Min<Node*>());
    }
    ::edge_node_prox.assign(ne, Min<Node*>());
    ::node_edge_prox.assign(nn, Min<Edge*>());

    for_overlapping_faces(accs, obs_accs, dmin, find_proximities, true, obs_only);

    for (size_t m = 0; m<meshes.size(); m++) {
    	Mesh& mesh = *meshes[m];

	    for (size_t n = 0; n < mesh.nodes.size(); n++) {
	    	int idx = mesh.nodes[n]->index;
	        for (int i = 0; i < 2; i++) {
	            Min<Face*> &m = ::node_prox[i][idx];
	            if (m.key < dmin)
	                cons.push_back(make_constraint(mesh.nodes[n], m.val, mu, mu_obs));
            }
            Min<Edge*> &me = ::node_edge_prox[idx];
            if (me.key < dmin)
                cons.push_back(make_constraint(me.val, mesh.nodes[n], mu, mu_obs));	        
	    }
	    for (size_t e = 0; e < mesh.edges.size(); e++) {
	        int idx = mesh.edges[e]->index;
	        for (int i = 0; i < 2; i++) {
	            Min<Edge*> &m = ::edge_prox[i][idx];
	            if (m.key < dmin)
	                cons.push_back(make_constraint(mesh.edges[e], m.val, mu, mu_obs));
            }
            Min<Node*> &me = ::edge_node_prox[idx];
            if (me.key < dmin)
                cons.push_back(make_constraint(mesh.edges[e], me.val, mu, mu_obs));            
	    }
	    for (size_t f = 0; f < mesh.faces.size(); f++) {
	        int idx = mesh.faces[f]->index;
	        for (int i = 0; i < 2; i++) {
	            Min<Node*> &m = ::face_prox[i][idx];
	            if (m.key < dmin)
	                cons.push_back(make_constraint(m.val, mesh.faces[f], mu, mu_obs));
	        }
	    }

        for (size_t i=0; i<obs_meshes.size(); i++) {
            if (obs_meshes[i]->proxy)
                make_proxy_constraints(mesh, * (CollisionProxy*)(obs_meshes[i]->proxy), cons);
        }
	}


    if (consistency_check) {
		cout << "> proximity2: " << endl;
		test_state(cons, "/tmp/cs");
	}

    destroy_accel_structs(accs);
    destroy_accel_structs(obs_accs);
    return cons;
}

void add_proximity (const Node *node, const Face *face);
void add_proximity (const Edge *edge0, const Edge *edge1);
void add_proximity (Node *node, Edge* edge);

void find_proximities (const Face *face0, const Face *face1) {
    for (int v = 0; v < 3; v++)
        add_proximity(face0->v[v]->node, face1);
    for (int v = 0; v < 3; v++)
        add_proximity(face1->v[v]->node, face0);
    for (int e0 = 0; e0 < 3; e0++)
        for (int e1 = 0; e1 < 3; e1++) {
            add_proximity(face0->adje[e0], face1->adje[e1]);
            add_proximity(face0->v[e0]->node, face1->adje[e1]);
            add_proximity(face1->v[e0]->node, face0->adje[e1]);
        }
}

static inline bool has_node(const Face* f, const Node* n) {
    return f && (f->v[0]->node == n || f->v[1]->node == n || f->v[2]->node == n);
}

static inline Vec3 get_outwards_normal(const Edge *e) {
    Vec3 p0 = e->n[0]->x, p1 = e->n[1]->x, ne = normalize(p1-p0);
    Vec3 opp = p0 - edge_opp_vert(e, e->adjf[0] ? 0 : 1)->node->x;
    Vec3 n = cross(e->adjf[0] ? e->adjf[0]->n : e->adjf[1]->n, ne);
    if (dot(n,opp) < 0) n= -n;
    return n;
}

// half-cylinder test
void add_proximity(Node *node, Edge* edge) {
    if (node == edge->n[0] || node == edge->n[1] ||
        !is_seam_or_boundary(node) || !is_seam_or_boundary(edge) ||
        (edge->adjf[0] && has_node(edge->adjf[0], node)) ||
        (edge->adjf[1] && has_node(edge->adjf[1], node))) return;
    
    Vec3 normal = get_outwards_normal(edge);
    Vec3 p0 = edge->n[0]->x, p1 = edge->n[1]->x, x = node->x;
    if (dot(x-p0,normal) < 0) 
        return;
    
    double d = dot(x-p0,p1-p0)/norm2(p1-p0);
    if (d < 0 || d > 1)
        return;
    
    double w0 = 1.0-d, w1 = d;
    double dist = norm(w0*p0 + w1*p1 - x);
    if (is_free(node))
        ::node_edge_prox[node->index].add(dist, edge);
    if (is_free(edge))
        ::edge_node_prox[edge->index].add(dist, node);
}

void add_proximity (const Node *node, const Face *face) {
    if (has_node(face,node)) return;
    Vec3 n;
    double w[4];
    double d = signed_vf_distance(node->x, face->v[0]->node->x,
                                  face->v[1]->node->x, face->v[2]->node->x,
                                  &n, w);
    d = abs(d);
    bool inside = (min(-w[1], -w[2], -w[3]) >= -1e-6);
    if (!inside)
        return;
    if (is_free(node)) {
        int side = dot(n, node->n)>=0 ? 0 : 1;
        ::node_prox[side][node->index].add(d, (Face*)face);
    }
    if (is_free(face)) {
        int side = dot(-n, face->n)>=0 ? 0 : 1;
        ::face_prox[side][face->index].add(d, (Node*)node);
    }
}

bool in_wedge (double w, const Edge *edge0, const Edge *edge1) {
    Vec3 x = (1-w)*edge0->n[0]->x + w*edge0->n[1]->x;
    bool in = true;
    for (int s = 0; s < 2; s++) {
        const Face *face = edge1->adjf[s];
        if (!face)
            continue;
        const Node *node0 = edge1->n[s], *node1 = edge1->n[1-s];
        Vec3 e = node1->x - node0->x, n = face->n, r = x - node0->x;
        in &= right_handed(e, n, r);
    }
    return in;
}

void add_proximity (const Edge *edge0, const Edge *edge1) {
    if (edge0->n[0] == edge1->n[0] || edge0->n[0] == edge1->n[1]
     || edge0->n[1] == edge1->n[0] || edge0->n[1] == edge1->n[1])
        return;
    Vec3 n;
    double w[4];
    double d = signed_ee_distance(edge0->n[0]->x, edge0->n[1]->x,
                                  edge1->n[0]->x, edge1->n[1]->x,
                                  &n, w);
    d = abs(d);
    bool inside = (min(w[0], w[1], -w[2], -w[3]) >= -1e-6
                   && in_wedge(w[1], edge0, edge1)
                   && in_wedge(-w[3], edge1, edge0));
    if (!inside)
        return;
    if (is_free(edge0)) {
        Vec3 edge0n = edge0->n[0]->n + edge0->n[1]->n;
        int side = dot(n, edge0n)>=0 ? 0 : 1;
        ::edge_prox[side][edge0->index].add(d, (Edge*)edge1);
    }
    if (is_free(edge1)) {
        Vec3 edge1n = edge1->n[0]->n + edge1->n[1]->n;
        int side = dot(-n, edge1n)>=0 ? 0 : 1;
        ::edge_prox[side][edge1->index].add(d, (Edge*)edge0);
    }
}

double area_cached (const Node *node);
double area_cached (const Edge *edge);
double area_cached (const Face *face);

Constraint *make_constraint (const Node *node, const Face *face,
                             double mu, double mu_obs) {
    IneqCon *con = new IneqCon;
    con->nodes[0] = (Node*)node;
    con->nodes[1] = (Node*)face->v[0]->node;
    con->nodes[2] = (Node*)face->v[1]->node;
    con->nodes[3] = (Node*)face->v[2]->node;
    for (int n = 0; n < 4; n++)
        con->free[n] = is_free(con->nodes[n]);
    double a = min(area_cached(node), area_cached(face));
    con->stiff = ::magic.collision_stiffness*a;
    double d = signed_vf_distance(con->nodes[0]->x, con->nodes[1]->x,
                                  con->nodes[2]->x, con->nodes[3]->x,
                                  &con->n, con->w);
    if (d < 0)
        con->n = -con->n;
    con->mu = (!is_free(node) || !is_free(face)) ? mu_obs : mu;
    return con;
}

Constraint *make_constraint (const Edge *edge0, const Edge *edge1,
                             double mu, double mu_obs) {
    IneqCon *con = new IneqCon;
    con->nodes[0] = (Node*)edge0->n[0];
    con->nodes[1] = (Node*)edge0->n[1];
    con->nodes[2] = (Node*)edge1->n[0];
    con->nodes[3] = (Node*)edge1->n[1];
    for (int n = 0; n < 4; n++)
        con->free[n] = is_free(con->nodes[n]);
    double a = min(area_cached(edge0), area_cached(edge1));
    con->stiff = ::magic.collision_stiffness*a;
    double d = signed_ee_distance(con->nodes[0]->x, con->nodes[1]->x,
                                  con->nodes[2]->x, con->nodes[3]->x,
                                  &con->n, con->w);
    if (d < 0)
        con->n = -con->n;
    con->mu = (!is_free(edge0) || !is_free(edge1)) ? mu_obs : mu;
    return con;
}

Constraint *make_constraint (const Edge* edge, const Node* node, double mu, double mu_obs) {
    IneqCon *con = new IneqCon;
    con->nodes[0] = (Node*)node;
    con->nodes[1] = (Node*)edge->n[0];
    con->nodes[2] = (Node*)edge->n[1];
    con->nodes[3] = 0;
    for (int n = 0; n < 4; n++)
        con->free[n] = con->nodes[n] ? is_free(con->nodes[n]) : false;
    
    double a = min(area_cached(edge), area_cached(node));
    con->stiff = ::magic.collision_stiffness*a;
    con->n = get_outwards_normal(edge);
    double d = signed_ve_distance(con->nodes[0]->x,con->nodes[1]->x,con->nodes[2]->x,
                                  &con->n, con->w);
    
    if (fabs(d) > 2.0* ::magic.repulsion_thickness) return 0;
    
    con->mu = (!is_free(node) || !is_free(edge)) ? mu_obs : mu;
    return con;
}

void make_proxy_constraints (Mesh& mesh, CollisionProxy& proxy, vector<Constraint*>& cons) {
    for (size_t i=0; i<mesh.nodes.size(); i++) {
        Constraint* c = proxy.constraint(mesh.nodes[i]);
        if (c)
            cons.push_back(c);
    }
}


double area_cached (const Node *node) {
    if (is_free(node))
    	return node->a;
    double a = 0;
    for (int v = 0; v < (int)node->verts.size(); v++)
        for (int f = 0; f < (int)node->verts[v]->adjf.size(); f++)
            a += area(node->verts[v]->adjf[f])/3;
    return a;
}

double area_cached (const Edge *edge) {
    double a = 0;
    if (edge->adjf[0])
        a += area(edge->adjf[0])/3;
    if (edge->adjf[1])
        a += area(edge->adjf[1])/3;
    return a;
}

double area_cached (const Face *face) {
    if (is_free(face))
        return face->a;
    const Vec3 &x0 = face->v[0]->node->x, &x1 = face->v[1]->node->x,
               &x2 = face->v[2]->node->x;
    return norm(cross(x1-x0, x2-x0))/2;
}

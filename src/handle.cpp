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

#include "handle.hpp"
#include "magic.hpp"
#include "display.hpp"
using namespace std;

static Vec3 directions[3] = {Vec3(1,0,0), Vec3(0,1,0), Vec3(0,0,1)};

void add_position_constraints (const Node *node, const Vec3 &x, double stiff,
                               vector<Constraint*> &cons);

Transformation normalize (const Transformation &T) {
    Transformation T1 = T;
    T1.rotation = normalize(T1.rotation);
    return T1;
}

vector<Constraint*> NodeHandle::get_constraints (double t) {
    double s = strength(t);
    if (!s)
        return vector<Constraint*>();
    if (!activated) {
        // handle just got started, fill in its original position
        x0 = motion ? inverse(normalize(motion->pos(t))).apply(node->x) : node->x;
        activated = true;
    }
    Vec3 x = motion ? normalize(motion->pos(t)).apply(x0) : x0;
    vector<Constraint*> cons;
    add_position_constraints(node, x, s*::magic.handle_stiffness, cons);
    return cons;
}

vector<Constraint*> CircleHandle::get_constraints (double t) {
    double s = strength(t);
    if (!s)
        return vector<Constraint*>();
    vector<Constraint*> cons;
    for (int n = 0; n < (int)mesh->nodes.size(); n++) {
        Node *node = mesh->nodes[n];
        if (node->label != label)
            continue;
        double theta = 2*M_PI*dot(reduce_xy(node->verts[0]->u), u)/c;
        Vec3 x = xc + (dx0*cos(theta) + dx1*sin(theta))*c/(2*M_PI);
        if (motion)
            x = motion->pos(t).apply(x);
        double l = 0;
        for (int e = 0; e < (int)node->adje.size(); e++) {
            const Edge *edge = node->adje[e];
            if (edge->n[0]->label != label || edge->n[1]->label != label)
                continue;
            l += norm(edge->n[1]->x - edge->n[0]->x);
        }
        add_position_constraints(node, x, s*::magic.handle_stiffness*l, cons);
    }
    return cons;
}

vector<Constraint*> GlueHandle::get_constraints (double t) {
    double s = strength(t);
    if (!s)
        return vector<Constraint*>();
    vector<Constraint*> cons;
    for (int i = 0; i < 3; i++) {
        GlueCon *con = new GlueCon;
        con->nodes[0] = nodes[0];
        con->nodes[1] = nodes[1];
        con->n = directions[i];
        con->stiff = s*::magic.handle_stiffness;
        cons.push_back(con);
    }
    return cons;
}

vector<Constraint*> SoftHandle::get_constraints (double t) {
    return vector<Constraint*>();
}
    
void SoftHandle::add_forces(double t, vector<Vec3> &fext, vector<Mat3x3>& Jext) {
	double s = strength(t);
	if (!s) return;
	Annotation::add((motion) ? motion->pos(t).apply(center) : center); 

    for (size_t v = 0; v < mesh->verts.size(); v++) {
        Vert *vert = mesh->verts[v];

        double d = norm(vert->u - center)/radius;
        if (d > 1.0 || !vert->node->active()) continue;
        int ix = v;
        
        Vec3 x0 = (motion) ? motion->pos(t).apply(vert->u) : vert->u; 
        Vec3 dx = vert->node->x-x0;

        double w = (1.0 + sq(d) * (-3.0 + 2.0*d)) * ::magic.handle_stiffness * s *vert->node->a;
        //w=0;
        fext[ix] += -w * dx;
        Jext[ix] += Mat3x3(-w);

        vert->node->flag |= Node::FlagResolveUni;
    }    
}

vector<Node*> SoftHandle::get_nodes() {
    std::vector<Node*> nodes;
    for (size_t v = 0; v < mesh->verts.size(); v++) {
        Vert *vert = mesh->verts[v];
        
        double d = norm2(vert->u - center)/sq(radius);
        if (d < 1.0)
            nodes.push_back(vert->node);
    }
    return nodes;
}

void add_position_constraints (const Node *node, const Vec3 &x, double stiff,
                               vector<Constraint*> &cons) {
    for (int i = 0; i < 3; i++) {
        EqCon *con = new EqCon;
        con->node = (Node*)node;
        con->x = x;
        con->n = directions[i];
        con->stiff = stiff;
        cons.push_back(con);
    }
}

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

 #include "breaking.hpp"
#include "sepstrength.hpp"
#include "util.hpp"
#include "cloth.hpp"
#include "magic.hpp"
#include "geometry.hpp"
#include "remesh.hpp"
#include "nearobs.hpp"
#include "physics.hpp"
#include "subset.hpp"
#include "proximity.hpp"
#include "dynamicremesh.hpp"
#include "collisionutil.hpp"
#include "display.hpp"
#include <map>
#include <queue>

using namespace std;

const int SUPPORT_RINGS = 2;

double edge_plane_intersect(const Vec3& e0, const Vec3& e1, const Vec3& u0, const Vec3& n) {
	double div = dot(e1-e0, n);
	if (div == 0) { cout << "parallel lines" << endl; exit(1);}
	return -dot(e0 - u0, n) / div;
}

Edge* split_edge_with_plane(Face* face, Node* node, Vec3& n) {
	int idx = (face->v[2]->node == node) ? 2 : (face->v[1]->node == node ? 1 : 0);
	Vert* e0 = face->v[PREV(idx)], *e1 = face->v[NEXT(idx)], *vert = face->v[idx];
	Edge* edge = get_edge(e0->node, e1->node);
	if (edge->n[0] != e0->node)
		swap(e0,e1);

	// get intersection point and angles
	double d = edge_plane_intersect(e0->u, e1->u, vert->u, cross(n,normal<MS>(face)));
    double l = norm(e1->u - e0->u);
    Vec3 split = (1.0-d)*e0->u + d*e1->u;
	double angle0 = get_angle(e0->u - vert->u, split - vert->u);
	double angle1 = get_angle(e1->u - vert->u, split - vert->u);

	const double min_length = node->mesh->parent->remeshing.size_min * 0.1;
    const double cut_min_length = node->mesh->parent->remeshing.size_min * 0.05;
    const double min_angle = 5 * M_PI/180;

    if (l < cut_min_length || area(face) < sqr(cut_min_length))
        return 0;
	
	// if close to existing edge, try to nudge vertex instead
	if (d < 0.2 && (d*l < min_length || angle0 < min_angle)) {
		try_move_node(e0->node, edge, d);			
	    return get_edge(node,e0->node);
	}
	else if (d > 0.8 && ((1-d)*l < min_length || angle1 < min_angle)) {
		try_move_node(e1->node, edge, 1-d);			
        return get_edge(node,e1->node);
	}
	else {
        RemeshOp op = split_edge(edge, d);
		Edge* new_edge = get_edge(node,op.added_nodes[0]);
		op.done();
        
		return new_edge;
	}
}

// TODO: displace MS to compensate collision size
void displace_sector(Node* node, double disp) {
	Edge* start = 0, *end = 0;
	for (size_t e=0; e<node->adje.size(); e++) {
		Edge* edge = node->adje[e];
		if (!next_face_cw(edge, node))
			start = edge;
		if (!next_face_ccw(edge, node))
			end = edge;
	}
	if (start || end) {
		Vec3 s0 = normalize(other_node(start, node)->x - node->x);
		Vec3 s1 = normalize(other_node(end, node)->x - node->x);
		node->x += normalize(cross(s1-s0,normal<WS>(node))) * disp;
	}
}

void split_sector(Node* node, Edge* start, Edge* end, MeshSubset& subset) {
	RemeshOp op;
	map<Vert*,Vert*> new_verts;

	// duplicate node
	Node* new_node = new Node(node->y,node->x,node->v,node->label,node->flag, true);
	node->preserve = true;
	op.added_nodes.push_back(new_node);
    //include(node, subset.active_nodes);
    //subset.active_nodes.push_back(new_node);

    // circle around split node, replacing edges and faces
	Vert* vert = 0;
	for(Edge* cur = start;;cur=next_edge_ccw(cur,node)) {
		Face* face = next_face_ccw(cur, node);
	
        // duplicate vertex if necessary
		if (face) {
			vert = get_vert(face, node);
			if (!new_verts[vert])
				new_verts[vert] = new Vert(vert->u);
		}

		// replace edge
		Edge* repl = new Edge(cur->n[0], cur->n[1], cur->theta_ideal, cur->preserve);
		for (int i=0; i<2; i++) if(repl->n[i] == node) repl->n[i] = new_node;
		op.added_edges.push_back(repl);

		if (cur == end)
			break;

		// don't remove first edge -> duplicate
		if (cur != start || !start->adjf[0] || !start->adjf[1])
			op.removed_edges.push_back(cur);
		
		if (!face)
			break;

		// replace face
		Face* new_face = new Face(face->v[0],face->v[1],face->v[2],face->Sp_str,face->Sp_bend,
								  face->material, face->damage);
		for (int i=0; i<3; i++) if(new_face->v[i]->node == node) new_face->v[i] = new_verts[vert];
		
		op.added_faces.push_back(new_face);
		op.removed_faces.push_back(face);
	}
	// add new vertices
	for (map<Vert*,Vert*>::iterator it = new_verts.begin(); it != new_verts.end(); ++it) {
		op.added_verts.push_back(it->second);
		connect(it->second, new_node);
	}
	op.apply(*node->mesh);
	op.done();
	
	// remove obsolete vertices
	for (map<Vert*,Vert*>::iterator it = new_verts.begin(); it != new_verts.end(); ++it) {
		if (it->first->adjf.empty())
			node->mesh->remove(it->first);
	}
	displace_sector(node, 1e-6);
	displace_sector(new_node, 1e-6);
}

void fix_nonmanifold(Node* node, MeshSubset& subset) {
	node->preserve = true;
	
	// search for the two crack start edges
	Edge *start = NULL;
	int num_start = 0;
	for (size_t e=0; e<node->adje.size(); e++) {
		Edge* edge = node->adje[e];
		if (!next_face_cw(edge, node)) {
			start = edge;
			num_start++;
		}
	}
	// check if really non-manifold
	if (num_start != 2) 
		return;

	split_sector(node, start, 0, subset);
}

// open a crack with the given splitplane
bool break_node(SplitNode& split, MeshSubset& subset) {
	Node* node = split.node;

	// generate edges to split along
	Edge* edges[2] = { 0, 0 };
	for (int i=0; i<2; i++) {
		if (split.faces[i]) {
			Edge* e = split_edge_with_plane(split.faces[i], node, split.normal[i]);
            if (!e) return false;
            
            edges[edges[0] ? 1:0] = e;
		}
	}
	if (!edges[0])
		return false;

	// first start edge: existing crack, or one of the new edges
	Edge *start = NULL;
	for (size_t e=0; e<node->adje.size(); e++) {
		Edge* edge = node->adje[e];
		if (!next_face_cw(edge, node) || edge == edges[1])
			start = edge;
	}
    if (is_seam_or_boundary(edges[0]) && is_seam_or_boundary(start))
        return false;

	Node *end0 = other_node(edges[0], node);
	Node *end1 = edges[1] ? other_node(edges[1], node) : 0;
           
    split_sector(node, start, edges[0], subset);

	if (end0) fix_nonmanifold(end0, subset);
	if (end1) fix_nonmanifold(end1, subset);
	if (end0) subset.active_nodes.push_back(end0);
    if (end1) subset.active_nodes.push_back(end1);
    
	return true;
}

void local_physics_step(MeshSubset& subset) {
    const double dt = min(1e-4, sim.step_time);

    // TODO: local constraints
    vector<Constraint*> cons;
    vector<Node*>& nodes = subset.active_nodes;
    vector<Edge*> edges = subset.get_edges();
    vector<Face*> faces = subset.get_faces();

    // could be optimized, so it only applies to nodes subset
    if (::magic.relax_method == 1)
        cons = proximity_constraints(sim.cloth_meshes, sim.obstacle_meshes,
                                    sim.friction, sim.obs_friction);

    // save old position, and advance boundaries
    for (size_t i=0; i<nodes.size(); i++)
        nodes[i]->x0 = nodes[i]->x;
    for (size_t i=0; i<subset.support_nodes.size(); i++) {
        subset.support_nodes[i]->x0 = subset.support_nodes[i]->x;
    }

    vector<Vec3> fext(nodes.size(), Vec3(0));
    vector<Mat3x3> Jext(nodes.size(), Mat3x3(0));

    activate_nodes(nodes);
    //add_external_forces(nodes, faces, sim.gravity, sim.wind, fext, Jext);
    //for (size_t h = 0; h < sim.handles.size(); h++)
    //    sim.handles[h]->add_forces(sim.time,fext,Jext);

    vector<Vec3> dv = implicit_update(nodes, edges, faces, fext, Jext, cons, dt);
    for (size_t n = 0; n < nodes.size(); n++) {
        nodes[n]->x += (dv[n]) * dt;
    }
    
    deactivate_nodes(nodes);

    // update stuff
    for (size_t i=0; i<faces.size(); i++)
        faces[i]->sigma = compute_sigma(faces[i]);
}

void perform_breaking(Mesh& mesh) {
    // compute sigma and separation strength for complete mesh
	for (size_t i=0; i< mesh.faces.size(); i++)
		mesh.faces[i]->sigma = compute_sigma(mesh.faces[i]);

	MeshSubset subset;
	for (size_t i=0; i< mesh.nodes.size(); i++) {		
        mesh.nodes[i]->flag &= ~Node::FlagResolveMax;
		if (separation_strength(mesh.nodes[i], 0, false) > 1)
			subset.active_nodes.push_back(mesh.nodes[i]);			
	}
	if (subset.active_nodes.empty())
		return;

	vector<AccelStruct*> obs_accs = create_accel_structs(sim.obstacle_meshes, false);
    
    // fracture sub-stepping
	for (int num=0; num < ::magic.max_cracks && !subset.active_nodes.empty(); num++) {
        subset.set_flag(Node::FlagMayBreak);
        subset.grow(SUPPORT_RINGS);		
		subset.set_flag(Node::FlagResolveMax);
    
		// dynamic remeshing
        map<Node*, Plane> planes = nearest_obstacle_planes(subset.get_all_nodes(), obs_accs);
        dynamic_remesh(subset, planes);
        subset.update_support();
    
        // physics update
        local_physics_step(subset);
    
        // recompute separation strength
        SplitNode cur_split(0);
        priority_queue<SplitNode> split;
        for (size_t i=0; i<subset.active_nodes.size(); i++)
            if ((subset.active_nodes[i]->flag & Node::FlagMayBreak) &&
                separation_strength(subset.active_nodes[i], &cur_split, false) > 1)
                split.push(cur_split);

        // rewind
        for (size_t i=0; i<subset.active_nodes.size(); i++)
            subset.active_nodes[i]->x = subset.active_nodes[i]->x0;
        for (size_t i=0; i<subset.support_nodes.size(); i++)
            subset.support_nodes[i]->x = subset.support_nodes[i]->x0;

        //subset.debug();
        subset.clear_flag(Node::FlagMayBreak);
        subset.active_nodes.clear();
        
        // break biggest one, put the rest into the zone for next round
        while (!split.empty()) {
            cur_split = split.top();
            split.pop();
            if (subset.active_nodes.empty()) {
                break_node(cur_split, subset);
            } else
                subset.active_nodes.push_back(cur_split.node);
        }
    }

    destroy_accel_structs(obs_accs);
    
	compute_ms_data(mesh);
}

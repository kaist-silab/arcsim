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

#include "simulation.hpp"

#include "breaking.hpp"
#include "collision.hpp"
#include "dynamicremesh.hpp"
#include "geometry.hpp"
#include "magic.hpp"
#include "display.hpp"
#include "nearobs.hpp"
#include "physics.hpp"
#include "plasticity.hpp"
#include "popfilter.hpp"
#include "proximity.hpp"
#include "separate.hpp"
#include "strainlimiting.hpp"
#include "io.hpp"
#include <iostream>
#include <fstream>
using namespace std;

bool consistency_check = false;
bool single_step = false;

static const bool verbose = false;
static const int proximity = Simulation::Proximity,
                 physics = Simulation::Physics,
                 strainlimiting = Simulation::StrainLimiting,
                 collision = Simulation::Collision,
                 remeshing = Simulation::Remeshing,
                 separation = Simulation::Separation,
                 popfilter = Simulation::PopFilter,
                 plasticity = Simulation::Plasticity,
                 fracture = Simulation::Fracture;

void physics_step (Simulation &sim, const vector<Constraint*> &cons);
void plasticity_step (Simulation &sim);
void strainlimiting_step (Simulation &sim, const vector<Constraint*> &cons);
void strainzeroing_step (Simulation &sim);
void equilibration_step (Simulation &sim);
void collision_step (Simulation &sim);
void remeshing_step (Simulation &sim, bool initializing=false);

void validate_handles (const Simulation &sim);

static void consistency(const char* text) {
	if (consistency_check) {
		cout << "> " << text << " ";
		test_state(sim, "/tmp/c");
	}
    if (single_step) {
        cout << "> " << text << endl;
        wait_key();
    }
}

void prepare (Simulation &sim) {
	for (size_t c = 0; c < sim.cloths.size(); c++) {
		Cloth& cloth = sim.cloths[c];
	    cloth.mesh.parent = &sim.cloths[c];
	    for (size_t i=0; i<cloth.mesh.faces.size(); i++)
            cloth.mesh.faces[i]->material = cloth.materials[cloth.mesh.faces[i]->flag];
        for (size_t i=0; i<cloth.mesh.edges.size(); i++)
            cloth.mesh.edges[i]->theta_ideal = dihedral_angle<MS>(cloth.mesh.edges[i]);
        compute_ms_data(cloth.mesh);
	}
    sim.cloth_meshes.resize(sim.cloths.size());
    for (size_t c = 0; c < sim.cloths.size(); c++) {
        sim.cloth_meshes[c] = &sim.cloths[c].mesh;
        update_x0(*sim.cloth_meshes[c]);
    }
    sim.obstacle_meshes.resize(sim.obstacles.size());
    for (size_t o = 0; o < sim.obstacles.size(); o++) {
        sim.obstacle_meshes[o] = &sim.obstacles[o].get_mesh();
        update_x0(*sim.obstacle_meshes[o]);
    }
}

void relax_initial_state (Simulation &sim) {
	validate_handles(sim);
    if (::magic.preserve_creases)
        for (int c = 0; c < (int)sim.cloths.size(); c++)
            reset_plasticity(sim.cloths[c]);
    
    bool equilibrate = true;
    if (equilibrate) {
        equilibration_step(sim);
        remeshing_step(sim, true);
        equilibration_step(sim);
    } else {
        remeshing_step(sim, true);
        strainzeroing_step(sim);
        remeshing_step(sim, true);
        strainzeroing_step(sim);
    }
    if (::magic.preserve_creases)
        for (int c = 0; c < (int)sim.cloths.size(); c++)
            reset_plasticity(sim.cloths[c]);
    ::magic.preserve_creases = false;
    if (::magic.fixed_high_res_mesh)
        sim.enabled[remeshing] = false;
}

void validate_handles (const Simulation &sim) {
    for (int h = 0; h < (int)sim.handles.size(); h++) {
        vector<Node*> nodes = sim.handles[h]->get_nodes();
        for (int n = 0; n < (int)nodes.size(); n++) {
            if (!nodes[n]->preserve) {
                cout << "Constrained node " << nodes[n]->index << " will not be preserved by remeshing" << endl;
                abort();
            }
        }
    }
}

vector<Constraint*> get_constraints (Simulation &sim, bool include_proximity);
void delete_constraints (const vector<Constraint*> &cons);
void update_obstacles (Simulation &sim, bool update_positions=true);

void advance_step (Simulation &sim);
void add_jitter (Simulation &sim);

void advance_frame (Simulation &sim) {
    for (int s = 0; s < sim.frame_steps; s++)
        advance_step(sim);
}

void advance_step (Simulation &sim) {	
	cout << "Sim frame " << sim.frame << " [" << sim.step << "]" << endl;
    sim.time += sim.step_time;
    sim.step++;

    if (::magic.add_jitter && sim.frame % 10 == 0 && sim.frame < 100)
    	add_jitter(sim);
    // super-simple manual timestep adaptivity
    if (sim.time >= sim.passive_time && sim.enabled[fracture]) {
        sim.enabled[fracture] = false;
        sim.step_time *= 5;
        magic.collision_stiffness *= 0.4;
    }
    Annotation::list.clear();
    update_obstacles(sim, false);
    vector<Constraint*> cons = get_constraints(sim, true);
    consistency("init step");
    physics_step(sim, cons);
    //cout << "phys" << endl;wait_key();
    consistency("physics");
    plasticity_step(sim);
    consistency("plasticity");
    strainlimiting_step(sim, cons);
    consistency("strainlimit");
    collision_step(sim);
    consistency("collision");
    //cout << "coll" << endl;wait_key();
    if (sim.step % sim.frame_steps == 0) {
        remeshing_step(sim);
        consistency("remeshing");
    	sim.frame++;
    }
    //cout << "rem" << endl;wait_key();
    
    delete_constraints(cons);
}

vector<Constraint*> get_constraints (Simulation &sim, bool include_proximity) {
    vector<Constraint*> cons;
    for (int h = 0; h < (int)sim.handles.size(); h++)
        append(cons, sim.handles[h]->get_constraints(sim.time));
    if (include_proximity && sim.enabled[proximity]) {
        sim.timers[proximity].tick();
        append(cons, proximity_constraints(sim.cloth_meshes,
                                           sim.obstacle_meshes,
                                           sim.friction, sim.obs_friction));
        sim.timers[proximity].tock();
    }
    return cons;
}

void delete_constraints (const vector<Constraint*> &cons) {
    for (int c = 0; c < (int)cons.size(); c++)
        delete cons[c];
}

// Steps

void update_velocities (vector<Mesh*> &meshes, vector<Vec3> &xold, double dt);

void step_mesh (Mesh &mesh, double dt);

void physics_step (Simulation &sim, const vector<Constraint*> &cons) {
    if (!sim.enabled[physics])
        return;
    sim.timers[physics].tick();
    for (size_t c = 0; c < sim.cloths.size(); c++) {
    	Mesh& mesh = sim.cloths[c].mesh;
        int nn = mesh.nodes.size();
        vector<Vec3> fext(nn, Vec3(0));
        vector<Mat3x3> Jext(nn, Mat3x3(0));

        activate_nodes(mesh.nodes);
        add_external_forces(mesh.nodes, mesh.faces, sim.gravity, sim.wind, fext, Jext);
        for (size_t h = 0; h < sim.handles.size(); h++)
            sim.handles[h]->add_forces(sim.time,fext,Jext);

        for (size_t m = 0; m < sim.morphs.size(); m++)
            if (sim.morphs[m].mesh == &sim.cloths[c].mesh)
                add_morph_forces(sim.cloths[c], sim.morphs[m], sim.time,
                                 sim.step_time, fext, Jext);
        vector<Vec3> dv = implicit_update(mesh.nodes, mesh.edges, mesh.faces, 
                                          fext, Jext, cons, sim.step_time);
        for (size_t n = 0; n < mesh.nodes.size(); n++) {
            mesh.nodes[n]->v += dv[n];
            mesh.nodes[n]->acceleration = dv[n] / sim.step_time;
        }
        //project_outside(mesh.nodes, cons);
        deactivate_nodes(mesh.nodes);
    }
    consistency("physics-pre");
    for (size_t c = 0; c < sim.cloth_meshes.size(); c++) {
        step_mesh(*sim.cloth_meshes[c], sim.step_time);
        compute_ws_data(*sim.cloth_meshes[c]);
    }
    for (size_t o = 0; o < sim.obstacle_meshes.size(); o++)
        step_mesh(*sim.obstacle_meshes[o], sim.step_time);
    sim.timers[physics].tock();
}

void step_mesh (Mesh &mesh, double dt) {
    for (int n = 0; n < (int)mesh.nodes.size(); n++)
        mesh.nodes[n]->x += mesh.nodes[n]->v*dt;
}

void plasticity_step (Simulation &sim) {
    if (!sim.enabled[plasticity])
        return;
    sim.timers[plasticity].tick();
    for (int c = 0; c < (int)sim.cloths.size(); c++) {
        plastic_update(sim.cloths[c]);
        optimize_plastic_embedding(sim.cloths[c]);
    }
    sim.timers[plasticity].tock();
}

vector<Vec3> node_positions (const vector<Mesh*> &meshes);

void strainlimiting_step (Simulation &sim, const vector<Constraint*> &cons) {
    if (!sim.enabled[strainlimiting])
        return;
    sim.timers[strainlimiting].tick();
    vector<Vec3> xold = node_positions(sim.cloth_meshes);
    strain_limiting(sim.cloth_meshes, get_strain_limits(sim.cloths), cons);
    update_velocities(sim.cloth_meshes, xold, sim.step_time);
    sim.timers[strainlimiting].tock();
}

void equilibration_step (Simulation &sim) {
    sim.timers[remeshing].tick();
    vector<Constraint*> cons;// = get_constraints(sim, true);
    // double stiff = 1;
    // swap(stiff, ::magic.handle_stiffness);
    for (int c = 0; c < (int)sim.cloths.size(); c++) {
        Mesh &mesh = sim.cloths[c].mesh;
        for (int n = 0; n < (int)mesh.nodes.size(); n++)
            mesh.nodes[n]->acceleration = Vec3(0);
        apply_pop_filter(sim.cloths[c], cons, 1);
    }
    // swap(stiff, ::magic.handle_stiffness);
    sim.timers[remeshing].tock();
    delete_constraints(cons);
    cons = get_constraints(sim, false);
    if (sim.enabled[collision]) {
        sim.timers[collision].tick();
        collision_response(sim.cloth_meshes, cons, sim.obstacle_meshes);
        sim.timers[collision].tock();
    }
    
    delete_constraints(cons);
}

void strainzeroing_step (Simulation &sim) {
    sim.timers[strainlimiting].tick();
    vector<StrainLimit> strain_limits(count_elements<Face>(sim.cloth_meshes), StrainLimit(1,1));
    vector<Constraint*> cons =
        proximity_constraints(sim.cloth_meshes, sim.obstacle_meshes,
                              sim.friction, sim.obs_friction);
    strain_limiting(sim.cloth_meshes, strain_limits, cons);
    delete_constraints(cons);
    sim.timers[strainlimiting].tock();
    if (sim.enabled[collision]) {
        sim.timers[collision].tock();
        collision_response(sim.cloth_meshes, vector<Constraint*>(),
                           sim.obstacle_meshes);
        sim.timers[collision].tock();
    }
}

void collision_step (Simulation &sim) {
    if (!sim.enabled[collision])
        return;
    sim.timers[collision].tick();
    vector<Vec3> xold = node_positions(sim.cloth_meshes);
    vector<Constraint*> cons = get_constraints(sim, false);
    collision_response(sim.cloth_meshes, cons, sim.obstacle_meshes);
    delete_constraints(cons);
    update_velocities(sim.cloth_meshes, xold, sim.step_time);
    sim.timers[collision].tock();
}

void remeshing_step (Simulation &sim, bool initializing) {
    if (!sim.enabled[remeshing])
        return;
    
    // remesh
    sim.timers[remeshing].tick();
    for (size_t c = 0; c < sim.cloths.size(); c++) {
        if (::magic.fixed_high_res_mesh)
            static_remesh(sim.cloths[c].mesh);
        else {
        	vector<AccelStruct*> obs_accs = create_accel_structs(sim.obstacle_meshes, false);
            map<Node*,Plane> planes = nearest_obstacle_planes(sim.cloths[c].mesh.nodes, obs_accs);
            destroy_accel_structs(obs_accs);
            dynamic_remesh(sim.cloths[c].mesh, planes);
        }
    }
    sim.timers[remeshing].tock();
    consistency("post mesh");

    // breaking
    if (sim.enabled[fracture] && sim.frame > 1) {
    	for (size_t c = 0; c < sim.cloths.size(); c++) {
    		perform_breaking(sim.cloths[c].mesh);
    	}
    	consistency("breaking");
    }
    
    // separate
    if (sim.enabled[separation]) {
        sim.timers[separation].tick();
        separate(sim.cloth_meshes, sim.obstacle_meshes);
        sim.timers[separation].tock();
    }
    consistency("separation");

    // apply pop filter
    if (sim.enabled[popfilter] && !initializing) {
        sim.timers[popfilter].tick();
        vector<Constraint*> cons = get_constraints(sim, true);
        for (size_t c = 0; c < sim.cloths.size(); c++)
            apply_pop_filter(sim.cloths[c], cons);
        delete_constraints(cons);
        sim.timers[popfilter].tock();
    }    
}

void update_velocities (vector<Mesh*> &meshes, vector<Vec3> &xold, double dt) {
    double inv_dt = 1/dt;
    int idx=0;
    for (size_t m = 0; m < meshes.size(); m++) {
    	vector<Node*>& nodes = meshes[m]->nodes;
    	for (size_t n = 0; n < nodes.size(); n++)
        	nodes[n]->v += (nodes[n]->x - xold[idx++]) * inv_dt;
    }
}

void update_obstacles (Simulation &sim, bool update_positions) {
    double decay_time = 0.1,
           blend = sim.step_time/decay_time;
    blend = blend/(1 + blend);
    for (int o = 0; o < (int)sim.obstacles.size(); o++) {
        sim.obstacles[o].get_mesh(sim.time);
        sim.obstacles[o].blend_with_previous(sim.time, sim.step_time, blend);
        if (!update_positions) {
            // put positions back where they were
            Mesh &mesh = sim.obstacles[o].get_mesh();
            for (int n = 0; n < (int)mesh.nodes.size(); n++) {
                Node *node = mesh.nodes[n];
                node->v = (node->x - node->x0)/sim.step_time;
                node->x = node->x0;
            }
        }
    }
}

void add_jitter (Simulation& sim) { 
	for (size_t c = 0; c<sim.cloth_meshes.size(); c++)
		for (size_t i = 0; i<sim.cloth_meshes[c]->nodes.size(); i++) 
			sim.cloth_meshes[c]->nodes[i]->x += sim.cloth_meshes[c]->nodes[i]->n * 
		                                        (1e-3*(rand() % 1000) * 1e-4);
}

// Helper functions

vector<Vec3> node_positions (const vector<Mesh*> &meshes) {
    vector<Vec3> xs(count_elements<Node>(meshes));
    int idx = 0;
    for (size_t m = 0; m < meshes.size(); m++) {
    	vector<Node*>& nodes = meshes[m]->nodes;
    	for (size_t n = 0; n < nodes.size(); n++)
        	xs[idx++] = nodes[n]->x;
    }
    return xs;
}

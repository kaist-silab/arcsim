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

#ifndef PHYSICS_HPP
#define PHYSICS_HPP

#include "cloth.hpp"
#include "geometry.hpp"
#include "simulation.hpp"
#include <vector>

template <Space s> Mat3x3 deformation_gradient (const Face *face);

Mat3x3 material_model (const Face *face, const Mat3x3& G);

template <Space s> 
double internal_energy (const std::vector<Face*>& faces, const std::vector<Edge*>& edges);

double constraint_energy (const std::vector<Constraint*> &cons);

double external_energy (const Cloth &cloth, const Vec3 &gravity,
                        const Wind &wind);

// A += dt^2 dF/dx; b += dt F + dt^2 dF/dx v
// also adds damping terms
// if dt == 0, just does A += dF/dx; b += F instead, no damping
template <Space s>
void add_internal_forces (const std::vector<Face*>& faces, const std::vector<Edge*>& edges,
						  SpMat<Mat3x3> &A, std::vector<Vec3> &b, double dt);

void add_constraint_forces (const std::vector<Constraint*> &cons,
                            SpMat<Mat3x3> &A, std::vector<Vec3> &b, double dt);

void add_external_forces (const std::vector<Node*>& nodes, const std::vector<Face*>& faces, 
						  const Vec3 &gravity, const Wind &wind, std::vector<Vec3> &fext,
                          std::vector<Mat3x3> &Jext);

void add_morph_forces (const Cloth &cloth, const Morph &morph, double t,
                       double dt,
                       std::vector<Vec3> &fext, std::vector<Mat3x3> &Jext);

std::vector<Vec3> implicit_update (std::vector<Node*>& nodes, const std::vector<Edge*>& edges, 
					  const std::vector<Face*>& faces,
					  const std::vector<Vec3> &fext, const std::vector<Mat3x3> &Jext,
                      const std::vector<Constraint*> &cons, double dt);

void project_outside (std::vector<Node*>& nodes, const std::vector<Constraint*> &cons);

#endif

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

#include "io.hpp"
#include "util.hpp"
#include <fstream>
#include <map>

using namespace std;

void tri2obj (const vector<string> &args) {
    if (args.size() != 2) {
        cout << "Converts output of Jonathan Shewchuk's Triangle program "
             << "to an OBJ file." << endl;
        cout << "Arguments:" << endl;
        cout << "    <tri>: Common prefix of <tri>.node and <tri>.ele"
             << endl;
        cout << "    <obj>: Output OBJ file" << endl;
        exit(EXIT_FAILURE);
    }
    assert(args.size() == 2);
    triangle_to_obj(args[0], args[1]);
}

void merge_meshes (const vector<string> &args) {
    assert(args.size() == 3 || args.size() == 4);
    Mesh meshm, meshw;
    load_obj(meshm, args[0]);
    load_obj(meshw, args[1]);
    assert(meshm.nodes.size() == meshw.nodes.size());
    Mesh mesh;
    double merge_dist = args.size()==4 ? atof(args[3].c_str()) : 1e-2;
    for (int n = 0; n < (int)meshm.nodes.size(); n++) {
        Vert *vert = new Vert(meshm.nodes[n]->x);
        mesh.add(vert);
        const Vec3 x = meshw.nodes[n]->x;
        Node *node = NULL;
        for (int n2 = 0; n2 < (int)mesh.nodes.size(); n2++) {
            Node *node2 = mesh.nodes[n2];
            if (norm2(node2->x - x) < sq(merge_dist)) {
                node = node2;
                break;
            }
        }
        if (node == NULL) {
            node = new Node(x, x, Vec3(0), 0, 0, false);
            mesh.add(node);
        }
        if (!node->label)
            node->label = meshw.nodes[n]->label;
        connect(vert, node);
    }
    for (int f = 0; f < (int)meshm.faces.size(); f++) {
        const Face *face = meshm.faces[f];
        mesh.add(new Face(mesh.verts[face->v[0]->node->index],
                          mesh.verts[face->v[1]->node->index],
                          mesh.verts[face->v[2]->node->index], 
                          Mat3x3(1), Mat3x3(0), face->material, face->damage));
    }
    save_obj(mesh, args[2]);
}

void split_meshes (const vector<string> &args) {
    assert(args.size() == 3);
    Mesh mesh;
    load_obj(mesh, args[0]);
    Mesh meshm, meshw;
    for (int v = 0; v < (int)mesh.verts.size(); v++) {
        const Vert *vert = mesh.verts[v];
        meshm.add(new Vert(vert->u));
        meshm.add(new Node(vert->u, vert->u, Vec3(0), 0, 0, false));
        connect(meshm.verts.back(), meshm.nodes.back());
        meshw.add(new Vert(vert->u));
        meshw.add(new Node(vert->node->x, vert->node->x, Vec3(0), 
        	               vert->node->label, vert->node->flag, false));
        connect(meshw.verts.back(), meshw.nodes.back());
    }
    for (int f = 0; f < (int)mesh.faces.size(); f++) {
        const Face *face = mesh.faces[f];
        meshm.add(new Face(meshm.verts[face->v[0]->index],
                           meshm.verts[face->v[1]->index],
                           meshm.verts[face->v[2]->index], 
                           Mat3x3(1), Mat3x3(0), face->material, face->damage));
        meshw.add(new Face(meshw.verts[face->v[0]->index],
                           meshw.verts[face->v[1]->index],
                           meshw.verts[face->v[2]->index], 
                           Mat3x3(1), Mat3x3(0), face->material, face->damage));
    }
    save_obj(meshm, args[1]);
    save_obj(meshw, args[2]);
}

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

#ifndef DISPLAY_HPP
#define DISPLAY_HPP

#include "simulation.hpp"
#include "timer.hpp"
#include <vector>

extern Timer fps;
extern int frame;

struct GlutCallbacks {
    void (*idle) ();
    void (*keyboard) (unsigned char, int, int);
    void (*special) (int, int, int);
    GlutCallbacks (): idle(NULL), keyboard(NULL), special(NULL) {}
};

struct Annotation {
	Vec3 color;
	Vec3 pos, dir;
	Face* face;
	Edge* edge;
	Node* node;

	static void add(Face* f, Vec3 c = Vec3(1,0,0)) { list.push_back(Annotation(f,0,0,c,Vec3(0))); }
	static void add(Edge* e, Vec3 c = Vec3(1,0,0)) { list.push_back(Annotation(0,e,0,c,Vec3(0))); }
	static void add(Node* n, Vec3 c = Vec3(1,0,0)) { list.push_back(Annotation(0,0,n,c,Vec3(0))); }
	static void add(Vec3 pos, Vec3 c = Vec3(1,0,0)) { list.push_back(Annotation(0,0,0,c,pos)); }
	static void add(Vec3 pos, Vec3 dir, Vec3 c = Vec3(1,0,0)) { list.push_back(Annotation(0,0,0,c,pos,dir)); }
	static std::vector<Annotation> list;
private:
	Annotation(Face* f, Edge* e, Node* n, Vec3 c, Vec3 p, Vec3 d=Vec3(0)) : 
		color(c),pos(p),dir(d),face(f),edge(e),node(n) {}
};

struct Pane {
	double lat, lon;
    Vec2 offset;
    Vec3 center;
    double scale;
    int window, parent;
    bool enabled;
    Vec3 pos(Vert *v);
    bool initialized;

    static Pane* current();
    static Pane panes[3];
    static Pane& material() { return panes[0]; }
	static Pane& world() { return panes[2]; }
	static Pane& plastic() { return panes[1]; }

    Pane (bool enable): lat(0), lon(0), offset(0), scale(0.5), enabled(enable), initialized(false) {}
};

void init_glut (const GlutCallbacks&);
void run_glut ();

void redisplay ();
void wait_key ();

#endif

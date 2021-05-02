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

#include "display.hpp"

#ifndef NO_OPENGL

#include "bvh.hpp"
#include "geometry.hpp"
#include "io.hpp"
#include "timer.hpp"
#include "opengl.hpp"
#include "physics.hpp"
#include "timer.hpp"
#include "dynamicremesh.hpp"
#include "util.hpp"
#include "sepstrength.hpp"
#include <sstream>

using namespace std;
string obj2png_filename;
extern string outprefix;

extern int frame;
extern Timer fps;

bool stepDebug;

vector<Annotation> Annotation::list;
Pane Pane::panes[3] = { false, false, true };

int display_mode = 0;
struct DisplayMode { string name; double scale; bool active; };
DisplayMode display_modes[] = { {"sigma", 1e4, true},
                                {"sigma_bend", 1e4, true},
                                {"sep strength", 1, true},
                                {"Sp_str", 1, true},
                                {"Sp_bend", 1, true},
                                {"-sentinel-", 0, false} };

Pane* Pane::current() { 
	int id = glutGetWindow();
	for (int i=0; i<3; i++)
		if (panes[i].window == id)
			return &panes[i];
	return 0;
}

Vec3 Pane::pos(Vert* v) {
	if (this == &world()) return v->node->x;
	if (this == &plastic()) return v->node->y;
	return v->u;
}

void reshape (int w, int h) {
    int npanes = 0;
    for (int i = 0; i < 3; i++)
        if (Pane::panes[i].enabled)
            npanes++;
    int j = 0;
    for (int i = 0; i < 3; i++) {
        glutSetWindow(Pane::panes[i].window);
        int x0 = w*j/npanes, x1 = Pane::panes[i].enabled ? w*(j+1)/npanes : x0+1;
        glutPositionWindow(x0, 0);
        glutReshapeWindow(x1-x0, h);
        glViewport(0, 0, x1-x0, h);
        if (Pane::panes[i].enabled)
            j = j + 1;
    }
}


void vertex (const Vec2 &x) {
    glVertex2d(x[0], x[1]);
}

void vertex (const Vec3 &x) {
    glVertex3d(x[0], x[1], x[2]);
}

void normal (const Vec3 &n) {
    glNormal3d(n[0], n[1], n[2]);
}

void color (const Vec3 &x) {
    glColor3d(x[0], x[1], x[2]);
}

Vec3 strain_color (const Face *face) {
    Mat3x3 F = deformation_gradient<WS>(face);
    Vec3 l = eigen_values(F.t()*F);
    double s0 = sqrt(l[0]) - 1, s1 = sqrt(l[1]) - 1;
    double tens = clamp(1e2*s0, 0., 0.5), comp = clamp(-1e2*s1, 0., 0.5);
    return Vec3(1-tens, (1-tens)*(1-comp), (1-comp));
}

Vec3 plasticity_color (const Face *face) {
    double s = norm_F(face->Sp_bend)/1000;
    double d = face->damage;
    s = min(s/2, 0.5);
    d = min(d/2, 0.5);
    return Vec3(1-s, (1-s)*(1-d), (1-d));
    // return colormap(trace(face->S_plastic)/500);
}

Vec3 origami_color (const Mat3x3& M) {
    double H = trace(M);
    return 0.9*Vec3(1 + H, 1 - abs(H)/2, 1 - H);
}

inline double matrix_mag(const Mat3x3& M) {
	Vec3 l = eigen_values(M);    
    return sqrt(sq(l[0])+sq(l[1])+sq(l[2])) * sgn(l[0]+l[1]+l[2]);
}

// http://www.sron.nl/~pault/colourschemes.pdf
// [-1,1]
Vec3 red_blue_colorscheme(double v) {
    v=clamp(0.5*(v+1.0),0.0,1.0);
    double v2=sq(v),v3=v2*v,v4=v3*v,v5=v4*v;
    return Vec3 (0.237-2.13*v + 26.92*v2-65.5*v3+63.5*v4-22.36*v5, 
                 sq((0.572+1.524*v-1.811*v2)/(1-0.291*v+0.1574*v2)),
                 1.0/(1.579-4.03*v+12.92*v2-31.4*v3+48.6*v4-23.36*v5));
}

Vec3 debug_color (Face *face, const Vert* vert) {
    double sc = 1.0/display_modes[display_mode].scale;
    switch(display_mode) {
        case 0: { // sigma        	
            compute_ms_data(face);
            Mat3x3 F = deformation_gradient<WS>(face);
    		Mat3x3 G = (F.t()*F - Mat3x3(1)) * 0.5;
    		Mat3x3 sigma = material_model(face, G);
			return red_blue_colorscheme(sc * matrix_mag(sigma));}
		case 1: { // sigma_bend
            Mat3x3 F_bend = Mat3x3(1);
            for (int i=0; i<3; i++) 
                F_bend += face->v[i]->node->curvature * (0.5/3.0 * face->material->fracture_bend_thickness);
            Mat3x3 G_bend = (F_bend.t()*F_bend - Mat3x3(1)) * 0.5;
            Mat3x3 sigma_bend = material_model(face, G_bend);
            return red_blue_colorscheme(sc * matrix_mag(sigma_bend));
        }
        case 2: { // sep strength
			//return red_blue_colorscheme(sc * separation_strength(vert->node,0,false));
            return red_blue_colorscheme(sc * vert->node->sep);
		}
        case 3: {// Sp_str
            double frob = 0;
            Mat3x3 S = face->Sp_str - Mat3x3(1);
            frob += norm2(S.col(0)) + norm2(S.col(1)) + norm2(S.col(2));
            return red_blue_colorscheme(sc * frob);
        }
        case 4: {// Sp_bend
            /*double frob = 0;
            Mat3x3 S = face->Sp_bend;
            frob += norm2(S.col(0)) + norm2(S.col(1)) + norm2(S.col(2));
            return red_blue_colorscheme(sc * frob);*/
            double v = abs(face->adje[0]->theta_ideal) +
                       abs(face->adje[1]->theta_ideal) +
                       abs(face->adje[2]->theta_ideal);
            return red_blue_colorscheme(sc * v);
        }
        default:
        	return Vec3(0);
    }
}

void draw_mesh_ms (Mesh &mesh, bool set_color=false) {
    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < mesh.faces.size(); i++) {
        if (i % 256 == 0) {
            glEnd();
            glBegin(GL_TRIANGLES);
        }
        Face *face = mesh.faces[i];
        for (int v = 0; v < 3; v++) {
        	if (set_color)
        	    color(debug_color(face,face->v[v]));        
            vertex(face->v[v]->u);
        }
    }
    glEnd();
    glLineWidth(2);
    glColor3d(0,0,0);
    glBegin(GL_LINES);
    for (size_t i=0; i<mesh.edges.size(); i++) {
        if (is_seam_or_boundary(mesh.edges[i])) {
            vertex(mesh.edges[i]->n[0]->verts[0]->u);
            vertex(mesh.edges[i]->n[1]->verts[0]->u);
        }
    }
    glEnd();
    glLineWidth(1);
    /*glPointSize(5);
    color(Vec3(1,0,0));
    glBegin(GL_POINTS);
    for (size_t i = 0; i < mesh.verts.size(); i++) {
    	if (mesh.verts[i]->node->preserve)
    	vertex(mesh.verts[i]->u);
    }   
    glEnd(); */
}

void draw_meshes_ms (bool set_color=false) {
    for (size_t m = 0; m < sim.cloth_meshes.size(); m++)
        draw_mesh_ms(*sim.cloth_meshes[m], set_color);
}

void shrink_face (const Face *face, double shrink_factor, double shrink_max,
                  Vec3 u[3]) {
    Vec3 u0 = face->v[0]->u, u1 = face->v[1]->u, u2 = face->v[2]->u;
    double a = face->a;
    double l = max(norm(u0 - u1), max(norm(u1 - u2), norm(u2 - u0)));
    double h = 2*a/l;
    double dh = min(h*shrink_factor, shrink_max);
    for (int v = 0; v < 3; v++) {
        Vec3 e1 = normalize(face->v[NEXT(v)]->u - face->v[v]->u),
             e2 = normalize(face->v[PREV(v)]->u - face->v[v]->u);
        Vec3 du = (e1 + e2)*dh / norm(cross(e1, e2));
        u[v] = face->v[v]->u + du;
    }
}

void draw_meshes_ms_fancy () {
    double shrink_factor = 0.1, shrink_max = 0.5e-3;
    for (size_t m = 0; m < sim.cloth_meshes.size(); m++) {
        const Mesh &mesh = *sim.cloth_meshes[m];
        glBegin(GL_TRIANGLES);
        glColor3f(0.5,0.5,0.5);
        for (int f = 0; f < (int)mesh.faces.size(); f++) {
            const Face *face = mesh.faces[f];
            for (int v = 0; v < 3; v++)
                vertex(face->v[v]->u);
        }
        glEnd();
        glBegin(GL_TRIANGLES);
        for (int f = 0; f < (int)mesh.faces.size(); f++) {
            const Face *face = mesh.faces[f];
            color(origami_color(face->Sp_bend)/1000.0);
            glColor3f(0.9,0.9,0.9);
            Vec3 u[3];
            shrink_face(face, shrink_factor, shrink_max, u);
            for (int v = 0; v < 3; v++)
                vertex(u[v]);
        }
        glEnd();
    }
}

void draw_mesh_ps (const Mesh &mesh, bool set_color=false) {
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < (int)mesh.faces.size(); i++) {
        Face *face = mesh.faces[i];
        if (i % 256 == 0) {
            glEnd();
            glBegin(GL_TRIANGLES);
        }
        normal(normal<PS>(face));
        for (int v = 0; v < 3; v++)
            vertex(face->v[v]->node->y);
    }
    glEnd();
}

void draw_meshes_ps (bool set_color=false) {
    for (size_t m = 0; m < sim.cloth_meshes.size(); m++)
        draw_mesh_ps(*sim.cloth_meshes[m], set_color);
}

template <Space s>
void draw_annotation(Annotation& a) {
	// glDisable(GL_COLOR_MATERIAL);
	glColor3f(a.color[0],a.color[1],a.color[2]);
	glLineWidth(2);
	glPointSize(6);
	if (a.face) {
		glBegin(GL_TRIANGLES);
		for (int i=0; i<3; i++)
			vertex(pos<s>(a.face->v[i]));
		glEnd();
	} else if (a.edge) {
		glBegin(GL_LINES);
		vertex(pos<s>(a.edge->n[0]->verts[0]));
		vertex(pos<s>(a.edge->n[1]->verts[0]));
		glEnd();
	} else if (a.node) {
		glBegin(GL_POINTS);		
		vertex(pos<s>(a.node->verts[0]));
		glEnd();
	} else {
		if (norm(a.dir) > 0) {
			glBegin(GL_LINES);
			vertex(a.pos);
			vertex(a.pos + a.dir);
			glEnd();
		} else {	
			glBegin(GL_POINTS);	
			vertex(a.pos);
			glEnd();
		}
	}
	glLineWidth(1);
    glPointSize(1);
	// glEnable(GL_COLOR_MATERIAL);
}

template <Space s>
void draw_mesh (const Mesh &mesh, bool set_color=false) {
    if (set_color)
        glDisable(GL_COLOR_MATERIAL);
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < (int)mesh.faces.size(); i++) {
        Face *face = mesh.faces[i];
        if (i % 256 == 0) {
            glEnd();
            glBegin(GL_TRIANGLES);
        }
        if (set_color) {
            int c = find((Mesh*)&mesh, sim.cloth_meshes);
            static const float phi = (1+sqrt(5))/2;
            double hue = c*(2 - phi)*2*M_PI; // golden angle
            hue = -0.6*M_PI + hue; // looks better this way :/
            //if (face->label % 2 == 1) hue += M_PI;
            static Vec3 a = Vec3(0.92, -0.39, 0), b = Vec3(0.05, 0.12, -0.99);
            Vec3 frt = Vec3(0.7,0.7,0.7) + (a*cos(hue) + b*sin(hue))*0.3,
                 bak = frt*0.5 + Vec3(0.5,0.5,0.5);
            float front[4] = {(float)frt[0], (float)frt[1], (float)frt[2], 1},
                  back[4] = {(float)bak[0], (float)bak[1], (float)bak[2], 1};
            glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, front);
            glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, back);
            // color(area_color(face));
        }
        normal(normal<s>(face));
        for (int v = 0; v < 3; v++)
            vertex(pos<s>(face->v[v]));
    }
    glEnd();
    if (set_color)
        glEnable(GL_COLOR_MATERIAL);
}

template <Space s>
void draw_meshes (bool set_color=false) {
    for (size_t m = 0; m < sim.cloth_meshes.size(); m++)
        draw_mesh<s>(*sim.cloth_meshes[m], set_color);
}

template <Space s>
void draw_seam_or_boundary_edges () {
    glColor3f(0,0,0);
    glBegin(GL_LINES);
    for (size_t m = 0; m < sim.cloth_meshes.size(); m++) {
        const Mesh &mesh = *sim.cloth_meshes[m];
        for (int e = 0; e < (int)mesh.edges.size(); e++) {
            const Edge *edge = mesh.edges[e];
            if (!is_seam_or_boundary(edge))
                continue;
            vertex(pos<s>(edge->n[0]));
            vertex(pos<s>(edge->n[1]));
        }
    }
    glEnd();
}

void draw_node_vels () {
    double dt = 0.01;
    glBegin(GL_LINES);
    for (size_t m = 0; m < sim.cloth_meshes.size(); m++) {
        const Mesh &mesh = *sim.cloth_meshes[m];
        for (int n = 0; n < (int)mesh.nodes.size(); n++) {
            const Node *node = mesh.nodes[n];
            glColor3d(0,0,1);
            vertex(node->x);
            vertex(node->x + dt*node->v);
            glColor3d(1,0,0);
            vertex(node->x);
            vertex(node->x - dt*node->v);
        }
    }
    for (size_t o = 0; o < sim.obstacles.size(); o++) {
        const Mesh &mesh = sim.obstacles[o].get_mesh();
        for (int n = 0; n < (int)mesh.nodes.size(); n++) {
            const Node *node = mesh.nodes[n];
            glColor3d(0,0,1);
            vertex(node->x);
            vertex(node->x + dt*node->v);
            glColor3d(1,0,0);
            vertex(node->x);
            vertex(node->x - dt*node->v);
        }
    }
    glEnd();
}

void draw_node_accels () {
    double dt2 = 1e-6;
    glBegin(GL_LINES);
    for (size_t m = 0; m < sim.cloth_meshes.size(); m++) {
        const Mesh &mesh = *sim.cloth_meshes[m];
        for (int n = 0; n < (int)mesh.nodes.size(); n++) {
            const Node *node = mesh.nodes[n];
            glColor3d(0,0,1);
            vertex(node->x);
            vertex(node->x + dt2*node->acceleration);
            glColor3d(1,0,0);
            vertex(node->x);
            vertex(node->x - dt2*node->acceleration);
        }
    }
    glEnd();
}

void directional_light (int i, const Vec3 &dir, const Vec3 &dif) {
    float diffuse[4] = {(float)dif[0], (float)dif[1], (float)dif[2], 1};
    float position[4] = {(float)dir[0], (float)dir[1], (float)dir[2], 0};
    glEnable(GL_LIGHT0+i);
    glLightfv(GL_LIGHT0+i, GL_DIFFUSE, diffuse);
    glLightfv(GL_LIGHT0+i, GL_POSITION, position);
}

void ambient_light (const Vec3 &a) {
    float ambient[4] = {(float)a[0], (float)a[1], (float)a[2], (float)1};
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
}

double aspect_ratio () {
    return (double)glutGet(GLUT_WINDOW_WIDTH)/glutGet(GLUT_WINDOW_HEIGHT);
}

void apply_view (const Pane &view) {
    glTranslatef(view.offset[0], view.offset[1], 0);
    glScalef(view.scale, view.scale, view.scale);
    glRotatef(view.lat, 1,0,0);
    glRotatef(view.lon, 0,1,0);    
    glTranslatef(view.center[0], view.center[1], view.center[2]);    
}

template<Space s>
void init_view(Pane& view) {
    if (view.initialized) 
        return;
    
    Vec3 c(0);
    int n=0;    
    for (size_t m=0; m<sim.cloth_meshes.size(); m++) {
        const vector<Vert*>& verts = sim.cloth_meshes[m]->verts;
        for (size_t i=0; i<verts.size(); i++) {
            c += pos<s>(verts[i]);
            n++;
        }            
    }
    view.center = Vec3(0);//-c / ((double)n);
    view.initialized = true;
}

void display_material () {
    init_view<MS>(Pane::material());
    glClearColor(1,1,1,1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1,1);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45, aspect_ratio(), 0.001, 1000);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0, 0, -1);
    // draw_meshes_ms_fancy();
    glColor3d(0.9,0.9,0.9);
    apply_view(Pane::material());
    
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glEnable(GL_CULL_FACE);
    draw_meshes_ms(true);
    glDisable(GL_CULL_FACE);
    glColor4d(0,0,0, 0.2);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    draw_meshes_ms();
    
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    for (size_t i=0; i<Annotation::list.size(); i++)
    	draw_annotation<MS>(Annotation::list[i]);
    // draw_ms_edges();
    glutSwapBuffers();
    GLenum errCode;
    const GLubyte *errString;
    if ((errCode = glGetError()) != GL_NO_ERROR) {
        errString = gluErrorString(errCode);
        cout << "OpenGL Error (ms): " << errString << endl;
    }
}

void display_plastic () {
    if (!Pane::plastic().enabled)
        return;
    init_view<PS>(Pane::plastic());
    glClearColor(1,1,1,1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1,1);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45, aspect_ratio(), 0.001, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0, 0, -1);
    glEnable(GL_LIGHTING);
    glEnable(GL_NORMALIZE);
    directional_light(0, Vec3(0,0,1), Vec3(0.5,0.5,0.5));
    ambient_light(Vec3(0.5));
    apply_view(Pane::plastic());
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    draw_meshes<PS>(true);
    glColor4d(0,0,0, 0.2);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    draw_meshes<PS>();
    draw_seam_or_boundary_edges<PS>();
    glutSwapBuffers();
    GLenum errCode;
    const GLubyte *errString;
    if ((errCode = glGetError()) != GL_NO_ERROR) {
        errString = gluErrorString(errCode);
        cout << "OpenGL Error (ps): " << errString << endl;
    }
}

void display_world () {
    init_view<WS>(Pane::world());
    glClearColor(1,1,1,1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1,1);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45, aspect_ratio(), 0.001, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0, 0, -1);
    glEnable(GL_LIGHTING);
    glEnable(GL_NORMALIZE);
    directional_light(0, Vec3(0,0,1), Vec3(0.5,0.5,0.5));
    ambient_light(Vec3(0.5));
    apply_view(Pane::world());
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    draw_meshes<WS>(true);
    glEnable(GL_CULL_FACE);
    glColor3f(0.8,0.8,0.8);
    for(int o = 0; o < (int)sim.obstacles.size(); o++)
        draw_mesh<WS>(sim.obstacles[o].get_mesh());
    glDisable(GL_CULL_FACE);
    glColor4d(0,0,0, 0.2);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    draw_meshes<WS>();
    draw_seam_or_boundary_edges<WS>();

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    for (size_t i=0; i<Annotation::list.size(); i++)
    	draw_annotation<WS>(Annotation::list[i]);
    
    glutSwapBuffers();
    GLenum errCode;
    const GLubyte *errString;
    if ((errCode = glGetError()) != GL_NO_ERROR) {
        errString = gluErrorString(errCode);
        cout << "OpenGL Error (ws): " << errString << endl;
    }

    // Update Title
    glutSetWindow(Pane::world().parent);
    char title[1024];
    sprintf(title, "%s frame %d time %.3f  debug: %s [scale %g]",outprefix.c_str(), sim.frame, sim.time,
    	display_modes[display_mode].name.c_str(), display_modes[display_mode].scale);
    glutSetWindowTitle(title);
}

struct MouseState {
    bool down;
    int x, y;
    int down_x, down_y;
    enum {ROTATE, TRANSLATE, SCALE} func;
} mouse_state;

void zoom (bool in) {
	Pane* pane = Pane::current();
	if (!pane) return;

    if (in)
        pane->scale *= 1.2;
    else
        pane->scale /= 1.2;
    glutPostRedisplay();
}

// click on element for debug information
void select_element(int x, int y, int button) {
	Pane* pane = Pane::current();    
	double modelview[16], project[16];
    int viewport[4];
    Vec3 p0,p1;
    
    glutSetWindow(pane->window);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
    glGetDoublev(GL_PROJECTION_MATRIX, project);
    glGetIntegerv(GL_VIEWPORT, viewport);
    // get 3d ray
    gluUnProject(x, viewport[3]-y, 0, modelview, project, viewport, &p0[0], &p0[1], &p0[2]);
    gluUnProject(x, viewport[3]-y, 1, modelview, project, viewport, &p1[0], &p1[1], &p1[2]);
    Face *face = 0; Vert *vert = 0;
    double maxZ = 1e100;
    double minD = 0.03;
    
    for (size_t c=0; c<sim.cloth_meshes.size(); c++) {
    	Mesh& mesh = *sim.cloth_meshes[c];
	    if (button == GLUT_LEFT_BUTTON) {
	        for (size_t t=0; t<mesh.faces.size(); t++) {
	            Face* cur = mesh.faces[t];            
	            double z;
	            if (triangle_ray_test(pane->pos(cur->v[0]), pane->pos(cur->v[1]), 
	                                  pane->pos(cur->v[2]), p0, p1-p0, z)) {
	                if (fabs(z)<maxZ) {
	                    maxZ = fabs(z);
	                    face = cur;
	                }
	            }
	        }
	    } else if(button == GLUT_RIGHT_BUTTON) {
	        for (size_t v=0; v<mesh.verts.size(); v++) {
	            Vert* cur = mesh.verts[v];
	            Vec3 ax = pane->pos(cur);
	            double d = norm(cross(ax-p0,ax-p1))/norm(p1-p0);
	            if (d < minD) {
	                minD = d;
	                vert = cur;
	            }
	        }
	    } 
	}
	Annotation::list.clear();
	if (face) {
		Annotation::add(face);
        cout << face << " index " << face->index << endl;
        map<Node*,Plane> planes;
        compute_face_sizing(face->v[0]->node->mesh->parent->remeshing,face, planes, true);
        Vec3 eig = eigen_values(face->sigma);
        cout << "max sigma " << eig[0] << " (toughness " << face->material->toughness << ")" << endl;        
        cout << "aspect " << aspect(face) << endl;
        cout << "Sp_str " << face->Sp_str << endl;
        cout << "theta ideal " << face->adje[0]->theta_ideal << " " << face->adje[1]->theta_ideal << " " << face->adje[2]->theta_ideal << endl;
	}
	if (vert) {
		Annotation::add(vert->node);
		cout << vert << " index " << vert->index << endl;
		cout << "preserve " << vert->node->preserve << " flags " << vert->node->flag << " label " 
             << vert->node->label << endl;
        cout << "sep " << separation_strength(vert->node,0,true) << endl;
        cout << "sizing " << vert->sizing << endl;
        cout << "x " << vert->node->x << endl;
        cout << "u " << vert->u << endl;
        cout << "vel " << vert->node->v << endl;
        cout << "acc " << vert->node->acceleration << endl;
        for (size_t i=0; i<vert->adjf.size(); i++)
            cout << vert->adjf[i] << " ";
        cout << endl;
	}
	redisplay();
}

void mouse (int button, int state, int x, int y) {
    mouse_state.down = (state == GLUT_DOWN);
    mouse_state.x = x;
    mouse_state.y = y;

    if (state == GLUT_UP && abs(mouse_state.down_x-x) < 2 && abs(mouse_state.down_y-y) < 2 &&
    	(button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON))
    	select_element(x, y, button);

    mouse_state.down_x = x;
	mouse_state.down_y = y;

    Pane* pane = Pane::current();
    if (!pane) return;
    
    if ((button == 3) || (button == 4)) { // It's a wheel event
        mouse_state.func = MouseState::SCALE;
        // Each wheel event reports like a button click, GLUT_DOWN then GLUT_UP
        if (state == GLUT_UP) return; // Disregard redundant GLUT_UP events
        
        double oldScale = pane->scale;        
        if (button == 3) {
            pane->scale *= 1.02;
        } else {
            pane->scale /= 1.02;
        }
        double s = pane->scale/oldScale;
        Vec2 mvVec = Vec2((double)x/glutGet(GLUT_WINDOW_WIDTH)-0.5, 
                          (double)-y/glutGet(GLUT_WINDOW_HEIGHT)+0.5);
        pane->offset = s * pane->offset + (1-s)*aspect_ratio()*mvVec;

        glutPostRedisplay();
    } else if (button == GLUT_LEFT_BUTTON) {
        mouse_state.func = MouseState::ROTATE;
    } else if (button == GLUT_MIDDLE_BUTTON) {
        mouse_state.func = MouseState::TRANSLATE;
    }
}

void motion (int x, int y) {
    if (!mouse_state.down)
        return;
    Pane* pane = Pane::current();
    if (!pane) return;

    if (mouse_state.func == MouseState::ROTATE) {
        double speed = 0.25;
        pane->lon += (x - mouse_state.x)*speed;
        pane->lat += (y - mouse_state.y)*speed;
        //pane->lat = clamp(pane->lat, -90., 90.);
    } else if (mouse_state.func == MouseState::TRANSLATE) {
        double speed = 1e-3;
        pane->offset[0] += (x - mouse_state.x)*speed;
        pane->offset[1] -= (y - mouse_state.y)*speed;
    }
    mouse_state.x = x;
    mouse_state.y = y;
    glutPostRedisplay();
}

void nop () {} // apparently needed by GLUT 3.0

void (*keyboard_sub)(unsigned char, int, int);
void keyboard_handler(unsigned char key, int x, int y) {
	if (key == 'z') {
        zoom(true);
    } else if (key == 'n') {
        stepDebug = false;
    } else if (key == 'x') {
        zoom(false);
    } else if (key == '[') {
    	display_modes[display_mode].scale *= 0.5;
    } else if (key == ']') {
    	display_modes[display_mode].scale *= 2.0;
    } else if (key == '\t') {
    	for (display_mode++; !display_modes[display_mode].active; display_mode++) {
    		if (display_modes[display_mode].scale == 0) // sentinel
    			display_mode = -1;
    	}
    } else
    	keyboard_sub(key,x,y);
    redisplay();
}

void init_glut (const GlutCallbacks &cb) {
    int argc = 1;
    char argv0[] = "";
    char *argv = argv0;
    glutInit(&argc, &argv);
    glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH|GLUT_MULTISAMPLE);
    glutInitWindowSize(1280,720);
    int mainWindow = glutCreateWindow("ARCSim");
    glutDisplayFunc(nop);
    glutReshapeFunc(reshape);
    glutIdleFunc(cb.idle);
    keyboard_sub = cb.keyboard; 
    glutKeyboardFunc(keyboard_handler);
    glutSpecialFunc(cb.special);
    double x[4] = {0, 1280/3, 1280*2/3, 1280};
    void (*display[3]) () = {display_material, display_plastic, display_world};
    for (int i = 0; i < 3; i++) {
        Pane::panes[i].window = glutCreateSubWindow(mainWindow, x[i],0, x[i+1],720);
        Pane::panes[i].parent = mainWindow;
        glutDisplayFunc(display[i]);
        glutKeyboardFunc(keyboard_handler);
        glutSpecialFunc(cb.special);
        glutMouseFunc(mouse);
        glutMotionFunc(motion);
    }
//    Pane::plastic().enabled = sim.enabled[Simulation::Plasticity];
}

void run_glut() {
    glutMainLoop();
}

void redisplay () {
    for (int i = 0; i < 3; i++) {
        glutSetWindow(Pane::panes[i].window);
        glutPostRedisplay();
    }
}

void wait_key () {   
    redisplay();
    stepDebug = true;

    while (stepDebug) { 
        glutMainLoopEvent();      
    };
    Annotation::list.clear();
}

#endif // NO_OPENGL

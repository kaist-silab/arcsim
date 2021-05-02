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

#include "display.hpp"
#include "opengl.hpp"
#include "util.hpp"
#include <boost/filesystem.hpp>
#include <cassert>
#include <cfloat>
#include <json/json.h>
#include <fstream>
#include <png.h>
#include <sstream>
using namespace std;


const int FILE_VERSION = 1;


// OBJ meshes

void get_valid_line (istream &in, string &line) {
    do
        getline(in, line);
    while (in && (line.length() == 0 || line[0] == '#'));
}

void triangle_to_obj (const string &inname, const string &outname) {
    fstream outfile(outname.c_str(), ios::out);
    { // nodes
        string filename = inname + ".node";
        fstream file(filename.c_str(), ios::in);
        string line;
        get_valid_line(file, line);
        stringstream linestream(line);
        int nv, dim, na, nb;
        linestream >> nv >> dim >> na >> nb;
        for (int i = 0; i < nv; i++) {
            get_valid_line(file, line);
            stringstream linestream(line);
            int index;
            linestream >> index;
            Vec2 u;
            linestream >> u[0] >> u[1];
            outfile << "v " << u[0] << " " << u[1] << " " << 0 << endl;
        }
    }
    { // eles
        string filename = inname + ".ele";
        fstream file(filename.c_str(), ios::in);
        string line;
        get_valid_line(file, line);
        stringstream linestream(line);
        int nt, nn, na;
        linestream >> nt >> nn >> na;
        for (int i = 0; i < nt; i++) {
            get_valid_line(file, line);
            stringstream linestream(line);
            int index;
            linestream >> index;
            int v0, v1, v2;
            linestream >> v0 >> v1 >> v2;
            outfile << "f " << v0+1 << " " << v1+1 << " " << v2+1 << endl;
        }
    }
}

vector<Face*> triangulate (const vector<Vert*> &verts);

void load_obj (Mesh &mesh, const string &filename) {
    delete_mesh(mesh);
    fstream file(filename.c_str(), ios::in);
    if(!file) {
        cout << "Error: failed to open file " << filename << endl;
        return;
    }
    while (file) {
        string line;
        get_valid_line(file, line);
        stringstream linestream(line);
        string keyword;
        linestream >> keyword;
        if (keyword == "vt") {
            Vec2 u;
            linestream >> u[0] >> u[1];
            mesh.add(new Vert(expand_xy(u)));
        } else if (keyword == "ms") {
            Vec3 u;
            linestream >> u[0] >> u[1] >> u[2];
            mesh.add(new Vert(u));
        } else if (keyword == "v") {
            Vec3 x;
            linestream >> x[0] >> x[1] >> x[2];
            mesh.add(new Node(x, x, Vec3(0), 0, 0, false));
        } else if (keyword == "ny") {
            Vec3 &y = mesh.nodes.back()->y;
            linestream >> y[0] >> y[1] >> y[2];
        } else if (keyword == "nv") {
            Vec3 &v = mesh.nodes.back()->v;
            linestream >> v[0] >> v[1] >> v[2];
        } else if (keyword == "nl") {
            linestream >> mesh.nodes.back()->label;
        } else if (keyword == "e") {
            int n0, n1;
            linestream >> n0 >> n1;
            mesh.add(new Edge(mesh.nodes[n0-1], mesh.nodes[n1-1], 0, 0));
        } else if (keyword == "ea") {
            linestream >> mesh.edges.back()->theta_ideal;
        } else if (keyword == "ed") {
            linestream >> mesh.edges.back()->damage;
        } else if (keyword == "ep") {
            linestream >> mesh.edges.back()->preserve;
        } else if (keyword == "f") {
            vector<Vert*> verts;
            vector<Node*> nodes;
            string w;
            while (linestream >> w) {
                stringstream wstream(w);
                int v, n;
                char c;
                wstream >> n >> c >> v;
                nodes.push_back(mesh.nodes[n-1]);
                if (wstream) {
                    verts.push_back(mesh.verts[v-1]);
                } else if (!nodes.back()->verts.empty()) {
                    verts.push_back(nodes.back()->verts[0]);
                } else {
                	verts.push_back(new Vert(nodes.back()->x));
                    mesh.add(verts.back());
                }
            }
            for (int v = 0; v < (int)verts.size(); v++)
                connect(verts[v], nodes[v]);
            vector<Face*> faces = triangulate(verts);
            for (int f = 0; f < (int)faces.size(); f++)
                mesh.add(faces[f]);
        } else if (keyword == "tm") {
            linestream >> mesh.faces.back()->flag;
        } else if (keyword == "tp") {
            Mat3x3 &S = mesh.faces.back()->Sp_bend;
            for (int i=0; i<3; i++) for (int j=0; j<3; j++)
            	linestream >> S(i,j);
        } else if (keyword == "ts") {
            Mat3x3 &S = mesh.faces.back()->Sp_str;
            for (int i=0; i<3; i++) for (int j=0; j<3; j++)
            	linestream >> S(i,j);
        } else if (keyword == "td") {
            linestream >> mesh.faces.back()->damage;
        }
    }
    mark_nodes_to_preserve(mesh);
    compute_ms_data(mesh);
}

void load_objs (vector<Mesh*> &meshes, const string &prefix) {
    for (int m = 0; m < (int)meshes.size(); m++)
        load_obj(*meshes[m], stringf("%s_%02d.obj", prefix.c_str(), m));
}

static double angle (const Vec3 &x0, const Vec3 &x1, const Vec3 &x2) {
    Vec3 e1 = normalize(x1 - x0);
    Vec3 e2 = normalize(x2 - x0);
    return acos(clamp(dot(e1, e2), -1., 1.));
}

vector<Face*> triangulate (const vector<Vert*> &verts) {
    int n = verts.size();
    double best_min_angle = 0;
    int best_root = -1;
    for (int i = 0; i < n; i++) {
        double min_angle = infinity;
        const Vert *vert0 = verts[i];
        for (int j = 2; j < n; j++) {
            const Vert *vert1 = verts[(i+j-1)%n], *vert2 = verts[(i+j)%n];
            min_angle=min(min_angle,
                          angle(vert0->node->x,vert1->node->x,vert2->node->x),
                          angle(vert1->node->x,vert2->node->x,vert0->node->x),
                          angle(vert2->node->x,vert0->node->x,vert1->node->x));
        }
        if (min_angle > best_min_angle) {
            best_min_angle = min_angle;
            best_root = i;
        }
    }
    int i = best_root;
    Vert* vert0 = verts[i];
    vector<Face*> tris;
    for (int j = 2; j < n; j++) {
        Vert *vert1 = verts[(i+j-1)%n], *vert2 = verts[(i+j)%n];
        tris.push_back(new Face(vert0, vert1, vert2, Mat3x3(1), Mat3x3(0), 0, 0));
    }
    return tris;
}

void save_obj (const Mesh &mesh, const string &filename) {
	set_indices((Mesh&)mesh);
    fstream file(filename.c_str(), ios::out);
    for (int v = 0; v < (int)mesh.verts.size(); v++) {
        const Vert *vert = mesh.verts[v];
        file << "ms " << vert->u[0] << " " << vert->u[1] << " " << vert->u[2] << endl;
        //if (vert->label)
        //    file << "vl " << vert->label << endl;
    }
    for (int n = 0; n < (int)mesh.nodes.size(); n++) {
        const Node *node = mesh.nodes[n];
        file << "v " << node->x[0] << " " << node->x[1] << " "
             << node->x[2] << endl;
        if (norm2(node->x - node->y))
            file << "ny " << node->y[0] << " " << node->y[1] << " "
                 << node->y[2] << endl;
        if (norm2(node->v))
            file << "nv " << node->v[0] << " " << node->v[1] << " "
                 << node->v[2] << endl;
        if (node->label)
            file << "nl " << node->label << endl;
    }
    for (int e = 0; e < (int)mesh.edges.size(); e++) {
        const Edge *edge = mesh.edges[e];
        if (edge->theta_ideal || edge->preserve) {
            file << "e " << edge->n[0]->index+1 << " " << edge->n[1]->index+1
                 << endl;
            if (edge->theta_ideal)
                file << "ea " << edge->theta_ideal << endl;
            if (edge->damage)
                file << "ed " << edge->damage << endl;
            if (edge->preserve)
                file << "ep " << edge->preserve << endl;
        }
    }
    for (int f = 0; f < (int)mesh.faces.size(); f++) {
        const Face *face = mesh.faces[f];
        file << "f " << face->v[0]->node->index+1 << "/" << face->v[0]->index+1
             << " " << face->v[1]->node->index+1 << "/" << face->v[1]->index+1
             << " " << face->v[2]->node->index+1 << "/" << face->v[2]->index+1
             << endl;
        if (face->material && mesh.parent) {
        	file << "tm " << find(face->material, mesh.parent->materials) << endl;
        }
        if (norm2_F(face->Sp_bend)) {
            const Mat3x3 &S = face->Sp_bend;
            file << "tp ";
            for (int i=0; i<3; i++) for (int j=0; j<3; j++)
            	file << S(i,j) << " ";
           	file << endl;
        }
        if (norm2_F(face->Sp_str)) {
            const Mat3x3 &S = face->Sp_str;
            file << "ts ";
            for (int i=0; i<3; i++) for (int j=0; j<3; j++)
            	file << S(i,j) << " ";
            file << endl;
        }
        if (face->damage)
            file << "td " << face->damage << endl;
    }
}

void save_objs (const vector<Mesh*> &meshes, const string &prefix) {
    for (int m = 0; m < (int)meshes.size(); m++)
        save_obj(*meshes[m], stringf("%s_%02d.obj", prefix.c_str(), m));
}

template<> void serializer<Simulation>(Simulation& sim, Serialize& s, const string& n) {
	for (size_t m = 0; m < sim.cloth_meshes.size(); m++) {
    	sim.cloth_meshes[m]->serializer(s);
    	if (s.load()) {
    		compute_ms_data(*sim.cloth_meshes[m]);
    	}
    }
	for (size_t m = 0; m < sim.obstacle_meshes.size(); m++)
    	sim.obstacle_meshes[m]->serializer(s);
}
void serialize_header(Serialize& s) {
	s.version = FILE_VERSION;
    if (s.load() || s.check()) {
		int sframe;    
		double stime;
		gzread(s.fp, &s.version, sizeof(int));
    	gzread(s.fp, &sframe, sizeof(int));
    	gzread(s.fp, &stime, sizeof(double));
    	if (s.load())
    		sim.time = stime;
    } else {
    	gzwrite(s.fp, &s.version, sizeof(int));
  		gzwrite(s.fp, &sim.frame, sizeof(int));
  		gzwrite(s.fp, &sim.time, sizeof(double));
	} 	
}

template<class T>
bool load_state (T& state, const string &prefix) {
    Serialize s;
    s.fp = gzopen( (prefix+".bin").c_str(), "r");
    if (!s.fp) {
    	return false;
    }    
    s.mode = Serialize::Load;
    serialize_header(s);
    serializer(state, s, "state");
    gzclose(s.fp);
    return true;
}

template<class T>
void save_state (T& state, const string &prefix) {
    Serialize s;
    s.fp = gzopen( (prefix+".bin").c_str(), "w3");
    if (!s.fp) {
    	cout << "can't open " << prefix << ".bin" << endl;
    	return;
    }
    s.mode = Serialize::Save;
	serialize_header(s);
    serializer(state, s, "state");
    gzclose(s.fp);
}
template bool load_state<Simulation>(Simulation&, const string&);
template void save_state<Simulation>(Simulation&, const string&);

string obtain_subframe_id() {
	// build consistent subframe indices
	static int id = 0, lastframe = -1;
	if (lastframe != sim.frame) {
		lastframe = sim.frame;
		id = 0;
	}
	return stringf("%d_%d",sim.frame, id++);
}

void save_transformation (const Transformation &tr, const string &filename) {
    FILE* file = fopen(filename.c_str(), "w");
    pair<Vec3, double> axis_angle = tr.rotation.to_axisangle();
    Vec3 axis = axis_angle.first;
    double angle = axis_angle.second * 180 / M_PI;
    fprintf(file, "<rotate angle=\"%f\" x=\"%f\" y=\"%f\" z=\"%f\"/>\n",
            angle, axis[0], axis[1], axis[2]);
    fprintf(file, "<scale value=\"%f\"/>\n", tr.scale);
    fprintf(file, "<translate x=\"%f\" y=\"%f\" z=\"%f\"/>\n",
            tr.translation[0], tr.translation[1], tr.translation[2]);
    fclose(file);
}

// Images

void flip_image (int w, int h, unsigned char *pixels);

void save_png (const char *filename, int width, int height,
               unsigned char *pixels, bool has_alpha = false);

#ifndef NO_OPENGL
void save_screenshot (const string &filename) {
    int w = 0, h = 0;
    for (int s = 0; s < 3; s++) {
        glutSetWindow(Pane::panes[s].window);
        w += glutGet(GLUT_WINDOW_WIDTH);
        h = max(h, glutGet(GLUT_WINDOW_HEIGHT));
    }
    unsigned char *pixels = new unsigned char[w*h*3];
    int x = 0;
    for (int s = 0; s < 3; s++) {
        glutSetWindow(Pane::panes[s].window);
        int wsub = glutGet(GLUT_WINDOW_WIDTH),
            hsub = glutGet(GLUT_WINDOW_HEIGHT);
        unsigned char *pixelsub = new unsigned char[wsub*hsub*3];
        glPixelStorei(GL_PACK_ALIGNMENT, 1);
        glReadPixels(0,0, wsub,hsub, GL_RGB, GL_UNSIGNED_BYTE, pixelsub);
        for (int j = 0; j < hsub; j++)
            for (int i = 0; i < wsub; i++)
                for (int c = 0; c < 3; c++)
                    pixels[(x+i+w*j)*3+c] = pixelsub[(i+wsub*j)*3+c];
        x += wsub;
        delete[] pixelsub;
    }
    flip_image(w,h, pixels);
    save_png(filename.c_str(), w,h, pixels);
    delete[] pixels;
}
#endif

void flip_image (int w, int h, unsigned char *pixels) {
    for (int j = 0; j < h/2; j++)
        for (int i = 0; i < w; i++)
            for (int c = 0; c < 3; c++)
                swap(pixels[(i+w*j)*3+c], pixels[(i+w*(h-1-j))*3+c]);
}

void save_png (const char *filename, int width, int height,
               unsigned char *pixels, bool has_alpha) {
#ifndef _WIN32
    FILE* file = fopen(filename, "wb");
    if (!file) {
        printf("Couldn't open file %s for writing.\n", filename);
        return;
    }
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL,
                                                  NULL, NULL);
    if (!png_ptr) {
        printf("Couldn't create a PNG write structure.\n");
        fclose(file);
        return;
    }
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        printf("Couldn't create a PNG info structure.\n");
        png_destroy_write_struct(&png_ptr, NULL);
        fclose(file);
        return;
    }
    if (setjmp(png_jmpbuf(png_ptr))) {
        printf("Had a problem writing %s.\n", filename);
        png_destroy_write_struct(&png_ptr, &info_ptr);
        fclose(file);
        return;
    }
    png_init_io(png_ptr, file);
    png_set_IHDR(png_ptr, info_ptr, width, height, 8,
                 has_alpha ? PNG_COLOR_TYPE_RGBA : PNG_COLOR_TYPE_RGB,
                 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT);
    int channels = has_alpha ? 4 : 3;
    png_bytep* row_pointers = (png_bytep*) new unsigned char*[height];
    for (int y = 0; y < height; y++)
        row_pointers[y] = (png_bytep) &pixels[y*width*channels];
    png_set_rows(png_ptr, info_ptr, row_pointers);
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
    delete[] row_pointers;
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(file);
#endif
}

void ensure_existing_directory (const std::string &path) {
    return;
    using namespace boost::filesystem;
    if (!exists(path))
        create_directory(path);
    if (!is_directory(path)) {
        cout << "Error: " << path << " is not a directory!" << endl;
        abort();
    }
}

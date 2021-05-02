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

#include "util.hpp"
#include "io.hpp"
#include "mesh.hpp"
#include <algorithm>
#include <cstdarg>
#include <cstdio>
#include <iomanip>
#include <limits>
#include <map>
#include <signal.h>
#include <sstream>
using namespace std;

void Stats::add (double x) {
    xs.push_back(x);
    sum += x;
    sorted = false;
}

void Stats::sort () const {
    if (sorted) return;
    std::sort(xs.begin(), xs.end());
    sorted = true;
}

double Stats::min () const {sort(); return xs.front();}
double Stats::max () const {sort(); return xs.back();}
double Stats::mean () const {return sum/xs.size();}
double Stats::median () const {return quantile(0.5);}
double Stats::quantile (double q) const {sort(); return xs[(int)(q*xs.size())];}

ostream &operator<< (ostream &out, const Stats &stats) {
    if (stats.xs.empty())
        out << "no data";
    else
        out << stats.min() << " " << stats.quantile(0.05) << " "
            << stats.quantile(0.25) << " " << stats.median() << " "
            << stats.quantile(0.75) << " " << stats.quantile(0.95) << " "
            << stats.max();
    return out;
}

inline string stringf (const string &format, ...) {
    char buf[256];
    va_list args;
    va_start(args, format);
    vsnprintf(buf, 256, format.c_str(), args);
    va_end(args);
    return std::string(buf);
}

template <typename T> string name (const T *p) {
    stringstream ss;
    ss << setw(3) << setfill('0') << hex << ((size_t)p/sizeof(T))%0xfff;
    return ss.str();
}

ostream &operator<< (ostream &out, const Vert *vert) {
	if (!vert)
		out << "v:none";
    else
    	out << "v:" << name(vert); 
    return out;
}

ostream &operator<< (ostream &out, const Node *node) {
	if (!node)
		out << "n:none";
	else
    	out << "n:" << name(node) << node->verts; 
    return out;
}

ostream &operator<< (ostream &out, const Edge *edge) {
    if (!edge)
    	out << "e:none";
    else 
    	out << "e:" << name(edge) << "(" << edge->n[0] << "-" << edge->n[1] << ")"; 
    return out;
}

ostream &operator<< (ostream &out, const Face *face) {
	if (!face)
		out << "f:none";
	else
    	out << "f:" << name(face) << "(" << face->v[0] << "-" << face->v[1] << "-" << face->v[2] << ")"; return out;
}

const double infinity = numeric_limits<double>::infinity();

int solve_quadratic (double a, double b, double c, double x[2]) {
    // http://en.wikipedia.org/wiki/Quadratic_formula#Floating_point_implementation
    double d = b*b - 4*a*c;
    if (d < 0) {
        x[0] = -b/(2*a);
        return 0;
    }
    double q = -(b + sgn(b)*sqrt(d))/2;
    int i = 0;
    if (abs(a) > 1e-12*abs(q))
        x[i++] = q/a;
    if (abs(q) > 1e-12*abs(c))
        x[i++] = c/q;
    if (i==2 && x[0] > x[1])
        swap(x[0], x[1]);
    return i;
}

bool is_seam_or_boundary (const Vert *v) {
    return is_seam_or_boundary(v->node);
}

bool is_seam_or_boundary (const Node *n) {
    for (int e = 0; e < (int)n->adje.size(); e++)
        if (is_seam_or_boundary(n->adje[e]))
            return true;
    return false;
}

bool is_seam_or_boundary (const Edge *e) {
    return !e->adjf[0] || !e->adjf[1]
        || edge_vert(e,0,0) != edge_vert(e,1,0)
        || edge_vert(e,0,1) != edge_vert(e,1,1);
}

bool is_seam_or_boundary (const Face *f) {
    return is_seam_or_boundary(f->adje[0])
        || is_seam_or_boundary(f->adje[1])
        || is_seam_or_boundary(f->adje[2]);
}

bool is_seam (const Edge* e) {
	return e->adjf[0] && e->adjf[1]
        && (edge_vert(e,0,0) != edge_vert(e,1,0)
            || edge_vert(e,0,1) != edge_vert(e,1,1));
}

void build_node_lookup(map<const Node*,Vec3>& nodemap, const vector<Mesh*>& meshes) {
	for (size_t i=0; i<meshes.size(); i++)
    	for (size_t j=0; j<meshes[i]->nodes.size(); j++)
    		nodemap[meshes[i]->nodes[j]] = meshes[i]->nodes[j]->x;    
}


void segfault() {
	raise(SIGSEGV);
}

void debug_save_meshes (const vector<Mesh*> &meshvec, const string &name,
                        int n) {
    static map<string,int> savecount;
    if (n == -1)
        n = savecount[name];
    else
        savecount[name] = n;
    save_objs(meshvec, stringf("tmp/%s%d", name.c_str(), n));
    savecount[name] = (savecount[name] + 1) % 10;
}

void debug_save_mesh (const Mesh &mesh, const string &name, int n) {
    static map<string,int> savecount;
    if (n == -1)
        n = savecount[name];
    else
        savecount[name] = n;
    save_obj(mesh, stringf("tmp/%s%d.obj", name.c_str(), n));
    savecount[name] = (savecount[name] + 1) % 10;
}

template<class T>
void serializer(T& v, Serialize& s, const string& name) { 
	if (s.load())
		gzread(s.fp, &v, sizeof(T));
	else if (s.save())
		gzwrite(s.fp, &v, sizeof(T));
	else if (s.check()) {
		T tmp;
		gzread(s.fp, &tmp, sizeof(T));
		if (tmp != v) {
			cout << "Serialization check failed: " << name << ": file " << 
					tmp << " mem " << v << std::endl;
			cout << "difference is " << (v-tmp) << endl;
			segfault();
		}
	}
}

template<class T> 
void serializer_dvec(vector<T>& v, Serialize& s, const string& name) { 
	serializer_array(v, s, name);
	for(size_t i=0; i<v.size(); i++) 
		serializer(v[i], s, name); 
}

template<int n> 
void serializer_vec(Vec<n>& x, Serialize& s, const string& name) { 
	for(int i=0; i<n; i++) 
		serializer(x[i], s, name); 
}

template<int n> 
void serializer_mat(Mat<n,n>& x, Serialize& s, const string& name) { 
	for(int i=0; i<n; i++) 
		serializer(x.col(i), s, name); 
}

template void serializer<int>(int&,Serialize&,const string&);
template void serializer<double>(double&,Serialize&,const string&);
template void serializer<bool>(bool&,Serialize&,const string&);
template<> void serializer<Vec2>(Vec2& x, Serialize& s,const string& n) { return serializer_vec(x,s,n); }
template<> void serializer<Vec3>(Vec3& x, Serialize &s,const string& n) { return serializer_vec(x,s,n); }
template<> void serializer<Mat2x2>(Mat2x2& x, Serialize& s,const string& n) { return serializer_mat(x,s,n); }
template<> void serializer<Mat3x3>(Mat3x3& x, Serialize& s,const string& n) { return serializer_mat(x,s,n); }
template<> void serializer<vector<Vec3> >(vector<Vec3>& x, Serialize& s, const string& n) { return serializer_dvec(x,s,n); }
template<> void serializer<vector<double> >(vector<double>& x, Serialize& s, const string& n) { return serializer_dvec(x,s,n); }

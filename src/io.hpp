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

#ifndef IO_HPP
#define IO_HPP

#include "mesh.hpp"
#include "util.hpp"
#include <zlib.h>

void triangle_to_obj (const std::string &infile, const std::string &outfile);

void load_obj (Mesh &mesh, const std::string &filename);
void load_objs (std::vector<Mesh*> &meshes, const std::string &prefix);

void save_obj (const Mesh &mesh, const std::string &filename);
void save_objs (const std::vector<Mesh*> &meshes, const std::string &prefix);

template<class T>
void test_state (T& state, const std::string &prefix);
template<class T>
void save_state (T& state, const std::string &prefix);
template<class T>
bool load_state (T& state, const std::string &prefix);

void save_transformation (const Transformation &tr,
                          const std::string &filename);

// w_crop and h_crop specify a multiplicative crop window
void save_screenshot (const std::string &filename);

// check that output directory exists; if not, create it
void ensure_existing_directory (const std::string &path);




// IMPLEMENTATION

std::string obtain_subframe_id();
void serialize_header(Serialize& s);

template<class T>
void test_state (T& state, const std::string& prefix) {
	// build consistent subframe indice
	std::string id = obtain_subframe_id();
	std::string cname = prefix + "_" + id + ".chk";
	
	Serialize s;
    s.mode = Serialize::Check;
    s.fp = gzopen(cname.c_str(), "r");
    if (!s.fp) {
    	s.mode = Serialize::Save;
    	s.fp = gzopen(cname.c_str(), "w3");
    	if (!s.fp) {
    		std::cout << "can't write check file " << cname << std::endl;
    		exit(1);
    	}
    }	
    std::cout << "CHECK: frame " << id << " [" << (s.save() ? "WRITE" : "READ") << "]\n";
    serialize_header(s);
    serializer(state, s);    
    gzclose(s.fp);
}



#endif

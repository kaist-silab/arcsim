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

#include "displayreplay.hpp"

#include "conf.hpp"
#include "display.hpp"
#include "io.hpp"
#include "misc.hpp"
#include "opengl.hpp"
#include <cstdio>
#include <fstream>
using namespace std;

#ifndef NO_OPENGL

extern string inprefix;
extern string outprefix;
static int frameskip;

static bool running = false;

static void reload () {
	Annotation::list.clear();
    int fullframe = ::frame*::frameskip;
    sim.frame = ::frame;
    sim.time = fullframe * sim.frame_time;
    if (!load_state(sim, stringf("%s/%05d",inprefix.c_str(), fullframe)) ||
        sim.cloth_meshes[0]->verts.empty()) {
        if (::frame == 0)
            exit(EXIT_FAILURE);
        if (!outprefix.empty())
            exit(EXIT_SUCCESS);
        ::frame = 0;
        reload();
    }
    for (int o = 0; o < (int)sim.obstacles.size(); o++)
        sim.obstacles[o].get_mesh(sim.time);
}

static void idle () {
    if (!running)
        return;
    fps.tick();
    if (!outprefix.empty()) {
        char filename[256];
        snprintf(filename, 256, "%s/%05d.png", outprefix.c_str(), ::frame );
        save_screenshot(filename);
    }
    ::frame += 1;
    reload();
    fps.tock();
    redisplay();
}

static void keyboard (unsigned char key, int x, int y) {
    unsigned char esc = 27, space = ' ';
    if (key == esc) {
        exit(0);
    } else if (key == space) {
        running = !running;
    }
}

static void special (int key, int x, int y) {
    bool shift = glutGetModifiers() & GLUT_ACTIVE_SHIFT,
         alt = glutGetModifiers() & GLUT_ACTIVE_ALT;
    int delta = alt ? 100 : shift ? 10 : 1;
    if (key == GLUT_KEY_LEFT) {
        ::frame -= delta;
        reload();
    } else if (key == GLUT_KEY_RIGHT) {
        ::frame += delta;
        reload();
    } else if (key == GLUT_KEY_HOME) {
        ::frame = 0;
        reload();
    }
    redisplay();
}

static void save_obstacle_transforms (const vector<Obstacle> &obs, int frame,
                                      double time) {
    if (!outprefix.empty() && frame < 100000) {
        for (int o = 0; o < (int)obs.size(); o++) {
            Transformation trans = identity();
            if (obs[o].transform_spline)
                trans = get_dtrans(*obs[o].transform_spline, time).first;
            save_transformation(trans, stringf("%s/%05dobs%02d.txt",
                                               outprefix.c_str(), frame, o));
        }
    }
}

void generate_obj (const vector<string> &args) {
    if (args.size() < 1 || args.size() > 2) {
        cout << "Generates output meshes from saved simulation state." << endl;
        cout << "Arguments:" << endl;
        cout << "    <out-dir>: Directory containing simulation output files"
             << endl;
        exit(EXIT_FAILURE);
    }
    ::inprefix = args[0];
    ::outprefix = args[0];
    ::frameskip = 1;
    if (!::outprefix.empty())
        ensure_existing_directory(::outprefix);
    char config_backup_name[256];
    snprintf(config_backup_name, 256, "%s/%s", inprefix.c_str(), "conf.json");
    load_json(config_backup_name, sim);
    prepare(sim);
    
    for (;;) {
        reload();
        char filename[256];
        snprintf(filename, 256, "%s/%05d", inprefix.c_str(), ::frame);
        save_objs(sim.cloth_meshes, filename);
        save_obstacle_transforms(sim.obstacles, frame, sim.time);
        ::frame += 1;
    }
}

void display_replay (const vector<string> &args) {
    if (args.size() < 1 || args.size() > 2) {
        cout << "Replays the results of a simulation." << endl;
        cout << "Arguments:" << endl;
        cout << "    <out-dir>: Directory containing simulation output files"
             << endl;
        cout << "    <sshot-dir> (optional): Directory to save images" << endl;
        exit(EXIT_FAILURE);
    }
    ::inprefix = args[0];
    ::outprefix = args.size()>1 ? args[1] : "";
    ::frameskip = 1;
    if (!::outprefix.empty())
        ensure_existing_directory(::outprefix);
    char config_backup_name[256];
    snprintf(config_backup_name, 256, "%s/%s", inprefix.c_str(), "conf.json");
    load_json(config_backup_name, sim);
    prepare(sim);
    reload();
    GlutCallbacks cb;
    cb.idle = idle;
    cb.keyboard = keyboard;
    cb.special = special;
    init_glut(cb);
    run_glut();
}

#else

void display_replay (const vector<string> &args) {opengl_fail();}

#endif // NO_OPENGL

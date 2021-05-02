#include "sepstrength.hpp"
#include "display.hpp"
#include "geometry.hpp"
#include "physics.hpp"

using namespace std;

struct FanPrecomp {
	Vec3 t0, t1;
	Vec3 frame_x, frame_y;
	double angle0, angle1;
	Face* face;
};

Mat3x3 compute_sigma(const Face* face) {/*
	Mat3x3 F_str = deformation_gradient<WS>(face);
    Mat3x3 G_str = (F_str.t()*F_str - Mat3x3(1)) * 0.5;
    Mat3x3 sigma_str = material_model(face, G_str);

    // add explicit bending
    Mat3x3 F_bend = Mat3x3(1);
    for (int i=0; i<3; i++) 
    	F_bend += face->v[i]->node->curvature * (0.5 * face->material->fracture_bend_thickness);
    Mat3x3 G_bend = (F_bend.t()*F_bend - Mat3x3(1)) * 0.5;
    Mat3x3 sigma_bend = material_model(face, G_bend);

    return get_positive(sigma_str) + get_positive(sigma_bend);*/
    Mat3x3 F_str = deformation_gradient<WS>(face);
    //Mat3x3 G_str = (F_str.t()*F_str - Mat3x3(1)) * 0.5;
    //Mat3x3 sigma_str = material_model(face, G_str);

    // add explicit bending
    Mat3x3 F_bend = Mat3x3(1);
    for (int i=0; i<3; i++) 
        F_bend += face->v[i]->node->curvature * (0.5 * face->material->fracture_bend_thickness);
    F_bend = F_str * F_bend;
    Mat3x3 G_bend = (F_bend.t()*F_bend - Mat3x3(1)) * 0.5;
    Mat3x3 sigma_bend = material_model(face, G_bend);

    return get_positive(sigma_bend);
}

static Vec3 get_ms_span (Face* f, Node* center, int offset) {
	int s = f->v[0]->node == center ? 0 : (f->v[1]->node == center ? 1 : 2);
	int s1 = NEXT(s);
	if (offset>0) s1 = NEXT(s1);
	return f->v[s1]->u - f->v[s]->u;
}

double separation_strength(Node* node, SplitNode* split, bool always_compute) {

	// first start/end edge
	Edge *start = NULL, *end = NULL;
	for (size_t e=0; e<node->adje.size(); e++) {
		Edge* edge = node->adje[e];
	    if (next_face_ccw(edge, node) == NULL)
	        end = edge;
		else if (next_face_cw(edge, node) == NULL)
	        start = edge;
	}
	bool has_crack = end || start;
	if (!start) start = node->adje[0];
	if (!end) end = node->adje[0];

	static vector<FanPrecomp> fan;
	if (fan.size() < node->adje.size())
		fan.resize(node->adje.size());

	// precompute tensor products and check upper bound
	double max_sigma = 0, min_toughness = infinity;
	double end_angle = 0;
	int num = 0;
	Vec3 n = normal<MS>(node);
    for(Edge* cur = start;;cur=next_edge_ccw(cur,node)) {
    	Face* face = next_face_ccw(cur, node);
    	if (!face || (num > 0 && cur==end)) break;
    	FanPrecomp& f = fan[num++];
    	f.face = face;
    	//face->sigma = Mat3x3(1);
    	Vec3 eig = eigen_values(face->sigma);
    	max_sigma = max(max_sigma, eig[0]);
    	min_toughness = min(min_toughness, face->material->toughness);
    	Vec3 u = get_ms_span(face, node, 0);
        Vec3 v = get_ms_span(face, node, 1);
        f.angle0 = end_angle;
        end_angle += get_angle(u,v);        
        f.angle1 = end_angle;
        f.t0 = face->sigma * normalize(cross(n,u));
        f.t1 = face->sigma * normalize(cross(n,v));
        f.frame_x = normalize(u);
        f.frame_y = normalize(v - dot(v,f.frame_x) * f.frame_x);
    }
    // upper bound
    if (has_crack) 
        min_toughness *= 0.5;
    if (2*max_sigma < min_toughness && !always_compute) 
        return 0;

    SplitNode maxSplit(node);

    // discretize (CCW)
    const int num_disc = 200;
    for (int i=0; i<num_disc; i++) {
        double angle_fwd = (double)i*end_angle/num_disc;
        double angle_bwd = angle_fwd+M_PI;
        if (angle_bwd > 2*M_PI) angle_bwd -= 2*M_PI;

        // integrate traction
        int sector = 0;
        Vec3 traction[2] = {Vec3(0), Vec3(0)};
        Vec3 normal[2] = {Vec3(0), Vec3(0)};
        Face* faces[2] = {NULL, NULL};
        		
        for (int f=0; f < num; f++) {
	        if (angle_fwd >= fan[f].angle0 && angle_fwd < fan[f].angle1) {
        		double phi = angle_fwd - fan[f].angle0;
        		Vec3 tn = cos(phi) * fan[f].frame_x + sin(phi) * fan[f].frame_y;
        		Vec3 tmult = fan[f].face->sigma * normalize(cross(n,tn));
        		traction[sector] += tmult - fan[f].t0;
        		sector = 1 - sector;
        		traction[sector] += fan[f].t1 - tmult;
        		normal[faces[0] ? 1:0] = tn;
        		faces[faces[0] ? 1:0] = fan[f].face;

        	} else if (!has_crack && angle_bwd >= fan[f].angle0 && angle_bwd < fan[f].angle1) {
        		double phi = angle_bwd - fan[f].angle0;
        		Vec3 tn = cos(phi) * fan[f].frame_x + sin(phi) * fan[f].frame_y;
        		Vec3 tmult = fan[f].face->sigma * normalize(cross(n,tn));
        		traction[sector] += tmult - fan[f].t0;
        		sector = 1 - sector;
        		traction[sector] += fan[f].t1 - tmult;	
        		normal[faces[0] ? 1:0] = tn;
        		faces[faces[0] ? 1:0] = fan[f].face;
        	} else {
        		traction[sector] += fan[f].t1 - fan[f].t0;
        	}
        	//if (node->index == I) wait_key();
        }
        if (!faces[0]) continue;
        if (!has_crack && (!faces[1])) continue;

        // project on common mid-vector            
        double sep = 0;
        if (dot(traction[0], traction[1]) < 0) {
            Vec3 mid = normalize(traction[0] - traction[1]);
            sep = 0.5 * min(fabs(dot(traction[0], mid)), fabs(dot(traction[1], mid)));
            if (has_crack)
            	sep *= 2.0;
        }

        // divide by toughness
        double toughness = faces[0]->material->toughness;
        if (faces[1]) toughness = 0.5*(faces[1]->material->toughness + toughness);
        sep /= toughness;

        if (sep > maxSplit.sep) {
        	maxSplit.sep = sep;
        	maxSplit.normal[0] = normal[0];
        	maxSplit.normal[1] = normal[1];
        	maxSplit.faces[0] = faces[0];
        	maxSplit.faces[1] = faces[1];        	
        }
   	}
    node->sep = maxSplit.sep;
   	if (split)
   		*split = maxSplit;
    
    return maxSplit.sep;
}

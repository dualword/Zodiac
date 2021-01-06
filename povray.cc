/*   Zodiac - molecular modelling environment. www.zeden.org
    Copyright (C) 2008  Nicola Zonta (nicola.zonta(at)zeden.org)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/
#include "ddwin.h"

#define WIREFRAME 0
#define STICK 1
#define CPK 2
#define BALLANDSTICK 3
#define BALLANDLINE 4

#define NO_ATOMS 0
#define SPHERES 1
#define CPK_SPHERES 2
#define SCALED_CPK_SPHERES 3

#define NO_BONDS 0
#define LINES 1
#define STICKS 2

#define BACKBONE_LINE 0
#define BACKBONE_STICK 1

#define AROMATIC_RINGS 0
#define KEKULE 1
#define AROMATIC_BONDS 2



string DDWin::POV_bond (Bond *bo) {
    stringstream out;
    bool vis = get_visible (bo);
    int ds = get_ds (bo);
    color col1 = get_color (bo -> GetBeginAtom ());
    color col2 = get_color (bo -> GetEndAtom   ());
    double vdw1 = get_vdw (bo -> GetBeginAtom ());
    double vdw2 = get_vdw (bo -> GetEndAtom   ());
    if (vis) {


    if (ds != NO_BONDS) {
        if (ds == STICKS || ds == LINES) {
            double a0size = vdw1 * gl->vdw_scale;
            double a1size = vdw2 * gl->vdw_scale;
            float stick_rad = gl->stick_rad;
            if (ds == LINES)  stick_rad = 0.01f;   
            int order;
            if (gl->aromatic_display_style==KEKULE) order=bo->GetBondOrder ();
            else order = bo -> GetBondOrder ();              
            if (order == 1 || order == 4 || order == 5) {
                vect c1 = get_coordinates(bo->GetBeginAtom ());
                vect c2 = get_coordinates(bo->GetEndAtom ());
				out <<POV_stick(c1, c2, col1, col2, stick_rad, a0size, a1size);
            }
            else if (order == 2) {
                float verts  [4][3];
                gl->compute_double_bond_vertexes (bo, verts);




                float x1 = verts[0][0];
                float y1 = verts[0][1];
                float z1 = verts[0][2];
                float x2 = verts[2][0];
                float y2 = verts[2][1];
                float z2 = verts[2][2];
                out << "merge {"<<endl;
                out << "    sphere {";
                out << "<"<<x1<<","<<y1<<","<<z1<<">,"<<stick_rad/2<<endl;
                out << "pigment { "<<POV_color (col1)<<"} finish {ambient 0.1 diffuse 0.7 phong 0.3 phong_size 30 }}"<<endl;    
                out << "    sphere {";
                out << "<"<<x2<<","<<y2<<","<<z2<<">,"<<stick_rad/2<<endl;
                out << "pigment { "<<POV_color (col2)<<"} finish {ambient 0.1 diffuse 0.7 phong 0.3 phong_size 30 }}"<<endl;    
                out << "     cylinder {"<<endl; 
                out << "<"<<x1<<","<<y1<<","<<z1<<">,";
                out << "<"<<x2<<","<<y2<<","<<z2<<">,"<<stick_rad/2<<endl;

                float gradx, grady, gradz, dist;
                gradx = x2 - x1;
                grady = y2 - y1;
                gradz = z2 - z1;
                dist = sqrt (gradx*gradx + grady*grady + gradz*gradz);
                out << "pigment {gradient <"<< gradx/dist<<","<< grady/dist<<","<< gradz/dist<<">\n  color_map {["<<a0size/dist<<" "<<POV_color (col1)<<"] ["<<1-a1size/dist<<" "<<POV_color (col2)<<"]} scale "<<dist<<" translate < "<<x1<<","<<y1<<","<<z1<<">}finish {ambient 0.1 diffuse 0.7 phong 0.3 phong_size 30 }}"<<endl;
  
     		out << "no_shadow" << endl;
                out << "}"<<endl;



                x1 = verts[1][0];
                y1 = verts[1][1];
                z1 = verts[1][2];
                x2 = verts[3][0];
                y2 = verts[3][1];
                z2 = verts[3][2];
                out << "merge {"<<endl;
                out << "    sphere {";
                out << "<"<<x1<<","<<y1<<","<<z1<<">,"<<stick_rad/2<<endl;
                out << "pigment { "<<POV_color (col1)<<"} finish {ambient 0.1 diffuse 0.7 phong 0.3 phong_size 30 }}"<<endl;    

                out << "    sphere {";
                out << "<"<<x2<<","<<y2<<","<<z2<<">,"<<stick_rad/2<<endl;
                out << "pigment { "<<POV_color (col2)<<"} finish {ambient 0.1 diffuse 0.7 phong 0.3 phong_size 30 }}"<<endl;    
                out << "     cylinder {"<<endl; 
                out << "<"<<x1<<","<<y1<<","<<z1<<">,";
                out << "<"<<x2<<","<<y2<<","<<z2<<">,"<<stick_rad/2<<endl;

//                float gradx, grady, gradz, dist;
                gradx = x2 - x1;
                grady = y2 - y1;
                gradz = z2 - z1;
                dist = sqrt (gradx*gradx + grady*grady + gradz*gradz);
                out << "pigment {gradient <"<< gradx/dist<<","<< grady/dist<<","<< gradz/dist<<">\n  color_map {["<<a0size/dist<<" "<<POV_color (col1)<<"] ["<<1-a1size/dist<<" "<<POV_color (col2)<<"]} scale "<<dist<<" translate < "<<x1<<","<<y1<<","<<z1<<">}finish {ambient 0.1 diffuse 0.7 phong 0.3 phong_size 30 }}"<<endl;
     		out << "no_shadow" << endl;
                out << "}"<<endl;





/*
                out << "cylinder {"<<endl; 
                out << "<"<<verts[0][0]<<","<<verts[0][1]<<","<<verts[0][2]<<">,";
                out << "<"<<verts[2][0]<<","<<verts[2][1]<<","<<verts[2][2]<<">,"<<stick_rad<<endl;
                out << "pigment { color rgbt <"<<bo->GetBeginAtom ()->color[0]<<","<<bo->GetBeginAtom ()->color[1]<<","<<bo->GetBeginAtom ()->color[2]<<","<<1-bo->GetBeginAtom ()->color[3]<<"> }finish {ambient 0.1 diffuse 0.7 phong 0.3 phong_size 30 }"<<endl;
                out << "}"<<endl;

                out << "cylinder {"<<endl; 
                out << "<"<<verts[1][0]<<","<<verts[1][1]<<","<<verts[1][2]<<">,";
                out << "<"<<verts[3][0]<<","<<verts[3][1]<<","<<verts[3][2]<<">,"<<stick_rad<<endl;


                out << "pigment { color rgbt <"<<bo->GetBeginAtom ()->color[0]<<","<<bo->GetBeginAtom ()->color[1]<<","<<bo->GetBeginAtom ()->color[2]<<","<<1-bo->GetBeginAtom ()->color[3]<<"> }finish {ambient 0.1 diffuse 0.7 phong 0.3 phong_size 30 }"<<endl;


                out << "}"<<endl;     */
            }
      //      else cout<<"order"<<endl;
        }
    }

    }
    return out.str ();
}

string DDWin::POV_stick (vect v1, vect v2, color c1, color c2, float rad, float a0size, float a1size) {
	stringstream out;
	out << "merge {"<<endl;
	out << "    sphere {";
	out << POV_vector (v1)<<","<<rad<<endl;
	out << "pigment { "<<POV_color (c1)<<"} finish {ambient 0.1 diffuse 0.7 phong 0.3 phong_size 30 }}"<<endl;
	
	out << "    sphere {";
	out << POV_vector (v2)<<","<<rad<<endl;
	out << "pigment { "<<POV_color (c2)<<"} finish {ambient 0.1 diffuse 0.7 phong 0.3 phong_size 30 }}"<<endl;
	out << "     cylinder {"<<endl; 
	out <<POV_vector (v1)<<","<<POV_vector (v2)<<","<<rad<<endl;

	float gradx, grady, gradz, dist;
	gradx = v2.x() - v1.x();
	grady = v2.y() - v1.y();
	gradz = v2.z() - v1.z();
	dist = sqrt (gradx*gradx + grady*grady + gradz*gradz);
	out << "pigment {gradient <"<< gradx/dist<<","<< grady/dist<<","<< gradz/dist<<">\n  color_map {["<<a0size/dist<<" "<<POV_color (c1)<<"] ["<<1-a1size/dist<<" "<<POV_color (c2)<<"]} scale "<<dist<<" translate "<< POV_vector (v1)<<"}finish {ambient 0.1 diffuse 0.7 phong 0.3 phong_size 30 }}"<<endl;
	
	out << "no_shadow" << endl;
	out << "}"<<endl;
	return out.str ();
	
} 


string DDWin::POV_color (color c) {
    stringstream out;
    float r = c.redF ();
    float g = c.greenF ();
    float b = c.blueF ();
    float a = 1.f - c.alphaF ();
    out <<"color rgbt <"<<r<< ","<<g<<","<<b<<","<<a<<">";
    return out.str ();
}


string DDWin::POV_vector (vect v) {
    stringstream out;
    float x = v.x();
    float y = v.y();
    float z = v.z();
    out << "<"<<x<<","<<y<<","<<z<<"> ";
    return out.str ();
    
}


string DDWin::POV_atom (Atom *at){
    stringstream out;
    bool vis = get_visible (at);
    int ds = get_ds (at);
    double vdw = get_vdw (at);
    color col = get_color (at);
    if (ds != NO_ATOMS && vis) {
         float rad;
         if (ds == CPK_SPHERES) rad = vdw;
         else if (ds == SCALED_CPK_SPHERES) rad = vdw*gl->vdw_scale;
         else rad = 0;
		vect v = get_coordinates(at);
         float x = v.x();
         float y = v.y();
         float z = v.z();
         out << "sphere {"<<endl;
         out << "<"<<x<<","<<y<<","<<z<<">,"<<rad<<endl;
         out << "pigment { "<< POV_color (col)<<"} finish {ambient 0.1 diffuse 0.7 phong 0.3 phong_size 30}"<<endl;
         out << "}"<<endl;
     
    }
    return out.str ();              
}


void DDWin::write_POV_source (string filename) {
    ofstream *file = new ofstream(filename.c_str());
    GLint viewport [4];


    glGetIntegerv(GL_VIEWPORT, viewport);
    float width = (float)viewport[2];
    float height = (float)viewport [3];
    *file << "//Generated by Zodiac version " <<VERSION<< endl<<endl;

    *file << "#declare Move = transform {matrix <" << endl;

 			GLdouble m[16];	
			glGetDoublev(GL_MODELVIEW_MATRIX, m);
			
			double norm = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
			*file << "\t\t" << m[0] / norm << ",  " << m[1] / norm << ", " << m[2] / norm << "," << endl;
	//		m_[0] = m[0] / norm;
	//		m_[1] = m[1] / norm;
	//		m_[2] = m[2] / norm;
			norm = sqrt(m[4] * m[4] + m[5] * m[5] + m[6] * m[6]);
			*file << "\t\t" << m[4] / norm << ",  " << m[5] / norm << ", " << m[6] / norm << "," << endl;

//			m_[3] = m[4] / norm;
//			m_[4] = m[5] / norm;
//			m_[5] = m[6] / norm;
			norm = sqrt(m[8] * m[8] + m[9] * m[9] + m[10] * m[10]);
			*file << "\t\t" << m[8] / norm << ",  " << m[9] / norm << ", " << m[10] / norm << "," << endl;

//			m_[6] = m[8] / norm;
//			m_[7] = m[9] / norm;
//			m_[8] = m[10] / norm;
			*file << "\t\t" << m[12] << ",  " << m[13] << ", " << m[14] << endl;
//			m_[9] = m[12];
//			m_[10] = m[13];
//			m_[11] = m[14];
		    *file << "\t\t>" << endl;
			*file << "\tinverse }" << endl;


    *file << "#declare Backbone_pigment = pigment{ color rgbt <0.2, 0.0, 0.8, 0>}";
    *file <<"#declare Backbone_finish = finish {ambient 0.1 diffuse 0.7 phong 0.3 phong_size 30 specular 0.8 crand 0.1}"<< endl;
    *file << "#declare Backbone_rad = 0.1;"<<endl;




	*file << "camera {" << endl
	      << "\tperspective" << endl
		  << "\tdirection <0.0, 0.0, -1.0>" << endl
		  << "\tright " << width / height << " * x" << endl
		  << "\tangle 83.0" << endl
		  << "\ttransform {Move}}" << endl;


    *file << "global_settings { ambient_light rgb <1,1,1> ";
//    *file << "radiosity { brightness 0.6 }";
    *file <<"}"<<endl<<endl;
    *file << "background {"<< POV_color (gl -> background_color)<<"}"<<endl;
  //  *file << "camera {  location <"<<camerax<<","<<cameray<<","<<cameraz<<"> look_at <"<<lookx<<","<<looky<<","<<lookz<<"> up <"<<upx<<","<<upy<<","<<upz<<"> right <"<<rightx<<","<<righty<<","<<rightz<<">}"<<endl;
    *file << "light_source {<0, 0, +30> color rgbt <1,1,1,0.3> transform {Move}}"<<endl<<endl;

    for (unsigned int i=0;i<molecules.size(); i++){
        Molecule *mol = molecules[i];
        FOR_BONDS_OF_MOL (b, mol) {
            *file <<POV_bond (&*b);
        }
        FOR_ATOMS_OF_MOL (a, mol) {
            *file <<POV_atom (&*a);        
        }
 /*   vector<OBRing*>::iterator i;
    vector<OBRing*> *rlist = (vector<OBRing*>*)mol -> GetData("RingList");
    for (i = rlist->begin();i != rlist->end();++i) {

            if (gl->aromatic_display_style == AROMATIC_RINGS && (*i) -> IsAromatic ()) {
                *file <<POV_ring ((*i));        
            }
        }
*/

  /*      for (unsigned int r=0; r<mol->residues.size(); r++){
            if (mol->residues[r]->backbone_list.size()>1) {
                Point a = mol->residues[r]->backbone_list[0];
                *file << "merge {"<<endl;
                *file << "sphere {<"<<a.x()<<","<<a.y() <<","<<a.z() <<">,Backbone_rad}"<<endl;
                for (unsigned int p=0; p<mol->residues[r]->backbone_list.size ()-1; p++) {
                    Point b, c;
                    b = mol->residues[r]->backbone_list[p];
                    c = mol->residues[r]->backbone_list[p+1];
                    *file << "cylinder {<"<<b.x()<<","<<b.y() <<","<<b.z() <<">,<"<<c.x()<<","<<c.y() <<","<<c.z() <<">,Backbone_rad}"<<endl;
                    *file << "sphere {<"<<c.x()<<","<<c.y() <<","<<c.z() <<">,Backbone_rad}"<<endl;
                }
                *file << "pigment { Backbone_pigment}"<<endl;
                *file << "finish { Backbone_finish}"<<endl;
                *file << "}"<<endl;
            }
        } */
    } 
    for (unsigned int i=0;i<graphical_objects.size(); i++){
        if (graphical_objects[i] -> is_surface ()) {
            *file << POV_surface ((Surface *) graphical_objects[i]);
        }
        else if  (graphical_objects[i] -> is_sphere ()) {
            *file << POV_sphere ((Sphere *) graphical_objects[i]);
        }
    }
    file->close ();

}

string DDWin::POV_ring (Ring *ring) {

    return 0;

/*
    stringstream out;
    vect center, pointO;
    center = ring -> center;
    pointO.null ();

    float rad = dist (ring->bonds[0]->GetBeginAtom ()->GetVector (), center) * cos (PI/ring->bonds.size ()) - gl -> aromatic_bond_inter_distance * 2;

        float color_a, angle_b;
        bool invert_c = false;
    vect pointA, pointB, pointC;
    pointC = ring -> center;
    pointA = subtract (ring->bonds[0]->GetBeginAtom ()->GetVector (), pointC);
    pointB = subtract (ring->bonds[0]->GetEndAtom ()->GetVector (), pointC);
    pointC = ring -> center;

 //   glTranslatef (pointC[0], pointC[1], pointC[2]);
    float A, B, C;

    Atom *next_atom, *last_atom;
    Bond *next_b, *last_b;
    float tot_a, rotaxang=0.f, rangle=0.f;
        int side = -1;
    //equation of plane Ax + By + Cz = 0


    last_b = ring->bonds[0];
    last_atom = last_b->GetEndAtom ();
    tot_a = angle (ring->bonds[0]->GetBeginAtom ()->GetVector (), center, ring->bonds[0]->GetEndAtom ()->GetVector ());
    







    if (!pointA.z() && !pointB.z()) { //all three points on z=0 plane...  we're already there
        vect orang (1.f, 0.f, 0.f);
        if (pointA.y()>0 && pointB.y()>0) {
            color_a = angle (pointA, pointO, orang);
            angle_b = angle (pointB, pointO, orang);
            if (angle_b>color_a) invert_c = true;
            color_a*=-1;
        }
        else if (pointA.y()<0 && pointB.y()<0) {
            color_a = angle (pointA, pointO, orang);
            angle_b = angle (pointB, pointO, orang);
            if (angle_b<color_a) invert_c = true;
        }
        else if (pointA.y()>0 && pointB.y()<0) {
            color_a = -angle (pointA, pointO, orang);
            if (pointA.x()<0) invert_c = true;
        }
        else {
            color_a = angle (pointA, pointO, orang);                
            if (pointA.x()>0) invert_c = true;
        }

    }
    else if (pointA.x() && pointB.x() && (pointA.y()*pointB.x()-pointA.x()*pointB.y())){
    //equation (A/C)*x + (B/C)*y +z = 0
        C = 1;
        B = (pointA.x()*pointB.z()-pointB.x()*pointA.z())/(pointA.y()*pointB.x()-pointA.x()*pointB.y());
        A = (-pointA.z()-B*pointA.y())/pointA.x();

        //rotation axis A*x + B*y



        if (B) {
            vect fixpoint (1.f, -A/B, 0.f);


            if (pointA.y()>=-A/B*pointA.x() && pointB.y()>=-A/B*pointB.x()) {
                color_a = angle (pointA, pointO, fixpoint);
                angle_b = angle (pointB, pointO, fixpoint);
                if (angle_b>color_a) invert_c = true;
                color_a*=-1;
                cout <<"case 1"<<endl;
            }
            else if (pointA.y()<-A/B*pointA.x() && pointB.y()<-A/B*pointB.x()) {
                color_a = angle (pointA, pointO, fixpoint);
                angle_b = angle (pointB, pointO, fixpoint);
                if (angle_b<color_a) invert_c = true;
                cout <<"case 2"<<endl;
            }
            else if (pointA.y()>=-A/B*pointA.x() && pointB.y()<-A/B*pointB.x()) {
                color_a = -angle (pointA, pointO, fixpoint);
                cout <<"case 3"<<endl;
                if (pointA.y()<pointA.x()*B/(A+0.001)) invert_c = true;
            }
            else if (pointA.y()<-A/B*pointA.x() && pointB.y()>=-A/B*pointB.x()){
                color_a = angle (pointA, pointO, fixpoint);                
                if (pointA.y()>=pointA.x()*B/(A+0.001)) invert_c = true;
                cout <<"case 4"<<endl;
            }
                


            if ((pointA.z()>0 && pointA.y()>-A/B*pointA.x()) || (pointA.z()<0 && pointA.y()<-A/B*pointA.x())) side = -1;
            else side = 1;


            float d = sqrt (1+(-A/B)*(-A/B));
            rotaxang = asin ((-A/B)/d);
      /////      if (rotaxang<0) rotaxang= PI-rotaxang;

        }

        //perpendicular aA + bB = 0


        //perpendicular -B*x + A*y = 0
        //point x=1 y=B/A
        //point on the plane z= -(A + B*B/A)
        float xp, yp, zp;
        xp = 1;
        yp = B/A;
        zp = -(A + B*B/A);

        float dist =  sqrt(xp*xp+ yp*yp+zp*zp);

        rangle = asin (zp / dist);

        if (rangle<0) {
            rangle = -rangle;

        }






    }






    out <<"torus {\n";
    out <<rad<<", "<<gl->stick_rad/2<<endl;
    out<<"rotate <90,0,0>"<<endl;
    out<<"rotate <0,0,"<<-rotaxang*180/PI<<">"<<endl;

    out <<"pigment { radial \n  color_map { "<<endl;


    out <<"[ 0 "<<  POV_color (ring -> bonds[0] -> GetBeginAtom () -> col) <<"]"<<endl;

    out <<"[ "<<tot_a/360<<" "<<POV_color (ring -> bonds[0] -> GetEndAtom () ->col)<<"]"<<endl;

    for (unsigned int i=0; i<ring->bonds.size ()-2; i++) {
        for (unsigned int b=0; b<last_atom->bonds.size (); b++) {
            next_b = last_atom->bonds[b];
            if (next_b->is_in_ring (ring) && next_b!=last_b) break;
        }
        if (next_b->GetBeginAtom ()!= last_atom) next_atom = next_b->GetBeginAtom ();
        else  next_atom = next_b->GetEndAtom ();

        tot_a += angle (last_atom->GetVector (), center, next_atom->GetVector ());
        out <<"[ "<<tot_a/360<<" "<<POV_color (next_atom ->col)<<"]"<<endl;

        last_atom = next_atom;
    }

    out <<"[ 1 color rgb <1,1,0>]"<<endl;
    out << "}"<<endl;
    out<<"rotate <90,0,0>"<<endl;
    if (invert_c)     out<<"rotate <180,0,0>"<<endl;
    out<<"rotate <0,0,"<<-color_a <<">"<<endl;
    out << "}"<<endl;


    out<<"rotate <"<<-rangle*side*180/PI<<",0,0>"<<endl;
    out<<"rotate <0,0,"<<rotaxang*180/PI<<">"<<endl;

    out<<"translate < "<<ring->center.x()<<","<<ring->center.y()<<","<<ring->center.z()<<">"<<endl;
    out <<"}"<<endl<<endl;
    return out.str ();
*/
}



string DDWin::POV_surface (Surface *surf) {

    stringstream out;
    out << "mesh2 {\n" <<endl;
    out <<  "    vertex_vectors {"<<endl;
    out << surf->vertices.size ()<<", "<<endl;
    for (unsigned int i=0; i<surf->vertices.size (); i++) {
		
        out <<"<"<<surf->vertices[i]->GetVector ().x()<<","<<surf->vertices[i]->GetVector ().y()<<","<<surf->vertices[i]->GetVector ().z()<<">,"<<endl;
    }
    out <<"}"<<endl;
    out <<  "    normal_vectors {"<<endl;
    out << surf->vertices.size ()<<", "<<endl;
    for (unsigned int i=0; i<surf->vertices.size (); i++) {
        out <<"<"<<surf->vertices[i]->normal.x()<<","<<surf->vertices[i]->normal.y()<<","<<surf->vertices[i]->normal.z()<<">,"<<endl;
    }
    out <<"}"<<endl;
    out <<"    texture_list {"<<endl;
    out << surf->vertices.size ()<<", "<<endl;  
    for (unsigned int i=0; i<surf->vertices.size (); i++) {
        out <<"texture {pigment {"<<POV_color (surf -> vertices[i] ->col)<<endl;
        out<<"}finish {ambient 0.1 diffuse 0.7 phong 0.3 phong_size 30 }}"<<endl;
    }  
    out <<"}"<<endl;
    out <<  "    face_indices {"<<endl;
    out << surf->faces.size ()<<", "<<endl;
    for (unsigned int i=0; i<surf->faces.size (); i++) {
        out <<"<"<<surf->faces[i]->v1->n<<","<<surf->faces[i]->v2->n<<","<<surf->faces[i]->v3->n<<">,"<<surf->faces[i]->v1->n<<","<<surf->faces[i]->v2->n<<","<<surf->faces[i]->v3->n<<","<<endl;
    }
    out <<"}"<<endl;
    out <<"}"<<endl;
    return out.str ();
}

string DDWin::POV_sphere (Sphere *sphere) {

    stringstream out;
    out << "    sphere {";
    out << POV_vector (sphere -> center)<<","<< sphere->radius << endl;
    out << "pigment { "<<POV_color (sphere -> col)<<"} finish {ambient 0.1 diffuse 0.7 phong 0.3 phong_size 30 }}"<<endl;    


    return out.str ();
}

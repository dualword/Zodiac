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

#include "graphical_object.h"
#include "constants.h"
#include <math.h>
#include "ddwin.h"

using namespace std;

GraphicalObject::GraphicalObject () : name ("no name"){
    mesh = false;
	set_clippable (false);
	
}

bool GraphicalObject::is_surface () {
    return false;
}

bool GraphicalObject::is_sphere () {
    return false;
}

bool GraphicalObject::is_grid () {
    return false;
}

bool GraphicalObject::is_map () {
    return false;
}


Map::Map () : Surface () {
	cube = 0;
}


bool Map::is_map () {
	return true;
}

void Map::clean () {
	for (unsigned int i=0; i < faces.size (); i++) {
		delete faces[i];
	}
	faces.clear ();
	for (unsigned int i=0; i < vertices.size (); i++) {
		delete vertices[i];
	}
	vertices.clear ();
	
}

void Map::render () {
    glNewList(lst,GL_COMPILE);
    for (unsigned int i=0; i<faces.size (); i++) {
        glBegin (GL_TRIANGLES);
        SurfVertex *v1 = faces[i]->v1;
        SurfVertex *v2 = faces[i]->v2;
        SurfVertex *v3 = faces[i]->v3;
        set_color (v1 -> col);
        glNormal3f (v1->normal.x(), v1->normal.y(), v1->normal.z());
        glVertex3f (v1-> GetVector ().x(), v1-> GetVector ().y(), v1-> GetVector ().z());
        set_color (v2 -> col);
        glNormal3f (v2->normal.x(), v2->normal.y(), v2->normal.z());
        glVertex3f (v2-> GetVector ().x(), v2-> GetVector ().y(), v2-> GetVector ().z());
        set_color (v3 -> col);
        glNormal3f (v3->normal.x(), v3->normal.y(), v3->normal.z());
        glVertex3f (v3-> GetVector ().x(), v3-> GetVector ().y(), v3-> GetVector ().z());
        glEnd ();
    }
    glEndList ();
}

Surface::Surface () {
    mesh = false;
	near_to_dist = 4.f;
	near_to = 0;
	molecule = 0;
	grid = 0;
}

bool Surface::is_surface () {
    return true;
}

void Surface::set_molecule (ZNMolecule *mo, ZNMolecule *near_t) {
    molecule = mo;
	set_clippable (get_clippable (mo));   
    vector <Atom *> atoms;
    FOR_ATOMS_OF_MOL (a, mo) {
		if (CountBonds (&*a) < 4)   //to speed things up atoms at center of tetrhaedrons octaedrons and so on are not considered
            atoms.push_back (&*a);
    }        
	
    grid = new cutoffGrid<Atom*>(atoms, 4);
	
	if (near_t) {
		near_to = near_t;   
		vector <Atom *> atoms2;
		FOR_ATOMS_OF_MOL (a, near_t) {
			if (CountBonds (&*a) < 4)   //to speed things up atoms at center of tetrhaedrons octaedrons and so on are not considered
				atoms2.push_back (&*a);
		}        
		
		near_to_grid = new cutoffGrid<Atom*>(atoms2, near_to_dist);
		
	}
}




void Surface::render_as_surface () {
    glNewList(lst,GL_COMPILE);
    for (unsigned int i=0; i<faces.size (); i++) {
        glBegin (GL_TRIANGLES);
        SurfVertex *v1 = faces[i]->v1;
        SurfVertex *v2 = faces[i]->v2;
        SurfVertex *v3 = faces[i]->v3;
        set_color (v1 -> col);
        glNormal3f (v1->normal.x(), v1->normal.y(), v1->normal.z());
        glVertex3f (v1-> GetVector ().x(), v1-> GetVector ().y(), v1-> GetVector ().z());
        set_color (v2 -> col);
        glNormal3f (v2->normal.x(), v2->normal.y(), v2->normal.z());
        glVertex3f (v2-> GetVector ().x(), v2-> GetVector ().y(), v2-> GetVector ().z());
        set_color (v3 -> col);
        glNormal3f (v3->normal.x(), v3->normal.y(), v3->normal.z());
        glVertex3f (v3-> GetVector ().x(), v3-> GetVector ().y(), v3-> GetVector ().z());
        glEnd ();
    }
    glEndList ();
	
	
}

void Surface::render () {
    if (mesh) render_as_mesh ();
    else render_as_surface ();
}

void Surface::render_as_mesh () {
    glNewList(lst,GL_COMPILE);
    for (unsigned int i=0; i<faces.size (); i++) {
        glBegin (GL_LINES);
        SurfVertex *v1 = faces[i]->v1;
        SurfVertex *v2 = faces[i]->v2;
        SurfVertex *v3 = faces[i]->v3;
        set_color (v1 -> col);
        glNormal3f (v1->normal.x(), v1->normal.y(), v1->normal.z());
        glVertex3f (v1-> GetVector ().x(), v1-> GetVector ().y(), v1-> GetVector ().z());
        set_color (v2 -> col);
        glNormal3f (v2->normal.x(), v2->normal.y(), v2->normal.z());
        glVertex3f (v2-> GetVector ().x(), v2-> GetVector ().y(), v2-> GetVector ().z());
        set_color (v3 -> col);
        glNormal3f (v3->normal.x(), v3->normal.y(), v3->normal.z());
        glVertex3f (v3-> GetVector ().x(), v3-> GetVector ().y(), v3-> GetVector ().z());
		
        set_color (v1 -> col);
        glNormal3f (v1->normal.x(), v1->normal.y(), v1->normal.z());
        glVertex3f (v1-> GetVector ().x(), v1-> GetVector ().y(), v1-> GetVector ().z());
        set_color (v2 -> col);
        glNormal3f (v2->normal.x(), v2->normal.y(), v2->normal.z());
        glVertex3f (v2-> GetVector ().x(), v2-> GetVector ().y(), v2-> GetVector ().z());
        set_color (v3 -> col);
        glNormal3f (v3->normal.x(), v3->normal.y(), v3->normal.z());
        glVertex3f (v3-> GetVector ().x(), v3-> GetVector ().y(), v3-> GetVector ().z());
		
		
		
        glEnd ();
    }
    glEndList ();
	
}

void Surface::color_by_mol (ZNMolecule *mol, float a) {
	delete grid;
	vector <Atom *> atoms;
    FOR_ATOMS_OF_MOL (at, mol) {
		atoms.push_back (&*at);
    }        
	
    grid = new cutoffGrid<Atom*>(atoms, 4);
	color_by_atom (a);
	
}

void Surface::color_by_map (Map *map, color c1, color c2, color c3, float f1, float f2, float f3) {
	for (unsigned int ii=0; ii< faces.size (); ii++) { 
        color_vertex_by_map (faces[ii]->v1,map, c1, c2, c3, f1, f2, f3);   
        color_vertex_by_map (faces[ii]->v2,map, c1, c2, c3, f1, f2, f3);    
        color_vertex_by_map (faces[ii]->v3,map, c1, c2, c3, f1, f2, f3);    
		
    }
}

void Surface::alpha_by_mol_distance (ZNMolecule *mol, float distance) {
	if (grid) delete grid;
	vector <Atom *> atoms;
    FOR_ATOMS_OF_MOL (a, mol) {
		atoms.push_back (&*a);
    }        
	
    grid = new cutoffGrid<Atom*>(atoms, distance);
    for (unsigned int ii=0; ii< faces.size (); ii++) { 
		fade_vertex_alpha_by_atom (faces[ii]->v1, distance);
		fade_vertex_alpha_by_atom (faces[ii]->v2, distance);
		fade_vertex_alpha_by_atom (faces[ii]->v3, distance);
    }
	
}


void Surface::color_by_atom (float a) {
    for (unsigned int ii=0; ii< faces.size (); ii++) { 
        color_vertex_by_atom_smooth (faces[ii]->v1,a);   
        color_vertex_by_atom_smooth (faces[ii]->v2,a); 
        color_vertex_by_atom_smooth (faces[ii]->v3,a); 
		
    }
}

void Surface::color_by_potential (ZNMolecule *mol, color lipo, color acc, color don, float threshold) {
	
	Chemscore *chemscore = new Chemscore;
	ZNMolecule *mol2 = new ZNMolecule;
	vector <ZNMolecule *> mols;
	mols.push_back(mol);
	chemscore ->load_environment (mols, mol2);
	
	float rl = lipo.redF ();
	float gl = lipo.greenF ();
	float bl = lipo.blueF ();
	float al = lipo.alphaF ();
	
	float ra = acc.redF ();
	float ga = acc.greenF ();
	float ba = acc.blueF ();
	float aa = acc.alphaF ();
	
	float rd = don.redF ();
	float gd = don.greenF ();
	float bd = don.blueF ();
	float ad = don.alphaF ();
	float r, g, b, a, tot, total_score;
	for (unsigned int ii=0; ii< faces.size (); ii++) {
		float scorelipo1 = chemscore ->Livalue (faces[ii]->v1 ->coordinates);
		float scoreacc1 = chemscore ->Acceptorvalue (faces[ii]->v1->coordinates);
		float scoredon1 = chemscore ->Donorvalue (faces[ii]->v1->coordinates);
		tot = 1 / (scorelipo1 + scoreacc1 + scoredon1);
		total_score = scorelipo1 + scoreacc1 + scoredon1;
		scorelipo1 *= tot;
		scoreacc1 *= tot;
		scoredon1 *= tot;
		r = (rl * scorelipo1 + ra * scoreacc1 + rd * scoredon1);
		g = (gl * scorelipo1 + ga * scoreacc1 + gd * scoredon1) ;
		b = (bl * scorelipo1 + ba * scoreacc1 + bd * scoredon1);
		a = (al * scorelipo1 + aa * scoreacc1 + ad * scoredon1) ;
		if (total_score < threshold) {
			color_vertex_by_color(faces[ii] ->v1, color (r, g, b, a));
		}
		else {
			float k1 = total_score /threshold ;
			if (k1 < 0.f) k1 = 0.f;
			if (k1 > 1.f) k1 = 1.f;
			float k2 = 1.f - k1;
			float lastr = faces[ii] ->v1 ->col.redF ();
			float lastg = faces[ii] ->v1 ->col.greenF ();
			float lastb = faces[ii] ->v1 ->col.blueF ();
			float lasta = faces[ii] ->v1 ->col.alphaF ();
			color_vertex_by_color(faces[ii] ->v1, color (r*k1 + lastr*k2, g*k1 + lastg*k2, b*k1 + lastg*k2, a*k1 + lasta*k2));
			
		}
		
		float scorelipo2 = chemscore ->Livalue (faces[ii]->v2->coordinates);
		float scoreacc2 = chemscore ->Acceptorvalue (faces[ii]->v2->coordinates);
		float scoredon2 = chemscore ->Donorvalue (faces[ii]->v2->coordinates);
		tot = 1 / (scorelipo2 + scoreacc2 + scoredon2);
		total_score = scorelipo2 + scoreacc2 + scoredon2;
		scorelipo2 *= tot;
		scoreacc2 *= tot;
		scoredon2 *= tot;
		r = (rl * scorelipo2 + ra * scoreacc2 + rd * scoredon2);
		g = (gl * scorelipo2 + ga * scoreacc2 + gd * scoredon2);
		b = (bl * scorelipo2 + ba * scoreacc2 + bd * scoredon2) ;
		a = (al * scorelipo2 + aa * scoreacc2 + ad * scoredon2);
		
		if (total_score < threshold) {
			color_vertex_by_color(faces[ii] ->v2, color (r, g, b, a));
		}
		else {
			float k1 = total_score /threshold ;
			if (k1 < 0.f) k1 = 0.f;
			if (k1 > 1.f) k1 = 1.f;
			float k2 = 1.f - k1;
			float lastr = faces[ii] ->v2 ->col.redF ();
			float lastg = faces[ii] ->v2 ->col.greenF ();
			float lastb = faces[ii] ->v2 ->col.blueF ();
			float lasta = faces[ii] ->v2 ->col.alphaF ();
			color_vertex_by_color(faces[ii] ->v2, color (r*k1 + lastr*k2, g*k1 + lastg*k2, b*k1 + lastg*k2, a*k1 + lasta*k2));
			
		}		
		
		float scorelipo3 = chemscore ->Livalue (faces[ii]->v3->coordinates);
		float scoreacc3 = chemscore ->Acceptorvalue (faces[ii]->v3->coordinates);
		float scoredon3 = chemscore ->Donorvalue (faces[ii]->v3->coordinates);
		tot = 1 / (scorelipo3 + scoreacc3 + scoredon3);
		total_score = scorelipo3 + scoreacc3 + scoredon3;
		scorelipo3 *= tot;
		scoreacc3 *= tot;
		scoredon3 *= tot;
		r = (rl * scorelipo3 + ra * scoreacc3 + rd * scoredon3);
		g = (gl * scorelipo3 + ga * scoreacc3 + gd * scoredon3) ;
		b = (bl * scorelipo3 + ba * scoreacc3 + bd * scoredon3) ;
		a = (al * scorelipo3 + aa * scoreacc3 + ad * scoredon3);
		if (total_score < threshold) {
			color_vertex_by_color(faces[ii] ->v3, color (r, g, b, a));
		}
		else {
			float k1 = total_score /threshold ;
			if (k1 < 0.f) k1 = 0.f;
			if (k1 > 1.f) k1 = 1.f;
			float k2 = 1.f - k1;
			float lastr = faces[ii] ->v3 ->col.redF ();
			float lastg = faces[ii] ->v3 ->col.greenF ();
			float lastb = faces[ii] ->v3 ->col.blueF ();
			float lasta = faces[ii] ->v3 ->col.alphaF ();
			color_vertex_by_color(faces[ii] ->v3, color (r*k1 + lastr*k2, g*k1 + lastg*k2, b*k1 + lastg*k2, a*k1 + lasta*k2));
			
		}
		
    }
	delete mol2;
	delete chemscore;
}

void Surface::color_by_color (color c) {
    for (unsigned int ii=0; ii< faces.size (); ii++) { 
        color_vertex_by_color (faces[ii]->v1,c);   
        color_vertex_by_color (faces[ii]->v2,c); 
        color_vertex_by_color (faces[ii]->v3,c); 
		
    }	
}




void Surface::color_vertex_by_color (SurfVertex *vert, color c) {
	vert -> col = c;
}

void Surface::color_vertex_by_map (SurfVertex *vert, Map *map, color c1, color c2, color c3, float f1, float f2, float f3) {
	float score = map ->get_value (vert ->coordinates.x(), vert ->coordinates.y(), vert ->coordinates.z());
	cerr << "score" << score<<endl;
	color c = average_3_colors(score, c1, c2, c3, f1, f2, f3);
	vert ->col = c;
}

void Surface::color_vertex_by_atom (SurfVertex *vert, float alp) {
	
    objectList<Atom*>* nbAtoms = grid->getNeighborObjects(vert-> GetVector ());
    if (nbAtoms) {
        vector <Atom *> neighbours = nbAtoms->objects;
        float rr, gg, bb, aa;
		float dis = dist (get_coordinates (neighbours[0]), get_coordinates (vert));
		//  float dis = sqrt((neighbours[0]-> GetVector ().x()-vert-> GetVector ().x())*(neighbours[0]-> GetVector ().x()-vert-> GetVector ().x())+
		// (neighbours[0]-> GetVector ().y()-vert-> GetVector ().y())*(neighbours[0]-> GetVector ().y()-vert-> GetVector ().y())+
		//(neighbours[0]-> GetVector ().z()-vert-> GetVector ().z())*(neighbours[0]-> GetVector ().z()-vert-> GetVector ().z()));
        color col = get_color (neighbours[0]);
        rr = col.redF ();
        gg = col.greenF ();
        bb = col.blueF ();
        aa = col.alphaF ();
		
        for (unsigned int i=0; i<neighbours.size (); i++) {
			float n_dist = dist (get_coordinates (neighbours[i]), get_coordinates (vert));
			//          float n_dist = sqrt((neighbours[i]-> GetVector ().x()-vert-> GetVector ().x())*(neighbours[i]-> GetVector ().x()-vert-> GetVector ().x())+
			//  (neighbours[i]-> GetVector ().y()-vert-> GetVector ().y())*(neighbours[i]-> GetVector ().y()-vert-> GetVector ().y())+
			//(neighbours[i]-> GetVector ().z()-vert-> GetVector ().z())*(neighbours[i]-> GetVector ().z()-vert-> GetVector ().z()));
            if (n_dist < dis) {
                dis = n_dist;
                color coli = get_color (neighbours[i]);
                rr = coli.redF ();    
                gg = coli.greenF ();
                bb = coli.blueF ();
                aa = coli.alphaF ();
            }
            
        }
        vert -> col = color (rr, bb, gg, aa*alp);
    }
}

void Surface::set_alpha (double d) {
	for (unsigned int ii=0; ii< faces.size (); ii++) { 
        set_vertex_alpha (faces[ii]->v1,d);   
        set_vertex_alpha (faces[ii]->v2,d); 
        set_vertex_alpha (faces[ii]->v3,d); 
		
    }	
}

void Surface::multiply_alpha (double d) {
	for (unsigned int ii=0; ii< faces.size (); ii++) { 
        multiply_vertex_alpha (faces[ii]->v1,d);   
        multiply_vertex_alpha (faces[ii]->v2,d); 
        multiply_vertex_alpha (faces[ii]->v3,d); 
		
    }	
}

void Surface::set_vertex_alpha (SurfVertex *vert, float a) {
	vert -> col.setAlphaF (a);
}

void Surface::multiply_vertex_alpha (SurfVertex *vert, float a) {
	float f = vert ->col.alphaF ();

	vert -> col.setAlphaF (a*f);
}


void Surface::fade_vertex_alpha_by_atom (SurfVertex *vert, float distance) {
	objectList<Atom*>* nbAtoms = grid->getNeighborObjects(vert-> GetVector ());
    if (nbAtoms) {
        vector <Atom *> neighbours = nbAtoms->objects;
		float dis = BIG;
		for (unsigned int i=0; i<neighbours.size (); i++) {
			float d = dist (get_coordinates (neighbours[i]), get_coordinates (vert));
			if (d < dis) dis = d;
        }
		float perc = dis / distance;
		float alpha = 1.f;
		if (perc < 0.4f) alpha = 1.f;
		else if (perc >0.7f) alpha = 0.f;
		else alpha = 1 -((perc - 0.4f) / (0.7f - 0.4f));
		vert ->col.setAlphaF (alpha);	
	}
	else vert ->col.setAlphaF (0.f);
}

void Surface::color_vertex_by_atom_smooth (SurfVertex *vert, float alp) { //blends 2 closest atoms color, or closest atom with previous color
	objectList<Atom*>* nbAtoms = grid->getNeighborObjects(vert-> GetVector ());
    if (nbAtoms) {
        vector <Atom *> neighbours = nbAtoms->objects;
        float r, g, b, a;
        r = g = b = a = 0.f;
		float total_val = 0.f;
		float threshold = 1.f;
		
        for (unsigned int i=0; i<neighbours.size (); i++) {
			float dis = dist (get_coordinates (neighbours[i]), get_coordinates (vert));
			float val = (1 - (dis / 4.f));
			if (val < 0) continue;
			else {
				color col = get_color (neighbours[i]);
				r += col.redF()*val;
				g += col.greenF()*val;
				b += col.blueF()*val;
				a += col.alphaF()*val;
				
				total_val += val;
			}
			
			
			
		}
		if (total_val > threshold) {
			r /= total_val;
			g /= total_val;
			b /= total_val;
			a /= total_val;
		}
		else {
			r += vert ->col.redF () * (threshold - total_val);
			g += vert ->col.greenF () * (threshold - total_val);
			b += vert ->col.blueF () * (threshold - total_val);
			a += vert ->col.alphaF () * (threshold - total_val);
			r /= threshold;
			g /= threshold;
			b /= threshold;
			a /= threshold;
		}
		vert -> col = color (r, g, b, a);
	}
	
	
	/*   objectList<Atom*>* nbAtoms = grid->getNeighborObjects(vert-> GetVector ());
	 if (nbAtoms) {
	 vector <Atom *> neighbours = nbAtoms->objects;
	 float r1, g1, b1, a1;
	 float r2, g2, b2, a2;
	 r1 = g1 = b1 = a1 = r2 = g2 = b2 = a2 = 0.f;
	 float dist1 = BIG;
	 float dist2 = BIG;
	 int idx1 = -1;
	 int idx2 = -1;
	 
	 for (unsigned int i=0; i<neighbours.size (); i++) {
	 float dis = dist (get_coordinates (neighbours[i]), get_coordinates (vert));
	 if (dis < dist1) {
	 dist1 = dis;
	 idx1 = i;
	 }
	 else if (dis < dist2) {
	 dist2 = dis;
	 idx2 = i;
	 }
	 
	 
	 }
	 
	 
	 
	 if (dist1 < BIG) {
	 color col = get_color (neighbours[idx1]);
	 float val = 1/ dist1;
	 //	float val = 1/(dist1*dist1*dist1);
	 r1 = col.redF()*val;
	 g1 = col.greenF()*val;
	 b1 = col.blueF ()*val;
	 a1 = col.alphaF ()*val;
	 if (dist2 < BIG) {
	 color col2 = get_color (neighbours[idx2]);
	 //		float val2 = 1/(dist1*dist1*dist1);
	 float val2 = 1/dist2;
	 r2 = col2.redF()*val2;
	 g2 = col2.greenF()*val2;
	 b2 = col2.blueF ()*val2;
	 a2 = col2.alphaF ()*val2;
	 float tv = val + val2;
	 vert -> col = color ((r1+r2)/tv, (g1+g2)/tv,(b1+b2)/tv,(a1+a2)/tv);
	 
	 }
	 else {
	 color col2 = vert ->col;
	 float val2 = 1-val;
	 r2 = col.redF()*val2;
	 g2 = col.greenF()*val2;
	 b2 = col.blueF ()*val2;
	 a2 = col.alphaF ()*val2;
	 }
	 
	 }
	 else {
	 // do nothing
	 }
	 
	 
	 
	 }
	 */
}


bool Sphere::is_sphere () {
    return true;
}

void Sphere::render () {
    render_as_surface ();
}

void Sphere::render_as_surface () {
    glNewList(lst,GL_COMPILE);
    glPushMatrix ();
    glTranslatef (center.x(), center.y(), center.z());
    glColor4f (col.redF (), col.greenF(), col.blueF(), col.alphaF());
    gluSphere (gluNewQuadric(), radius, 30, 30);
    glPopMatrix ();
    glEndList ();
	
}





Grid::Grid () {
    cube = new MarchingCubes (); 
    res = 1.f;
    threshold = 1.f; 
    ff = new MMFF ();
	//    function = &null;    
}

float Grid::null (float x, float y, float z) {
    cerr << "function not defined" << endl;
    return 0.f;
}

bool Grid::is_grid () {
    return true;
}

void Grid::init_ff () {
    probe = new ZNMolecule ();
    Atom *at = new Atom ();
    at -> SetAtomicNum (6);
	//    at -> getVdw ();
    at -> SetPartialCharge (0);
	//   at -> MMFFstring = "CR";
	//   at -> MMFFtype = 1;
    probe->ZNAddAtom (at);
    ff -> clear_nonbonded_interactions ();
    ff -> initialize_interaction (probe, env);  
}


void Grid::load () {
	
    int xm = (int)((end.x() - origin.x())* res);
    int ym = (int)((end.y() - origin.y())* res);
    int zm = (int)((end.z() - origin.z())* res);
	//   cout <<xm<<" "<<ym<<" "<<zm<<endl;
	
	
    cube->set_resolution( xm, ym, zm) ;
    cube->set_method (false); //use original MC algo?
    cube->set_limits (origin.x(), origin.y(), origin.z(), end.x(), end.y(), end.z());
    cube->init_all() ;
    float x, y, z;
	
    for (unsigned int k=0; k<zm; k++) {
        z =  cube->to_real_z (k);
        for (unsigned int j=0; j<ym; j++) {
            y =  cube->to_real_y (j);
            for (unsigned int i=0; i<xm; i++) {
                x = cube->to_real_x (i);
                float value = function (x, y, z);
                cube->set_data(value, i, j, k );
            }    
        }
    }
    cube->run(threshold) ;
}



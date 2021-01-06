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
#include <QtOpenGL/qgl.h>
#include "constants.h"
#include <math.h>
#include "ddwin.h"




using namespace std;





GraphicalObject::GraphicalObject () {
    mesh = false;

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


Surface::Surface () {
    mesh = false;
}

bool Surface::is_surface () {
    return true;
}

void Surface::set_molecule (Molecule *mo) {
    molecule = mo;   
    vector <Atom *> atoms;
    FOR_ATOMS_OF_MOL (a, mo) {
        Atom * at = &*a;
        atoms.push_back (at);
    }        

    grid = new cutoffGrid<Atom*>(atoms, 4);};


void Surface::render_as_surface () {

    glNewList(list,GL_COMPILE);
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
    glNewList(list,GL_COMPILE);
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

void Surface::color_by_atom (float a) {
    for (unsigned int ii=0; ii< faces.size (); ii++) { 
        color_vertex_by_atom_smooth (faces[ii]->v1,a);   
        color_vertex_by_atom_smooth (faces[ii]->v2,a); 
        color_vertex_by_atom_smooth (faces[ii]->v3,a); 

    }
}




void Surface::color_vertex_by_atom (SurfVertex *vert, float alp) {

    objectList<Atom*>* nbAtoms = grid->getNeighborObjects(vert-> GetVector ());
    if (nbAtoms) {
        vector <Atom *> neighbours = nbAtoms->objects;
        float rr, gg, bb, aa;
        float dist = sqrt((neighbours[0]-> GetVector ().x()-vert-> GetVector ().x())*(neighbours[0]-> GetVector ().x()-vert-> GetVector ().x())+
    (neighbours[0]-> GetVector ().y()-vert-> GetVector ().y())*(neighbours[0]-> GetVector ().y()-vert-> GetVector ().y())+
(neighbours[0]-> GetVector ().z()-vert-> GetVector ().z())*(neighbours[0]-> GetVector ().z()-vert-> GetVector ().z()));
        color col = get_color (neighbours[0]);
        rr = col.redF ();
        gg = col.greenF ();
        bb = col.blueF ();
        aa = col.alphaF ();

        for (unsigned int i=0; i<neighbours.size (); i++) {
            float n_dist = sqrt((neighbours[i]-> GetVector ().x()-vert-> GetVector ().x())*(neighbours[i]-> GetVector ().x()-vert-> GetVector ().x())+
    (neighbours[i]-> GetVector ().y()-vert-> GetVector ().y())*(neighbours[i]-> GetVector ().y()-vert-> GetVector ().y())+
(neighbours[i]-> GetVector ().z()-vert-> GetVector ().z())*(neighbours[i]-> GetVector ().z()-vert-> GetVector ().z()));
            if (n_dist < dist) {
                dist = n_dist;
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

void Surface::color_vertex_by_atom_smooth (SurfVertex *vert, float alp) {

    objectList<Atom*>* nbAtoms = grid->getNeighborObjects(vert-> GetVector ());
    if (nbAtoms) {
        vector <Atom *> neighbours = nbAtoms->objects;
        float r, g, b, a;
        r = g = b = a = 0.f;
        float tot_val = 0.f;
        for (unsigned int i=0; i<neighbours.size (); i++) {
            float dist = sqrt((neighbours[i]-> GetVector ().x()-vert-> GetVector ().x())*(neighbours[i]-> GetVector ().x()-vert-> GetVector ().x())+
    (neighbours[i]-> GetVector ().y()-vert-> GetVector ().y())*(neighbours[i]-> GetVector ().y()-vert-> GetVector ().y())+
(neighbours[i]-> GetVector ().z()-vert-> GetVector ().z())*(neighbours[i]-> GetVector ().z()-vert-> GetVector ().z()));
            float val =1/(dist*dist*dist);
            color col = get_color (neighbours[i]);
            r += col.redF()*val;
            g += col.greenF()*val;
            b += col.blueF ()*val;
            a += col.alphaF ()*val;
            tot_val+=val;
            
        }
        vert -> col = color (r/tot_val, g/tot_val,b/tot_val, a/tot_val * alp);
    }
}


bool Sphere::is_sphere () {
    return true;
}

void Sphere::render () {
    render_as_surface ();
}

void Sphere::render_as_surface () {
    glNewList(list,GL_COMPILE);
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
    treshold = 1.f; 
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
    probe = new Molecule ();
    Atom *at = new Atom ();
    at -> SetAtomicNum (6);
//    at -> getVdw ();
    at -> SetPartialCharge (0);
 //   at -> MMFFstring = "CR";
 //   at -> MMFFtype = 1;
    probe -> AddAtom (*at);
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
    cube->run(treshold) ;
}



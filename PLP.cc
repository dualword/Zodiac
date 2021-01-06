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
#include "PLP.h"

#define DONOR 0
#define ACCEPTOR 1
#define BOTH 2
#define NONPOLAR 3

#define STERIC 0
#define HBOND 1
#define REPULSIVE 2


float PLPInteraction::value () {
    float r = dist (at1-> GetVector (), at2-> GetVector ());
    float a, b, c, d, e, f;
    float hb = 1;
    if (type == STERIC) {
        
        b = etab.GetVdwRad(at1 -> GetAtomicNum ()) + etab.GetVdwRad(at2 -> GetAtomicNum ());
        a = 0.93f * b;
        c = 1.25f * b;
        d = 1.5f * b;
        e = -0.4f;
        f = 15.f;
    }
    else if (type == HBOND) {
        double ang = angle (root, at2-> GetVector (), at1-> GetVector ());
     //   if (ang>PI) ang = 2*PI-ang;
     //   cout <<ang;
        if (ang < 90.f) hb = 0.f;
        else if (ang < 120.f) hb = (ang-90.f)/(30.f);
        else hb = 1.f;
   //     cout <<"hb "<<hb<<endl;
        a = 2.3f;
        b = 2.6f;
        c = 3.1f;
        d = 3.4f;
        e = -2.f;
        f = 15.f;
    }
    else if (type == REPULSIVE) {
        a = 3.4f;
        e = 0.f;
        f = 15.f;
    }
    if (r < a) return f*(a-r)/a*hb;    
    else if (type != REPULSIVE) {
        if (r<b) return e*(r-a)/(b-a)*hb;
        else if (r<c) return e*hb;
        else if (r<d) return e*(d-r)/(d-c)*hb;
        else return 0.f;
    }
    else return 0.f;
}
    




PLP::PLP () {
    is_initialised = TRUE;
}



void PLP::clear_nonbonded_interactions () {
    for (unsigned int i=0; i<NBInteractions.size (); i++) delete NBInteractions[i];
    NBInteractions.clear ();
}





void PLP::compute_forces () {
    for (unsigned int i=0; i<NBInteractions.size (); i++) {
        NBInteractions[i] -> set_forces ();
    }
}

/*
void PLP::compute_force (PLPInteraction *plpint) {
    for (unsigned int i=0; i<3; i++) {
        plpint->at1-> GetVector ()[i]-= DX;
        float E = compute_interaction (plpint);
        plpint->at1-> GetVector ()[i]+= 2*DX;
        float E2 = compute_interaction (plpint);
        plpint->at1-> GetVector ()[i]-= DX;
        plpint->at1->force[i]-=  (E2-E)/(2*DX); //3* force
        plpint->at2->force[i]+=  (E2-E)/(2*DX);

    }
}

float PLP::compute_interaction (PLPInteraction *plpint) {
    float r = distance (plpint->at1-> GetVector (), plpint->at2-> GetVector ());
    float a, b, c, d, e, f;
    float hb = 1;
    if (plpint->type == STERIC) {
        b = plpint->at1->vdw + plpint->at2->vdw;
        a = 0.93*b;
        c = 1.25*b;
        d = 1.5*b;
        e = -0.4;
        f = 15;
    }
    else if (plpint->type == HBOND) {
        double ang = angle (plpint->root, plpint->at2-> GetVector (), plpint->at1-> GetVector ());
     //   if (ang>PI) ang = 2*PI-ang;
     //   cout <<ang;
        if (ang<90) hb = 0;
        else if (ang<120) hb = (ang-90)/(30);
        else hb = 1;
   //     cout <<"hb "<<hb<<endl;
        a = 2.3;
        b = 2.6;
        c = 3.1;
        d = 3.4;
        e = -2;
        f = 15;
    }
    else if (plpint->type == REPULSIVE) {
        a = 3.4;
        e = 0;
        f = 15;
    }
    if (r < a) return f*(a-r)/a*hb;    
    else if (plpint->type != REPULSIVE) {
        if (r<b) return e*(r-a)/(b-a)*hb;
        else if (r<c) return e*hb;
        else if (r<d) return e*(d-r)/(d-c)*hb;
        else return 0;
    }
    else return 0;
}

*/


void PLP::load_internal_interactions () {}

void PLP::load_nonbonded_interactions () {

    clear_nonbonded_interactions ();
    Molecule *mol = target_mol;
    FOR_ATOMS_OF_MOL (a, mol) {
        if (a -> IsHydrogen ()) continue;
        objectList<Atom*>* nbAtoms = far_grid->getNeighborObjects(a -> GetVector ());
        if (nbAtoms) {
            vector <Atom *> neighbours = nbAtoms->objects;      
            for (unsigned int j=0; j<neighbours.size (); j++) {
                if (neighbours[j]-> IsHydrogen ()) continue;
                float x = 0.f, y = 0.f, z = 0.f;
                if (neighbours[j]-> GetAtomicNum () == 7 || neighbours[j]-> GetAtomicNum () == 8) {
                    FOR_NBORS_OF_ATOM (n, neighbours[j]) {
                        if (n -> IsHydrogen ()) {
                            x += n -> GetVector ().x();
                            y += n -> GetVector ().y();
                            z += n -> GetVector ().z();
                        }
                    }
                    x/= CountBonds (neighbours[j]);
                    y/= CountBonds (neighbours[j]);
                    z/= CountBonds (neighbours[j]);
                }
                PLPInteraction *plpint = new PLPInteraction;
                plpint -> at1 = &*a;
                plpint -> at2 = neighbours[j];
                plpint -> I = getPLPtype (plpint->at1);
                plpint -> J = getPLPtype (plpint->at2);
                plpint -> type = getPLPinteractiontype (plpint->I, plpint->J);
                plpint -> root.x() = x;
                plpint -> root.y() = y;
                plpint -> root.z() = z;
                NBInteractions.push_back (plpint);    
            }        
        }
    }
}


int PLP::getPLPtype (Atom *at) {
    Molecule *mol = (Molecule *) at -> GetParent ();
    if (at->IsOxygen ()) {
        if (mol->bonded_to (at, -2, 1)) return BOTH;
        else return ACCEPTOR;
    }
    else if (at -> IsNitrogen () == 7) {
        if (mol->bonded_to (at, -2, 1)) return DONOR;
        else return ACCEPTOR;
    }
    else return NONPOLAR;
}

 
int PLP::getPLPinteractiontype (int I, int J) {
    if (I==NONPOLAR || J==NONPOLAR) return STERIC;
    else if (I == BOTH || J ==BOTH)  {return HBOND;}
    else if (I != J) { return HBOND;}
    else return REPULSIVE;
}

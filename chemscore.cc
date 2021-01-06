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

#include "chemscore.h"
#include "math.h"

#define DONOR 0
#define ACCEPTOR 1
#define BOTH 2
#define POLAR 3
#define NONPOLAR 4
#define METAL 5

using namespace OpenBabel;

float f (float x, float x1, float x2) {
    if (x<= x1) return 1;
    else if (x>x2) return 0;
    else return (x2-x)/(x2-x1);    

}


float ChemscoreHBInteraction::value () {
    float r = dist (at1-> GetVector (), at2-> GetVector ());
    float ang;
    if (at1->GetAtomicNum () == 1) ang = angle (at2-> GetVector (), at1-> GetVector (), root-> GetVector ());
    else ang = angle (at1-> GetVector (), at2-> GetVector (), root-> GetVector ());
    float dghbond = -3.34f;
    float r0 = 1.85f;
    float a0 = 180.f;
    float drab = r-r0;
    if (drab < 0) drab = -drab;
    float daab = ang-a0;
    if (daab < 0) daab = -daab;
    float dr1 = 0.25f;
    float dr2 = 0.65f;
    float da1 = 30.f;
    float da2 = 80.f;

    return dghbond * f(drab, dr1, dr2) * f(daab, da1, da2);
}


float ChemscoreLiInteraction::value () {
   float r1, r2, dg;
    if (metal) {
        dg = -6.03f;
        r1 = 2.6f;
        r2 = 3.0f;
    }
    else {
        dg = -0.117f;
        r1 = 4.1f;
        r2 = 7.1f;
    }
    float r = dist (at1-> GetVector (), at2-> GetVector ());

    return dg * f(r, r1, r2);

}




float ChemscoreClInteraction::value () {
    float r = dist (at1-> GetVector (), at2-> GetVector ());
    if (type ==0) {
        float rcl = 1.6f;
        float dghb = -3.34f;
        if (r>rcl) return 0.f;
        else return (20.f/dghb)*(rcl-r)/rcl;
    }
    else if (type ==1) {
        float rcl = 1.38f;
        float dgm = -6.03f;
        if (r>rcl) return 0.f;
        else return (20.f/dgm)*(rcl-r)/rcl;

    }
    else {
        float rcl;
        if (at1->GetAtomicNum () == 16 || at2->GetAtomicNum () == 16) rcl=3.35f;
        else rcl = 3.1f; 
        if (r>rcl) return 0.f;
        else return 1.f+4.f*(rcl-r)/rcl-1.f;
    }
}




Chemscore::Chemscore () {
    is_initialised = TRUE;
}

int Chemscore::getChemscoretype (Atom *at) {
    Molecule *mol = (Molecule *) at -> GetParent ();

    if (at->GetAtomicNum () == 8) { 
        return ACCEPTOR;
    }
    if (at->GetAtomicNum () == 7) {
        if (CountBonds (at)<3) return ACCEPTOR;
        else return POLAR;
    }
    if (at->GetAtomicNum () == 1 && (mol -> bonded_to (at, -2, 7) || mol -> bonded_to (at, -2, 8))) return DONOR;
    if (at->GetAtomicNum () == 15) return POLAR;
    if (at->GetAtomicNum () == 9) return POLAR;
    if (at->GetAtomicNum () == 16 && (mol->bonded_to (at, -2,7)+mol ->bonded_to (at, -2,8))) return POLAR;
    
    if (at->GetAtomicNum () == 6) {
        if (mol ->bonded_to (at, 3, 7)) return POLAR;
        else if (mol ->bonded_to (at, 2, 8)) return POLAR;
        else if ((mol ->bonded_to (at, -2,7)+mol ->bonded_to (at, -2,8))>2) return POLAR;
        else return NONPOLAR;
    } 
    if (at->GetAtomicNum () == 26) return METAL;  //iron other metals are missing
    return NONPOLAR;








/*  implicit H parametrization
    if (at->GetAtomicNum () == 8) { 
        if (mol ->bonded_to (at, -2, 1)) return BOTH;
        else return ACCEPTOR;
    }
    if (at->GetAtomicNum () == 7) {
        if (mol ->bonded_to (at, -2, 1) ) && (mol ->bonded_to (at, 2, 6) )) return BOTH; //imine
        else if (mol ->bonded_to (at, -2, 1) )) return DONOR; 
        else if (bound.size ()<3) return ACCCEPTOR;
        else return POLAR;
    }
    if (at->GetAtomicNum () == 1 && (mol ->bonded_to (at, -2, 7) || mol ->bonded_to (at, -2, 8))) return DONOR;
    if (at->GetAtomicNum () == 15) return POLAR;
    if (at->GetAtomicNum () == 9) return POLAR;
    if (at->GetAtomicNum () == 16 && (mol ->bonded_to (at, -2,7)+mol ->bonded_to (at, -2,8)) return POLAR;
    
    if (at->GetAtomicNum () == 6) {
        if (mol ->bonded_to (at, 3, N)) return POLAR;
        else if (mol ->bonded_to (at, 2, O)) return POLAR;
        else if ((mol ->bonded_to (at, -2,7)+mol ->bonded_to (at, -2,8))>2) return POLAR;
        else return NONPOLAR;
    } 
    if (at->GetAtomicNum () == 26) return METAL;  //iron other metals are missing
    return NONPOLAR;

*/
}






void Chemscore::clear_nonbonded_interactions () {
    for (unsigned int i=0; i<HBInteractions.size (); i++) delete HBInteractions[i];

    for (unsigned int i=0; i<LiInteractions.size (); i++) delete LiInteractions[i];

    for (unsigned int i=0; i<ClInteractions.size (); i++) delete ClInteractions[i];

    HBInteractions.clear ();
    ClInteractions.clear ();
    LiInteractions.clear ();
}





void Chemscore::compute_forces () {
    for (unsigned int i=0; i<HBInteractions.size (); i++) {
        HBInteractions[i] -> set_forces ();
    }
    for (unsigned int i=0; i<LiInteractions.size (); i++) {
        LiInteractions[i] -> set_forces ();
    }
    for (unsigned int i=0; i<ClInteractions.size (); i++) {
        ClInteractions[i] -> set_forces ();
    }
}


/*void Chemscore::compute_force (ChemscoreHBInteraction *cint) {
    vect dir = subtract (cint -> at1 -> GetVector (), cint -> at2 -> GetVector ());
    float force_x = cint -> derive (cint -> at1 -> GetVector ().x);
    dir.multiply (force_x / dir.x);
    cin -> at2 -> force = sum (cin -> at2 -> force, dir);
    dir.multiply (-1.f);
    cin -> at1 -> force = sum (cin -> at2 -> force, dir);
}
*/

/*
void Chemscore::compute_force (ChemscoreHBInteraction *cint) {
    vect dir = subtract (cint -> at1 -> GetVector (), cint -> at2 -> GetVector ())
    cint->at1-> GetVector ().x -= DX;
    float E = compute_HB_interaction (cint);
    cint->at1-> GetVector ().x += 2*DX;
    float E2 = compute_HB_interaction (cint);
    cint->at1-> GetVector ().x -= DX;
    float force = - 3* (E2-E)/(2*DX); //3* force
    float x = dir.x;
    vect force_vec = dir;
    force_vec.scale_to (force/x);
    cint->at1->force = force_vec;
    force_vec.multiply (-1.f);
    cint->at2->force = force_vec;
    }
}


void Chemscore::compute_force (ChemscoreLiInteraction *cint) {
    for (unsigned int i=0; i<3; i++) {
        cint->at1-> GetVector ()[i]-= DX;
        float E = compute_Li_interaction (cint);
        cint->at1-> GetVector ()[i]+= 2*DX;
        float E2 = compute_Li_interaction (cint);
        cint->at1-> GetVector ()[i]-= DX;
        cint->at1->force[i]-= 3* (E2-E)/(2*DX); //3* force
        cint->at2->force[i]+= 3* (E2-E)/(2*DX);
    }
}

void Chemscore::compute_force (ChemscoreClInteraction *cint) {
    for (unsigned int i=0; i<3; i++) {
        cint->at1-> GetVector ()[i]-= DX;
        float E = compute_Cl_interaction (cint);
        cint->at1-> GetVector ()[i]+= 2*DX;
        float E2 = compute_Cl_interaction (cint);
        cint->at1-> GetVector ()[i]-= DX;
        cint->at1->force[i]-= 3* (E2-E)/(2*DX); //3* force
        cint->at2->force[i]+= 3* (E2-E)/(2*DX);
    }
}


float Chemscore::compute_HB_interaction (ChemscoreHBInteraction *hbint) {
    float r = distance (hbint->at1-> GetVector (), hbint->at2-> GetVector ());
    float ang;
    if (hbint->at1->GetAtomicNum () ==1) ang = angle (hbint->at2-> GetVector (), hbint->at1-> GetVector (), hbint->root-> GetVector ());
    else ang = angle (hbint->at1-> GetVector (), hbint->at2-> GetVector (), hbint->root-> GetVector ());
    float dghbond = -3.34;
    float r0 = 1.85;
    float a0 = 180;
    float drab = r-r0;
    if (drab<0) drab = -drab;
    float daab = ang-a0;
    if (daab<0) daab = -daab;
    float dr1 = 0.25;
    float dr2 = 0.65;
    float da1 = 30;
    float da2 = 80;

    return dghbond*f(drab, dr1, dr2)*f(daab, da1, da2);
}

float Chemscore::compute_Li_interaction (ChemscoreLiInteraction *liint) {
    float r1, r2, dg;
    if (liint->metal) {
        dg = -6.03;
        r1 = 2.6;
        r2 = 3.0;
    }
    else {
        dg = -0.117;
        r1 = 4.1;
        r2 = 7.1;
    }
    float r = distance (liint->at1-> GetVector (), liint->at2-> GetVector ());

    return dg*f(r, r1, r2);
}

float Chemscore::compute_Cl_interaction (ChemscoreClInteraction *clint) {
    float r = distance (clint->at1-> GetVector (), clint->at2-> GetVector ());
    if (clint->type ==0) {
        float rcl = 1.6;
        float dghb = -3.34;
        if (r>rcl) return 0;
        else return (20/dghb)*(rcl-r)/rcl;
    }
    else if (clint->type ==1) {
        float rcl = 1.38;
        float dgm = -6.03;
        if (r>rcl) return 0;
        else return (20/dgm)*(rcl-r)/rcl;

    }
    else {
        float rcl;
        if (clint->at1->GetAtomicNum () == 16 || clint->at2->GetAtomicNum () == 16) rcl=3.35;
        else rcl = 3.1; 
        if (r>rcl) return 0;
        else return 1+4*(rcl-r)/rcl-1;
    }
}


*/




void Chemscore::load_internal_interactions () {}

void Chemscore::load_nonbonded_interactions () {

    clear_nonbonded_interactions ();
    Molecule *mol = target_mol;

    FOR_ATOMS_OF_MOL(a, mol)
    {
        objectList<Atom*>* nbAtoms = far_grid->getNeighborObjects(a -> GetVector ());
        if (nbAtoms) {
            vector <Atom *> neighbours = nbAtoms->objects;      
            for (unsigned int j=0; j<neighbours.size (); j++) {
                Atom * root = NULL;
                bool hb = false;
                bool li = false;
                bool me = false;
                int ltype = getChemscoretype (&*a);
                int ptype = getChemscoretype (neighbours[j]);
                if (ltype == DONOR) {
                    if (ptype == ACCEPTOR) {
                        hb = true;
                        OBBondIterator i;
                        root = a -> BeginNbrAtom (i);
                    }
                }
                else if (ltype == ACCEPTOR) {
                    if (ptype == DONOR) {
                        hb = true;
                        OBBondIterator i;
                        root = neighbours[j] -> BeginNbrAtom (i);
                    }
                    else if (ptype == METAL) {
                        me = true;
                    }
                }
                else if (ltype == NONPOLAR) {

                    if (ptype == NONPOLAR) {
                        li = true;
                    }
                }
                if (hb) {
                    ChemscoreHBInteraction *hbint = new ChemscoreHBInteraction;
                    hbint->at1 = &*a;
                    hbint->at2=neighbours[j];
                    hbint->root = root;
                    HBInteractions.push_back (hbint);  
                }    
                else if (li || me) {
                    ChemscoreLiInteraction *liint = new ChemscoreLiInteraction;
                    liint->at1 = &*a;
                    liint->at2 = neighbours[j];
                    liint->metal = me;

                    LiInteractions.push_back (liint);  
                }
                //clash
                if (neighbours[j] -> GetAtomicNum ()==1 || a -> GetAtomicNum () ==1) continue;
                ChemscoreClInteraction *clint = new ChemscoreClInteraction;
                clint->at1 = &*a;
                clint->at2=neighbours[j];
                if (ptype == METAL && ltype == ACCEPTOR) clint->type = 1;
                else if ((clint->at1->GetAtomicNum () ==7 || clint->at1->GetAtomicNum () ==8) && (clint->at2->GetAtomicNum () ==7 || clint->at2->GetAtomicNum () ==8)) {
                    if (mol -> bonded_to (clint->at1,1, 1) || mol -> bonded_to (clint->at2,1,1)) clint->type = 0;
                }
                else clint->type = 2;
                ClInteractions.push_back (clint);  
                


            }        
        }
    }
}


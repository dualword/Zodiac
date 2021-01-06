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



using namespace OpenBabel;

float f (float x, float x1, float x2) {
    if (x<= x1) return 1;
    else if (x>x2) return 0;
    else return (x2-x)/(x2-x1);    

}


float ChemscoreHBInteraction::value () {
    float r = dist (get_coordinates (at1), get_coordinates (at2));
    float ang;
    if (at1->GetAtomicNum () == 1) ang = angle (get_coordinates (at2), get_coordinates (at1), get_coordinates (root));
    else ang = angle (get_coordinates (at1), get_coordinates (at2), get_coordinates (root));
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
//		if (isnan (dghbond * f(drab, dr1, dr2) * f(daab, da1, da2))) cerr << "isnan chemscore" << drab << " "<< dr1 << " " << f(daab, da1, da2) << endl;
//	cerr << "value" <<dghbond * f(drab, dr1, dr2) * f(daab, da1, da2)<< endl;
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
    float r = dist (get_coordinates (at1), get_coordinates (at2));
	//if (isnan (dg * f(r, r1, r2))) cerr << "isnan chemscore" << dg << " "<< r << " " << f(r, r1, r2) << endl;

    return dg * f(r, r1, r2);

}




float ChemscoreClInteraction::value () {
    float r = dist (get_coordinates (at1), get_coordinates (at2));
    if (type == 0) {
        float rcl = 1.6f;
        float dghb = -3.34f;
        if (r>rcl) return 0.f;
        else return (20.f/dghb)*(rcl-r)/rcl;
    }
    else if (type == 1) {
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




Chemscore::Chemscore () : ForceField () {
    is_initialised = TRUE;
	LiGrid = NULL;
	ClGrid = NULL;
}

int Chemscore::getChemscoretype (Atom *at) {
    ZNMolecule *mol = (ZNMolecule *) at -> GetParent ();

    if (at->GetAtomicNum () == 8) { 
        return ZN_CS_ACCEPTOR;
    }
    if (at->GetAtomicNum () == 7) {
        if (CountBonds (at) <3) return ZN_CS_ACCEPTOR;
        else return ZN_CS_POLAR;
    }
    if (at->GetAtomicNum () == 1 && (mol -> bonded_to (at, -2, 7) || mol -> bonded_to (at, -2, 8))) return ZN_CS_DONOR;
    if (at->GetAtomicNum () == 15) return ZN_CS_POLAR;
    if (at->GetAtomicNum () == 9) return ZN_CS_POLAR;
    if (at->GetAtomicNum () == 16 && (mol->bonded_to (at, -2,7)+mol ->bonded_to (at, -2,8))) return ZN_CS_POLAR;
    
    if (at->GetAtomicNum () == 6) {
        if (mol ->bonded_to (at, 3, 7)) return ZN_CS_POLAR;
        else if (mol ->bonded_to (at, 2, 8)) return ZN_CS_POLAR;
        else if ((mol ->bonded_to (at, -2,7)+mol ->bonded_to (at, -2,8))>2) return ZN_CS_POLAR;
        else return ZN_CS_NONPOLAR;
    } 
    if (at->GetAtomicNum () == 26) return ZN_CS_METAL;  //iron other metals are missing
    return ZN_CS_NONPOLAR;








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
	HBInteractions.clear ();
    for (unsigned int i=0; i<ClInteractions.size (); i++) delete ClInteractions[i];
	ClInteractions.clear ();
    for (unsigned int i=0; i<LiInteractions.size (); i++) delete LiInteractions[i];
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


void Chemscore::load_grids (vect cent, double rad) {
	vect min_corner = vect (cent.x()-rad,cent.y()-rad,cent.z()-rad); 
	vect max_corner = vect (cent.x()+rad,cent.y()+rad,cent.z()+rad); 
	LiGrid = new DataGrid (min_corner, max_corner, 0.3);
	ClGrid = new DataGrid (min_corner, max_corner, 0.3);	
	for (int k = 0; k < LiGrid -> z_size (); k++) {
		for (int j = 0; j < LiGrid -> y_size (); j++) {
			for (int i = 0; i < LiGrid -> x_size (); i++) {
			float xx = LiGrid ->get_x (i);
			float yy = LiGrid ->get_y (j);
			float zz = LiGrid ->get_z (k);
				LiGrid ->set_value (i, j, k, Livalue (vect (xx, yy, zz)));
				ClGrid ->set_value (i, j, k, Clvalue (vect (xx, yy, zz)));
			}
		}
	}
	
	
}




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

#include "database.h"
#include "obabel_includes.h"
#include "FF.h"
#include "constants.h"
#include "cutoffGrid.h"

#ifdef WIN32
#include <float.h>
#define isnan _isnan
#endif // WIN32



ForceFieldInteraction::ForceFieldInteraction () {
    at1 = NULL;
    at2 = NULL;
}

/*
double ForceFieldInteraction::derive (double &variable, double *function_value) {
    double ref_value = variable;
    variable = ref_value - DX;
	cerr << "var-dx"<<variable << endl;
    double E = value ();
    variable = ref_value + DX;
	cerr << "var+dx"<<variable<<endl;
    double E2 = value ();
    variable = ref_value;
    if (function_value) *function_value = (E2 + E) / 2.f;
    assert (!isnan ((E2-E)/(2*DX)));
//	cerr << "derive" << (E2-E)/(2*DX)<<endl; 
    return  (E2-E)/(2*DX);
}*/

double ForceFieldInteraction::derive (Atom *at, double *function_value) {
    vect v = get_coordinates(at);
	double ref = v.x ();
	set_coordinates (at,vect (ref-DX, v.y (), v.z ()));
    double E = value ();
	set_coordinates (at,vect (ref+DX, v.y (), v.z ()));
    double E2 = value ();
	set_coordinates (at,vect (ref, v.y (), v.z ()));
    if (function_value) *function_value = (E2 + E) / 2.f;	
    assert (!isnan ((E2-E)/(2*DX)));
//	cerr << "derive" << (E2-E)/(2*DX)<<endl; 
    return  (E2-E)/(2*DX);
}

double ForceFieldInteraction::derive_y (Atom *at, double *function_value) {
    vect v = get_coordinates(at);
	double ref = v.y ();
	set_coordinates (at, vect (v.x (), ref-DX, v.z ()));
    double E = value ();
	set_coordinates (at, vect (v.x (), ref+DX, v.z ()));
    double E2 = value ();
	set_coordinates (at, vect (v.x (), ref, v.z ()));
    if (function_value) *function_value = (E2 + E) / 2.f;
    assert (!isnan ((E2-E)/(2*DX)));
//	cerr << "derive" << (E2-E)/(2*DX)<<endl; 
    return  (E2-E)/(2*DX);
}

double ForceFieldInteraction::derive_z (Atom *at, double *function_value) {
    vect v = get_coordinates(at);
	double ref = v.z ();
	set_coordinates (at, vect (v.x (), v.y (), ref-DX));
    double E = value ();
	set_coordinates (at, vect (v.x (), v.y (), ref+DX));
    double E2 = value ();
	set_coordinates (at, vect (v.x (), v.y (), ref));
    if (function_value) *function_value = (E2 + E) / 2.f;
    assert (!isnan ((E2-E)/(2*DX)));
//	cerr << "derive" << (E2-E)/(2*DX)<<endl; 
    return  (E2-E)/(2*DX);
}



void ForceFieldInteraction::set_forces (bool score, double mult) {
    assert (at1);
    assert (at2);
	vect v1 = get_coordinates (at1);
	double x1 = v1.x();
	//cerr << "v1" << v1 << endl;
	//cerr << "x1" << x1 << endl;
    assert (at1 != at2);
    vect dir = subtract (get_coordinates(at1), get_coordinates(at2));
	if (!dir.module ()) dir = vect (1., 0., 0.);
    assert (dir.module ());
    double score_val;
    double force_x;
    force_x = derive (at1, &score_val);
    dir.multiply (force_x / dir.x());
	dir.multiply(mult);
    vect force1 = get_back_force (at1);
    vect force2 = get_back_force (at2);
    force2 += dir;
    force1 -= dir;

    set_back_force (at1, force1);
	set_back_force (at2, force2);


    if (score) {
		double score1 = get_back_score (at1);
		double score2 = get_back_score (at2);
	//	cerr <<score1;
		score1 += score_val;
		score2 += score_val;
	//	cerr << score1 << score_val<< endl;

		set_back_score (at1, score1);
		set_back_score (at2, score2);

    }

}

ElasticRestrain::ElasticRestrain () : ForceFieldInteraction () { k = 1.f; dist0 = 1.5;}

float ElasticRestrain::value () {
	float d = dist (get_coordinates(at1), get_coordinates(at2)) - dist0;
	return  k * d * d;
}



ForceField::ForceField () {
    is_initialised = false;
    target_mol = NULL;
    near_grid = NULL;
    far_grid = NULL;
}

ForceField::~ForceField() {}

void ForceField::clear_nonbonded_interactions () {
// cerr <<"general clearnonb"<<endl;
}

void ForceField::clear () {
    clear_nonbonded_interactions ();
    environment.clear ();
    target_mol = NULL;
    delete near_grid;
    delete far_grid;
    near_grid = NULL;
    far_grid = NULL;
}


void ForceField::load_mol (ZNMolecule *mol) {
    assert (mol);
    target_mol = mol;
}


void ForceField::load_environment (vector<ZNMolecule *> envir, ZNMolecule *mol, vect cent, double rad) {
    environment.clear ();
    for (unsigned int i=0; i<envir.size (); i++) {
        FOR_ATOMS_OF_MOL (a, envir[i]) {
			ZNMolecule *other_mol = (ZNMolecule *) a->GetParent ();
			if (other_mol == mol) continue;
			if (other_mol ->multi && mol ->multi) {
				
				Database_molecule *dbm1 = (Database_molecule *) mol;
				Database_molecule *dbm2 = (Database_molecule *) other_mol;
				if (dbm1 -> database == dbm2 -> database) continue;				
			}
			environment.push_back (&*a); 
		}
	}
    delete near_grid;
    delete far_grid;
    near_grid  = new cutoffGrid<Atom*> (environment, FF_NEAR_NONBONDED_CUTOFF);
    far_grid   = new cutoffGrid<Atom*> (environment, FF_FAR_NONBONDED_CUTOFF );
}


void ForceField::update () {
   load_nonbonded_interactions (); 
}

void ForceField::initialize_internal (ZNMolecule* mol, vector<ZNMolecule *> envir) {
    clear ();
    load_environment (envir, mol);
    load_mol (mol);
    load_internal_interactions (0);
}

void ForceField::initialize_interaction (ZNMolecule* mol, vector<ZNMolecule *> envir, vect cent, float rad) {
    clear ();
    load_environment (envir, mol);
    load_mol (mol);
  //  load_nonbonded_interactions ();
	load_grids (cent, rad);


}

void ForceField::initialize (ZNMolecule* mol, vector<ZNMolecule *> envir) {
    clear ();
    load_environment (envir, mol);
    load_mol (mol);
    load_internal_interactions (0);
    load_nonbonded_interactions ();

}



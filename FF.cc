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
#include <openbabel/forcefield.h>
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
    vect &v = (vect &) at -> GetVector ();
	double ref = v.x ();
	at -> SetVector (vect (ref-DX, v.y (), v.z ()));
    double E = value ();
	at -> SetVector (vect (ref+DX, v.y (), v.z ()));
    double E2 = value ();
	at -> SetVector (vect (ref, v.y (), v.z ()));
    if (function_value) *function_value = (E2 + E) / 2.f;	
    assert (!isnan ((E2-E)/(2*DX)));
//	cerr << "derive" << (E2-E)/(2*DX)<<endl; 
    return  (E2-E)/(2*DX);
}

double ForceFieldInteraction::derive_y (Atom *at, double *function_value) {
    vect &v = (vect &) at -> GetVector ();
	double ref = v.y ();
	at -> SetVector (vect (v.x (), ref-DX, v.z ()));
    double E = value ();
	at -> SetVector (vect (v.x (), ref+DX, v.z ()));
    double E2 = value ();
	at -> SetVector (vect (v.x (), ref, v.z ()));
    if (function_value) *function_value = (E2 + E) / 2.f;
    assert (!isnan ((E2-E)/(2*DX)));
//	cerr << "derive" << (E2-E)/(2*DX)<<endl; 
    return  (E2-E)/(2*DX);
}

double ForceFieldInteraction::derive_z (Atom *at, double *function_value) {
    vect &v = (vect &) at -> GetVector ();
	double ref = v.z ();
	at -> SetVector (vect (v.x (), v.y (), ref-DX));
    double E = value ();
	at -> SetVector (vect (v.x (), v.y (), ref+DX));
    double E2 = value ();
	at -> SetVector (vect (v.x (), v.y (), ref));
    if (function_value) *function_value = (E2 + E) / 2.f;
    assert (!isnan ((E2-E)/(2*DX)));
//	cerr << "derive" << (E2-E)/(2*DX)<<endl; 
    return  (E2-E)/(2*DX);
}



void ForceFieldInteraction::set_forces (bool score) {
    assert (at1);
    assert (at2);
	vect &v1 = (vect &) at1 -> GetVector ();
	double &x1 = v1.x();
	//cerr << "v1" << v1 << endl;
	//cerr << "x1" << x1 << endl;
    assert (at1 != at2);
    vect dir = subtract ((vect&) at1 -> GetVector (), (vect&) at2 -> GetVector ());
	if (!dir.module ()) dir = vect (1., 0., 0.);
    assert (dir.module ());
    double score_val;
    double force_x;
    force_x = derive (at1, &score_val);
    dir.multiply (force_x / dir.x());
    vect force1 = get_force (at1);
    vect force2 = get_force (at2);
    force2 += dir;
    force1 -= dir;
    set_force (at1, force1);
	set_force (at2, force2);

    if (score) {
		double score1 = get_old_score (at1);
		double score2 = get_old_score (at2);
	//	cerr <<score1;
		score1 += score_val;
		score2 += score_val;
	//	cerr << score1 << score_val<< endl;
		set_old_score (at1, score1);
		set_old_score (at2, score2);
    }

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


void ForceField::load_mol (Molecule *mol) {
    assert (mol);
    target_mol = mol;
}


void ForceField::load_environment (vector<Molecule *> envir, Molecule *mol) {
    environment.clear ();
    for (unsigned int i=0; i<envir.size (); i++) {
        FOR_ATOMS_OF_MOL (a, envir[i]) {
			Molecule *other_mol = (Molecule *) a->GetParent ();
			if (other_mol == mol) continue;
			if (other_mol ->multi && mol ->multi) {
				
				Database_molecule *dbm1 = (Database_molecule *) mol;
				Database_molecule *dbm2 = (Database_molecule *) other_mol;
				if (dbm1 -> database == dbm2 -> database) continue;				
			}
			environment.push_back (&*a); 
		}
	}
//cout <<"envsize"<<environment.size ()<<endl;
    delete near_grid;
    delete far_grid;
    near_grid  = new cutoffGrid<Atom*> (environment, FF_NEAR_NONBONDED_CUTOFF);
    far_grid   = new cutoffGrid<Atom*> (environment, FF_FAR_NONBONDED_CUTOFF );
}


void ForceField::update () {
   load_nonbonded_interactions (); 
}

void ForceField::initialize_internal (Molecule* mol, vector<Molecule *> envir) {
    clear ();
    load_environment (envir, mol);
    load_mol (mol);
    load_internal_interactions ();
}

void ForceField::initialize_interaction (Molecule* mol, vector<Molecule *> envir) {
    clear ();
    load_environment (envir, mol);
    load_mol (mol);
    load_nonbonded_interactions ();

}

void ForceField::initialize (Molecule* mol, vector<Molecule *> envir) {
    clear ();
    load_environment (envir, mol);
    load_mol (mol);
    load_internal_interactions ();
    load_nonbonded_interactions ();

}



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
#ifndef MINIMIZE_H
#define MINIMIZE_H

#include "constants.h"
#include "ZNdata.h"
#include "FF.h"
#include "thread.h"



class HapticThread;

typedef struct {
	vect coordinates;
    Atom *atom;
} Fragment_Atom;

typedef struct {
    vector <Fragment_Atom> atoms;
    vect translation;
  //  float rotation [3];
    float rotation_quat [4];
} pFragment;

typedef struct {
	unsigned int type;
    float *value;
    float der;
} Dof;

class Data; 
class Minimize  {

public:
float counter;
float step;

Molecule *minimising_molecule;
Molecule *haptic_molecule;

ForceField *internal_ff;
ForceField *interaction_ff;
Minimize (Data *dat);


bool automove;
bool color_by_score;

float total_E;
float total_interaction_E;
float total_internal_E;

void start_haptic_mode ();
void stop_haptic_mode () ;
//void haptic_step ();
void initialize (Molecule *mol);
void initialize_6 (Molecule *mol);

void minimize_energy_step ();
float compute_energy ();
void apply_forces (Molecule *mol, float trunc = 0.f);
void apply_force_to_atom (Atom *a, float trunc = 0.f);


//void initialise_minimisation ();
void deinitialise_minimisation ();


void clear ();
void clear_fragments ();
int haptic_dof_mode;
int haptic_number_of_threads;
vector <pFragment> fragments;
HapticThread *haptic_thread;

void update_fragment_position (pFragment&);
private:
int forcefields_sanity_check ();

int iterations;
float last_E;
Data *data;


#ifdef HAPTICS
public: void update_molecule_position_with_haptic_pointer (Molecule *min_mol);
#endif //HAPTICS

};

#endif

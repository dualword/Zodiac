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
#include "ils.h"
#include "pso.h"
#include "nms.h"


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
	Minimize (Data *dat, ForceField *int_ff = 0, ForceField *inter_ff = 0);
	
	void set_molecule (ZNMolecule *mol) {target_molecule = mol;}
	
	
	ZNMolecule *target_molecule;
	Optimiser *optimiser;
	
	ZNMolecule *minimising_molecule;
	ZNMolecule *haptic_molecule;
	
	ForceField *internal_ff;
	vector <ForceField *> interaction_ffs;
	float counter;
	float step;
	
	bool automove;
	bool color_by_score;
	
	float total_E;
	float total_interaction_E;
	float total_internal_E;
	
	
	float score ();
	void start_haptic_mode ();

	void initialize (ZNMolecule *mol);
	void initialize_6 (ZNMolecule *mol);
	
	void minimize_energy_step ();
	float compute_energy ();
	void apply_forces (ZNMolecule *mol, float trunc = 0.f);
	void apply_force_to_atom (Atom *a, float trunc = 0.f);
	
	
	//void initialise_minimisation ();
	void deinitialise_minimisation ();
	
	
	void clear ();
	void clear_fragments ();
	int haptic_dof_mode;
	int haptic_number_of_threads;
	vector <pFragment> fragments;

	
	void update_fragment_position (pFragment&);
private:
	int forcefields_sanity_check ();
	
	int iterations;
	float last_E;
	Data *data;
	
	
//#ifdef HAPTICS
public: void update_molecule_position_with_haptic_pointer (ZNMolecule *min_mol);
//#endif //HAPTICS
	
};


class ScoreMolecule : public Function {
public:
	ScoreMolecule (Minimize *min, bool move_bool = true, vect bs_center = vect (), float bs_rad = 10) : Function (), minimise (min), moving (move_bool)
	{

		float r = bs_rad;
		int gap = 0;
		ZNMolecule *mol = min -> target_molecule;
		vect curr_c = get_root_fragment_center (mol);
		//cerr << curr_c<<bs_center<<r<< endl;
		vector <double> vv = get_dihedrals (mol);
		values.resize (gap + vv.size ());
		if (moving) {
			gap = 6; 
			values.resize (gap + vv.size ());
			values[0] = bs_center.x() - curr_c.x();
			Variable* v = new Variable;
			v->value = &values[0];
			*v->value = values[0];
			v ->min_val = bs_center.x() - curr_c.x() -r;
			v ->max_val = bs_center.x() - curr_c.x() +r;
			_variables.push_back(v);
			
			values[1] = bs_center.y() - curr_c.y();
			Variable* v2 = new Variable;
			v2->value = &values[1];
			*v2->value = values[1];
			v2 ->min_val = bs_center.y() - curr_c.y() -r;
			v2 ->max_val = bs_center.y() - curr_c.y() +r;
			_variables.push_back(v2);
			
			values[2] = bs_center.z() - curr_c.z();
			Variable* v3 = new Variable;
			v3->value = &values[2];
			*v3->value = values[2];
			v3 ->min_val = bs_center.z() - curr_c.z() -r;
			v3 ->max_val = bs_center.z() - curr_c.z() +r;
			_variables.push_back(v3);
			
			//rotation
			values[3] = 0.f;
			Variable* v4 = new Variable;
			v4->value = &values[3];
			*v4->value = values[3];
			v4 ->min_val = -M_PI;
			v4 ->max_val = M_PI;
			_variables.push_back(v4);
			
			values[4] = 0.f;
			Variable* v5 = new Variable;
			v5->value = &values[4];
			*v5->value = values[4];
			v5 ->min_val = -M_PI;
			v5 ->max_val = M_PI;
			_variables.push_back(v5);
			
			values[5] = 0.f;
			Variable* v6 = new Variable;
			v6->value = &values[5];
			*v6->value = values[5];
			v6 ->min_val = -M_PI;
			v6 ->max_val = M_PI;
			_variables.push_back(v6);
		}
		
		
		
		for (unsigned int i = 0; i < vv.size(); i++) {
			
			values[i+gap] = (float) vv[i];
			Variable* v = new Variable;
			v->value = &values[i+gap];
			*v->value = values[i+gap];
			v ->min_val = -M_PI;
			v ->max_val = M_PI;
			_variables.push_back(v);
		}
		
	};
	
	
	float evaluate ()
	{
		
		build_molecule_from_dofs (minimise ->target_molecule, values, moving);
		set_needs_redraw (minimise ->target_molecule, true);
	//	cerr << "score " << minimise ->score ();
		return minimise -> score ();
	};
	
private:
	vector<float> values;
	Minimize *minimise;
	bool moving;
};










/*

class LinkerFunction : public Function {
public:
	LinkerFunction (Minimize *min, Database *frags, int ) : Function (), minimise (min), _fragments_source (frags)
	{
		int max_fragments = 4;
		float r = 2.;

		values.resize (6 + max_fragments * 3);
		values[0] = 0.f;
		Variable* v = new Variable;
		v->value = &values[0];
		*v->value = values[0];
		v ->min_val = -r;
		v ->max_val = r;
		_variables.push_back(v);
			
		values[1] = 0.f;
		Variable* v2 = new Variable;
		v2->value = &values[1];
		*v2->value = values[1];
		v2 ->min_val = -r;
		v2 ->max_val = r;
		_variables.push_back(v2);
		
		values[2] = 0.f;
		Variable* v3 = new Variable;
		v3->value = &values[2];
		*v3->value = values[2];
		v3 ->min_val = -r;
		v3 ->max_val = r;
		_variables.push_back(v3);
			
		//rotation
		Variable* v4 = new Variable;
		v4->value = &values[3];
		*v4->value = values[3];
		v4 ->min_val = -M_PI;
		v4 ->max_val = M_PI;
		_variables.push_back(v4);
		
		values[1] = 0.f;
		Variable* v5 = new Variable;
		v5->value = &values[4];
		*v5->value = values[4];
		v5 ->min_val = -M_PI;
		v5 ->max_val = M_PI;
		_variables.push_back(v5);
		
		values[2] = 0.f;
		Variable* v6 = new Variable;
		v6->value = &values[5];
		*v6->value = values[5];
		v6 ->min_val = -M_PI;
		v6 ->max_val = M_PI;
		_variables.push_back(v6);
		
		
		for (unsigned int i = 0; i < max_fragments; i++) {
			
			values[i+6] = 0;
			Variable* v = new Variable;
			v->value = &values[i+6];
			*v->value = values[i+6];
			v ->min_val = 0;
			v ->max_val = _fragments_source -> count_molecules ();
			_variables.push_back(v);


			
			values[i+6+1] = 0;
			Variable* rot = new Variable;
			rot->value = &values[i+6+1];
			*rot->value = values[i+6+1];
			rot ->min_val = 0;
			rot ->max_val =  2 * M_PI;
			_variables.push_back(rot);
		}
		
	};
	
	
	float evaluate ()
	{
		
		arrange_fragments (fragments, values);
	//	set_needs_redraw (minimise ->target_molecule, true);
		float score = 0;
//		for (int i = 0; i < fragments.size ()) 
		return minimise -> score ();
	};
	
private:
	vector<float> values;
	Minimize *minimise;
	vector <ZNMolecule *> fragments;
	Database *_fragments_source;
	
	void arrange_fragments (vector <ZNMolecule *> frags, vector <float>) {};
};


*/


#endif

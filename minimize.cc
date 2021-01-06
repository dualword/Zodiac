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
#include "minimize.h"

#ifdef HAPTICS
#include "haptics.h"
#endif //HAPTICS

#ifdef WIN32
#include <float.h>
#define isnan _isnan
#endif // WIN32

Minimize::Minimize (Data *dat) {
	data = dat;
	internal_ff = new MMFF ();
	interaction_ff = new MMFF ();
	haptic_dof_mode = 2;
	total_E = 0.f;
	total_internal_E = 0.f;
	total_interaction_E = 0.f;
	clear ();
    minimising_molecule = NULL;
    haptic_molecule = NULL;
 //   minimise_reference_molecule = NULL;


}



int Minimize::forcefields_sanity_check () {
	string error_string = "";
	int out = 1;
	if (!internal_ff -> is_initialised) {
		error_string += "Internal Forcefield not correctly initialised.\n";
		out = 0;
	}
	if (!interaction_ff -> is_initialised) {
		error_string += "Interaction Forcefield not correctly initialised.\n";
		out = 0;
	}
	if (!out)     QMessageBox::critical(data->ddwin, "Forcefields sanity check failed" , QString (error_string.c_str ()) );
	return out;
}


/*
void Minimize::initialise_minimisation (Molecule *mol, Forcefield *ff) {
    clear ();

    ff -> clear_internal_interactions ();
    ff -> clear_nonbonded_interactions (); 
    minimising_molecule = data -> ddwin -> target_molecule;  


    data->mmff->initialize (mol, data -> ddwin -> molecules);  



    data -> undo_stack -> beginMacro ("Energy Minimisation");
    MoveAtomsCommand *command = new MoveAtomsCommand (data -> ddwin -> gl, 1);
    FOR_ATOMS_OF_MOL(a, minimising_molecule) {
        command -> add (&*a, (vect &) a -> GetVector ());
    }

    data -> ddwin -> execute (command);
	


}
*/

void Minimize::deinitialise_minimisation () {

	clear ();
    minimising_molecule = NULL;

}


void Minimize::start_haptic_mode () {

//	if (forcefields_sanity_check ()) {
		initialize (data->ddwin->target_molecule);
		data->ddwin->haptic_mode = true;
		data -> ddwin -> lock_editing ();
		total_E = 0.;
		total_internal_E = 0.;
		total_interaction_E = 0.;
		if (haptic_dof_mode!=0) {
			internal_ff->initialize_internal (data->ddwin->target_molecule, data->ddwin->molecules);
		}
		interaction_ff->initialize_interaction (data->ddwin->target_molecule, data->ddwin->molecules);
        data -> undo_stack -> beginMacro ("Haptic Simulation");
        MoveAtomsCommand *command = new MoveAtomsCommand (data -> ddwin -> gl, 1);
        FOR_ATOMS_OF_MOL(a, haptic_molecule) {
            command -> add (&*a, (vect &) a -> GetVector ());
        }
//        command -> name ("haptic simulation");
        data -> ddwin -> execute (command);
//	}
        haptic_thread -> initialise (this);
        haptic_thread -> start ();

}

void Minimize::stop_haptic_mode () {
	if (haptic_thread-> isRunning ()) {
		haptic_thread -> stop ();
		haptic_thread -> wait ();
		MoveAtomsCommand *command = new MoveAtomsCommand (data -> ddwin -> gl, 0);
		FOR_ATOMS_OF_MOL(a, haptic_molecule) {
			command -> add (&*a, a -> GetVector ());
		}
  //  command -> name ("haptic simulation");
		data -> ddwin -> execute (command);	
		data -> undo_stack -> endMacro ();
		data->ddwin->haptic_mode = false;
		haptic_molecule = NULL;
		data -> ddwin -> unlock_editing ();
	}
}

void Minimize::initialize (Molecule *mol) {
    FOR_ATOMS_OF_MOL(a, mol) {
  //      vect *force = (vect *) a -> GetData ("force");
 //       force -> null ();
	}

	if (haptic_dof_mode ==0) { 
		initialize_6 (mol);
	}
    haptic_molecule = mol;
}


void Minimize::initialize_6 (Molecule *mol) { 
	pFragment frag;
    frag.translation = mol -> center;
	frag.rotation_quat [0] = 1.f;
	frag.rotation_quat [1] = 0.f;
	frag.rotation_quat [2] = 0.f;
	frag.rotation_quat [3] = 0.f;

        FOR_ATOMS_OF_MOL(a, mol) {
		Fragment_Atom fa;
        fa.coordinates = subtract (a -> GetVector (), mol -> center);
		fa.atom = &*a;
		frag.atoms.push_back (fa);
	}
	fragments.push_back (frag);

}


void Minimize::clear () {
	iterations = 0;
	step = 1;
	last_E = 0;
	clear_fragments ();
}

void Minimize::clear_fragments () {
	fragments.clear ();
}
/*
void Minimize::haptic_step () {
	const float maxforce = 200.0f;

	counter ++;

	Molecule * min_mol = haptic_molecule;
    assert (min_mol);

#ifdef HAPTICS
	update_molecule_position_with_haptic_pointer (min_mol);
#endif //HAPTICS

    vect haptic_force (0.f, 0.f, 0.f);

    for (unsigned int i=0; i<min_mol->atoms.size (); i++) {
		min_mol->atoms[i]->score = 0.f;
        min_mol -> atoms [i] -> force.null ();
	}



	if (haptic_dof_mode ==0) { // 6 dofs model


        vect tot_force, tot_torque;
        tot_force.null ();
        tot_torque.null ();
		interaction_ff->update ();
		interaction_ff->compute_forces ();




		for (unsigned int fi=0; fi<fragments.size (); fi++) { 

			for (unsigned int i=0; i<fragments[fi].atoms.size (); i++) {
			    tot_force = sum (tot_force, fragments[fi].atoms[i].atom->force);
				vect to = torque (fragments[fi].atoms[i].atom->force, fragments[fi].atoms[i].atom-> GetVector (), fragments[fi].translation);
                tot_torque = sum (tot_torque, to);
			}

			haptic_force = sum (haptic_force, tot_force);
			if (automove) {
                vect force = tot_force;
                force.trunc_at (500.f);
                force.multiply (STEP_SIZE);
				fragments[fi].translation = sum (fragments[fi].translation, force); 
			}
		    float cut = 10.f;
            tot_torque.trunc_at (cut);

			float new_quat [4], mult_quat [4];
			axis_angle_to_quaternion (tot_torque, tot_torque.module ()*0.001, new_quat);

			multiply_quaternions (new_quat, fragments[fi].rotation_quat, mult_quat);
			normalize_quaternion (mult_quat);
			fragments[fi].rotation_quat[0] = mult_quat [0];
			fragments[fi].rotation_quat[1] = mult_quat [1];
			fragments[fi].rotation_quat[2] = mult_quat [2];
			fragments[fi].rotation_quat[3] = mult_quat [3];
			update_fragment_position (fragments[fi]);
		}

	}



	else if (haptic_dof_mode ==1) { // 6+r dofs model
	}
	else { //3N model
		vect lastCenter;
        lastCenter = min_mol -> center;
        assert (min_mol->atoms.size ());


		internal_ff->compute_forces ();
		interaction_ff->update ();
		interaction_ff->compute_forces ();
        apply_forces (min_mol, maxforce);

		for (unsigned int i=0; i<min_mol->atoms.size (); i++) {
            vect force = min_mol -> atoms[i] -> force;
    		haptic_force = sum (haptic_force, force);
		}




		if (!automove) {
			vect cent = find_mass_center (min_mol->atoms);
			for (unsigned int i=0; i<min_mol->atoms.size (); i++) {
	    		min_mol->atoms[i]-> GetVector () = subtract (min_mol->atoms[i]-> GetVector (), cent);
				min_mol->atoms[i]-> GetVector () = sum (min_mol->atoms[i]-> GetVector (), lastCenter);

			}
		}
    //	total_interaction_E = interaction_ff->total_energy;
    //	total_internal_E = internal_ff->total_energy;
	//    total_E = total_interaction_E + total_internal_E;

//	    data->ddwin->haptic_menu->update ();

        #ifdef HAPTICS


	    // HAPTIC FUNCTION needs some force cut off and scaling 
	    vect actual_force;
	    data->ddwin->gl->world_to_haptic_coordinates (haptic_force, actual_force, -5, 5, -5, 5, -5, 5);
	    gHaptics.setForce(
		haptic_force.x() / maxforce, /////////////why not actual force???
		haptic_force.y() / maxforce,
		haptic_force.z() / maxforce
		);
        #endif //HAPTICS

    }
    if (color_by_score) {
    //    cout << "color_bu"<<endl;
        vector <color_mask> masks;
        color_mask mask;
        mask.intensity = 1.0f;
        mask.only_to = 0;
        mask.excluding = 0;
        mask.type =  2; //score //see menu.cc
        masks.push_back (mask);
        data -> ddwin->gl->apply_color_masks (masks, min_mol, false);
    }
}
*/
#ifdef HAPTICS
void Minimize::update_molecule_position_with_haptic_pointer (Molecule *min_mol) {


	vect haptic_coords, new_coords;
	float x, y, z;
	double pos[3];
	gHaptics.getPosition(pos);
	x = pos[0]; y = pos[1]; z = pos[2];

	/*
	x = 0; 
	y = 0;
	z = 5 * sin (counter/100);
	*/
	haptic_coords.x() = x;
	haptic_coords.y() = y;
	haptic_coords.z() = z;

	float minx, maxx, miny, maxy, minz, maxz;

	minx = -5; maxx = 5; 
	miny = -5; maxy = 5;
	minz = -5; maxz = 5;

	data->ddwin->gl->haptic_to_world_coordinates (haptic_coords, new_coords, minx, maxx, miny, maxy, minz, maxz);

    FOR_ATOMS_OF_MOL(a, min_mol) {
        a -> GetVector () += new_coords;
        a -> GetVector () -= center;
	}

    min_mol -> center = new_coords;
}

#endif //HAPTICS


void Minimize::update_fragment_position (pFragment& frag) {
	for (unsigned int i=0; i<frag.atoms.size (); i++) {
		vect rotated_coords;
		rotated_coords = rotate_vector_using_quaternion (frag.atoms[i].coordinates, frag.rotation_quat);
		frag.atoms[i].atom-> GetVector () = sum (rotated_coords, frag.translation);
	}
}


void Minimize::apply_forces (Molecule *mol, float trunc) {
  //  assert (mol -> atoms.size ());
	FOR_ATOMS_OF_MOL(a, mol) {
        vect force = get_force (&*a);
	//	cerr << force << "force"<<endl;
        assert (!isnan (force.x()));
        assert (!isnan (force.y()));
        assert (!isnan (force.z()));
        if (trunc > 0.f) 		force.trunc_at (trunc);
        force.multiply (STEP_SIZE);
	//	cerr << force << endl;
        assert (!isnan (force.x()));
        assert (!isnan (force.y()));
        assert (!isnan (force.z()));
		vect coord = (vect &) a -> GetVector ();
		a -> SetVector (coord + force);
    }
}


void Minimize::minimize_energy_step () { 
	bool converged = false;
    FOR_ATOMS_OF_MOL(a, data->mmff->target_mol) {
        vect *force = (vect *) a -> GetData ("force");
		force -> null ();
	}
	iterations ++;
	data -> mmff -> update ();
	data -> mmff -> compute_forces ();
    apply_forces (data -> mmff -> target_mol, 500.f);

	float this_E = compute_energy ();

	if (last_E-this_E<MIN_ENERGY) converged=true;
	if (iterations ==1) {
		converged = false; 
	}
	if (last_E-this_E<0) {
		converged = false;
		step /=1.1f;
	}
	else step *=1.1f;
	last_E = this_E;
	data->ddwin->gl->draw_molecule (data->ddwin->target_molecule);


	if (converged || iterations > MAX_ITERATIONS) {
        deinitialise_minimisation ();
	}
} 



/*


void Minimize::minimize_energy () {
int iterations = 0;
float step =1;
bool converged = false;
float last_E, this_E;
last_E =0;
this_E =0;
//   cout <<compute_energy ()<<"energy"<<endl;
while (!converged && iterations<MAX_ITERATIONS) {
for (unsigned int i=0; i<data->mmff->target_mol->atoms.size (); i++) {
//           ddwin->data->mmff->target_mol->atoms[i]->score = 0.;
data->mmff->target_mol->atoms[i]->force[0] = 0.;
data->mmff->target_mol->atoms[i]->force[1] = 0.;
data->mmff->target_mol->atoms[i]->force[2] = 0.;
}

iterations ++;
data->mmff->update ();
data->mmff->compute_forces ();

for (unsigned int i=0; i<data->mmff->target_mol->atoms.size (); i++) {
for (unsigned int j=0; j<3; j++) {
float force = data->mmff->target_mol->atoms[i]->force[j];
if (force > 500) force = 500;
if (force < -500) force = -500;

data->mmff->target_mol->atoms[i]-> GetVector ()[j]+=force*STEP_SIZE;
}
}


this_E = compute_energy ();

if (last_E-this_E<MIN_ENERGY) converged=true;
if (iterations ==1) {
converged = false; 
}
if (last_E-this_E<0) {
converged = false;
step /=1.1;
}
else step *=1.1;
last_E = this_E;
data->ddwin->gl->draw_molecule (data->ddwin->target_molecule);
} 
}
*/
/*

float Minimize::derive (Dof *dof) {
*dof->value-=DX;
float E1 = compute_energy ();
*dof->value+=2*DX;
float E2 = compute_energy ();
*dof->value-=DX;    
return (E2-E1)/(2*DX);
}
*/

float Minimize::compute_energy () {
	//   cout <<data->mmff->compute_total_energy ();
	return data->mmff->compute_total_energy ();
}
